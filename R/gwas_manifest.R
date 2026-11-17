#' Download, review, and ingest GWAS studies from a manifest
#'
#' Each row represents one study. New inputs are downloaded or read from the
#' local cache and checked against the supplied MD5. Automatic import requires
#' a verified source and header match from [inspectGWAS()] and a passing
#' preview. Review the per-study inspection files for rows marked for review.
#' Files downloaded by gact are removed after a successful import and descriptor
#' save. Preview runs, failed imports, and inputs supplied directly to the
#' cache retain their files.
#'
#' @param manifest A TSV path or data frame with `accession`, `filename`,
#'   `url`, `md5`, `source`, `trait`, `type`, `build`, and `ancestry`.
#'   Use unique accessions and filenames, a basename with a supported file
#'   extension, and an HTTPS URL. URL can be blank for a cached file.
#'   Shared study metadata can be supplied through `defaults`. Other metadata
#'   fields from [inspectGWAS()] are optional.
#' @param descriptor Path to the existing database's GAlist RDS. Successful
#'   imports update this file. Previously imported studies are detected by
#'   accession, including legacy GWAS Catalog study URLs in `reference`.
#' @param data_dir Download and cache directory. Defaults to the database's
#'   download directory.
#' @param out_dir New or empty report directory outside `data_dir`. Use a
#'   fresh report directory on subsequent runs to preserve earlier reports.
#'   Defaults to a new batch directory under the database's recipes/reviews.
#' @param knowledge Path to verified `formats.json` from [inspectGWAS()].
#'   With no verified match, the workflow produces a preview for review.
#' @param defaults Named list of shared study metadata for this release.
#'   Nonempty manifest values take precedence. A GWAS Catalog accession
#'   supplies the study URL when `reference` is missing.
#' @param import Import automatically after a verified preview. Use `FALSE`
#'   to check readiness without writing to the database.
#' @param preview_rows Maximum number of rows per preview.
#' @param engine Select the R importer or the C++ batch engine for supported
#'   verified gzip TSV profiles.
#' @param workers Number of C++ staging workers on Unix. Windows uses one worker.
#' @param scratch_dir Parent directory for temporary C++ working files.
#'   Defaults to SLURM_TMPDIR when available, or R's temporary directory.
#'   Each run removes its own working directory when it exits.
#' @return Invisibly, the per-study report written to `manifest_review.tsv`.
#'   Status is `imported`, `already_imported`, `ready`, or `needs_review`.
#'   Preview QC warnings are reported separately. Full-data QC counts are
#'   included for newly imported studies. The `cleanup` column records the
#'   outcome for managed downloads.
#'   The out_dir attribute gives the batch report directory.
#' @export
runGWASManifest <- function(manifest, descriptor, data_dir = NULL, out_dir = NULL,
                            knowledge = NULL, defaults = list(), import = TRUE,
                            preview_rows = 1000L, engine = c("R", "cpp"), workers = 1L,
                            scratch_dir = NULL) {
  engine <- match.arg(engine)
  if (length(workers) != 1L || !is.finite(workers) || workers < 1 ||
      workers > .Machine$integer.max || workers != floor(workers)) stop("workers must be a positive integer")
  required <- c("accession", "filename", "url", "md5", "source", "trait",
                "type", "build", "ancestry")
  required_values <- setdiff(required, "url")
  if (is.character(manifest) && length(manifest) == 1L) {
    manifest <- data.table::fread(file = manifest, sep = "\t", data.table = FALSE)
  }
  if (!is.data.frame(manifest) || !nrow(manifest) || anyDuplicated(names(manifest)) ||
      !all(required %in% union(names(manifest), names(defaults)))) {
    stop("Supply a nonempty manifest with columns: ", paste(required, collapse = ", "))
  }
  if (!is.list(defaults) || (length(defaults) &&
      (is.null(names(defaults)) || anyDuplicated(names(defaults)) ||
       !all(names(defaults) %in% .gact_gwas_metadata_fields) ||
       !all(lengths(defaults) == 1L)))) {
    stop("defaults must be a named list of scalar study metadata")
  }
  if (!identical(import, TRUE) && !identical(import, FALSE)) {
    stop("import must be TRUE or FALSE")
  }
  factors <- vapply(manifest, is.factor, logical(1))
  manifest[factors] <- lapply(manifest[factors], as.character)
  for (field in names(defaults)) {
    if (!field %in% names(manifest)) manifest[[field]] <- NA
    value <- manifest[[field]]
    value[is.na(value) | !nzchar(trimws(as.character(value)))] <- defaults[[field]]
    manifest[[field]] <- value
  }
  manifest[required] <- lapply(manifest[required], function(x) trimws(as.character(x)))
  if (!is.null(knowledge)) knowledge <- normalizePath(knowledge, mustWork = TRUE)
  descriptor <- normalizePath(descriptor, winslash = "/", mustWork = TRUE)
  GAlist <- readRDS(descriptor)
  gstat <- GAlist$dirs[["gstat"]]
  if (length(gstat) != 1L || is.na(gstat) || !dir.exists(gstat)) {
    stop("GAlist must point to an existing gstat directory")
  }
  inventory <- file.path(gstat, "GWAS_information.csv")
  inventory_ids <- if (file.exists(inventory)) {
    ids <- read.csv2(inventory, colClasses = "character")$id
    if (is.null(ids)) stop("GWAS_information.csv is missing the id column")
    ids
  } else character()
  study_ids <- as.character(GAlist$study$id)
  if (anyDuplicated(study_ids) || anyDuplicated(inventory_ids) ||
      !setequal(study_ids, as.character(inventory_ids))) {
    stop("GAlist study IDs differ from GWAS_information.csv. Load the current descriptor.")
  }

  database <- dirname(gstat)
  recipes <- file.path(database, "recipes")
  if (is.null(data_dir)) {
    if ("download" %in% names(GAlist$dirs)) data_dir <- GAlist$dirs[["download"]]
    if (!length(data_dir) || is.na(data_dir)) data_dir <- file.path(database, "download")
  }
  if (is.null(out_dir)) {
    out_dir <- tempfile("batch-", tmpdir = file.path(recipes, "reviews"))
  }
  data_dir <- .gact_gwas_output_path(path.expand(data_dir))
  out_dir <- .gact_gwas_output_path(path.expand(out_dir))
  input_key <- .gact_gwas_path_key(data_dir)
  output_key <- .gact_gwas_path_key(out_dir)
  if (output_key == input_key || startsWith(output_key, paste0(input_key, "/"))) {
    stop("out_dir must be outside data_dir")
  }
  if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) {
    stop("out_dir must be new or empty")
  }
  for (path in c(data_dir, out_dir)) {
    if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
      stop("Cannot create directory: ", path)
    }
  }

  data.table::fwrite(manifest, file.path(out_dir, "manifest.tsv"), sep = "\t", na = "NA")
  accessions <- toupper(manifest$accession)
  filenames <- .gact_gwas_path_key(manifest$filename)
  duplicate <- duplicated(accessions) | duplicated(accessions, fromLast = TRUE) |
    duplicated(filenames) | duplicated(filenames, fromLast = TRUE)
  jobs <- vector("list", nrow(manifest))
  report <- data.frame(accession = manifest$accession, filename = manifest$filename,
    status = "needs_review", study_id = NA_character_, recipe = NA_character_,
    issues = NA_character_, warnings = NA_character_, cleanup = NA_character_,
    n_input = NA_integer_, n_final = NA_integer_, n_dropped = NA_integer_,
    n_flipped = NA_integer_, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, , drop = FALSE]
    tryCatch({
      accession <- accessions[[i]]
      if (!grepl("^[A-Z0-9][A-Z0-9._-]*$", accession) || duplicate[[i]]) {
        stop("Invalid or duplicate accession or duplicate filename")
      }
      hits <- .gact_gwas_manifest_matches(GAlist$study, accession)
      if (length(hits) > 1L) stop("Accession matches multiple imported studies")
      if (length(hits)) {
        study_id <- GAlist$study$id[hits]
        study_file <- unname(GAlist$studyfiles[study_id])
        if (length(study_file) != 1L || !file.exists(study_file)) {
          stop("Stored summary statistics are missing for ", study_id)
        }
        old_md5 <- GAlist$study$input_md5[hits]
        if (length(old_md5) && !is.na(old_md5) &&
            !identical(tolower(old_md5), tolower(row$md5))) {
          stop("Accession was imported with a different input checksum")
        }
        report$status[i] <- "already_imported"
        report$study_id[i] <- study_id
        next
      }
      if (anyNA(row[required_values]) || any(!nzchar(unlist(row[required_values])))) {
        stop("Missing required study metadata")
      }
      if (!row$type %in% c("binary", "quantitative") ||
          !row$build %in% c("GRCh37", "GRCh38")) {
        stop("Specify type (binary or quantitative) and build (GRCh37 or GRCh38)")
      }
      build <- .gact_catalog_build(GAlist)
      if (length(build) != 1L || is.na(build) || !identical(build, row$build)) {
        stop("Verify that the study and marker catalog use the same genome build")
      }
      if (!grepl("^[[:alnum:]][[:alnum:]_. -]*\\.(tsv|txt|csv|sumstats|imp)(\\.(gz|bgz|bz2|xz))?$",
                 row$filename, ignore.case = TRUE)) {
        stop("filename must be a basename with a supported summary-statistics extension")
      }
      if (grepl("^GCST[0-9]+$", accession) &&
          (is.null(row$reference) || is.na(row$reference) || !nzchar(row$reference))) {
        row$reference <- .gact_gwas_catalog_url(accession)
      }
      file <- .gact_gwas_manifest_file(row, data_dir)
      pattern <- paste0("^", gsub(".", "[.]", row$filename, fixed = TRUE), "$")
      inspection <- inspectGWAS(data_dir, row, file.path(out_dir, sprintf("study_%05d", i)),
                                file_column = "filename", knowledge = knowledge,
                                preview_rows = preview_rows,
                                pattern = pattern, recursive = FALSE, reader = engine)
      entry <- inspection$review
      if (nrow(entry) != 1L || !identical(entry$file, file) || entry$metadata_matches != 1L) {
        stop("Could not isolate the manifest input in the download directory")
      }
      report$recipe[i] <- entry$recipe
      recipe <- jsonlite::read_json(entry$recipe)
      report$warnings[i] <- paste(unlist(recipe$validation$warnings), collapse = " | ")
      issues <- recipe$inspection_issues
      if (!isTRUE(entry$preview_ok) || !isTRUE(recipe$inspection$profile_reused) || length(issues)) {
        stop(paste(c("Review the preview or verified format", issues), collapse = " | "))
      }
      recipe$accession <- accession
      recipe$input_md5 <- tolower(row$md5)
      recipe$input_url <- row$url
      if (engine == "cpp") .gact_cpp_gwas_mapping(recipe)
      if (!import) {
        report$status[i] <- "ready"
        next
      }
      if (engine == "cpp") {
        jobs[[i]] <- list(recipe = recipe, file = file)
        next
      }
      supplied_fields <- setdiff(.gact_gwas_metadata_fields, c("study_n_total", "study_n_case", "study_n_control"))
      supplied <- Filter(Negate(is.null), recipe[intersect(supplied_fields, names(recipe))])
      result <- do.call(ingestStatDB, c(list(GAlist = GAlist,
        stat = .gact_gwas_read_source(file, skip = 0L), schema = recipe$mapping,
        recipe = recipe, descriptor = descriptor), supplied))
      GAlist <- result$GAlist
      report$status[i] <- "imported"
      report$study_id[i] <- tail(GAlist$study$id, 1L)
      if (!is.null(result$cleanup)) {
        report$cleanup[i] <- result$cleanup$status
        warnings <- c(recipe$validation$warnings, result$cleanup$warning)
        report$warnings[i] <- paste(unlist(warnings), collapse = " | ")
      }
      for (field in c("n_input", "n_final", "n_dropped", "n_flipped")) {
        report[[field]][i] <- result$qc[[field]]
      }
      rm(result)
    }, error = function(e) {
      report$issues[i] <<- conditionMessage(e)
    })
  }
  if (any(lengths(jobs) > 0L)) {
    work_dir <- .gact_gwas_work_dir(scratch_dir)
    on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
    report <- .gact_cpp_gwas_batch(GAlist, descriptor, jobs, report, recipes, work_dir, as.integer(workers))
  }
  data.table::fwrite(report, file.path(out_dir, "manifest_review.tsv"), sep = "\t", na = "NA")
  attr(report, "out_dir") <- out_dir
  invisible(report)
}

.gact_gwas_catalog_url <- function(accession) {
  paste0("https://www.ebi.ac.uk/gwas/studies/", accession)
}

.gact_gwas_manifest_matches <- function(study, accession) {
  explicit <- which(toupper(study$accession) == accession)
  catalog_url <- .gact_gwas_catalog_url(accession)
  legacy <- which(sub("/$", "", as.character(study$reference)) == catalog_url)
  union(explicit, legacy)
}

.gact_gwas_manifest_file <- function(row, data_dir) {
  expected <- tolower(row$md5)
  if (!grepl("^[a-f0-9]{32}$", expected)) stop("Supply the provider's MD5 checksum")
  file <- file.path(data_dir, row$filename)
  cached <- file.exists(file)
  if (!cached) {
    lock <- paste0(file, ".download-lock")
    if (!dir.create(lock, showWarnings = FALSE)) {
      stop("Cannot lock download cache for ", row$filename)
    }
    on.exit(unlink(lock, recursive = TRUE), add = TRUE)
    cached <- file.exists(file)
  }
  candidate <- file
  if (cached) {
    candidate <- normalizePath(file, winslash = "/", mustWork = TRUE)
    if (!identical(dirname(candidate), data_dir) || dir.exists(candidate)) {
      stop("Cached input must be a file inside data_dir")
    }
  } else {
    if (is.na(row$url) || !grepl("^https://", row$url, ignore.case = TRUE)) {
      stop("Missing cached input or HTTPS download URL")
    }
    candidate <- tempfile(".gact-download-", tmpdir = data_dir)
    on.exit(unlink(candidate), add = TRUE)
    curl::curl_download(row$url, candidate, mode = "wb")
  }
  if (!identical(unname(tools::md5sum(candidate)), expected)) {
    stop("Input MD5 checksum mismatch: ", row$filename)
  }
  if (!cached) {
    receipt <- tempfile(".gact-download-receipt-", tmpdir = data_dir)
    on.exit(unlink(receipt, expand = FALSE), add = TRUE)
    jsonlite::write_json(list(url = row$url, md5 = expected), receipt, auto_unbox = TRUE)
    .publish_gwas_import_file(candidate, file)
    .publish_gwas_import_file(receipt, paste0(file, ".gact-download.json"), replace = TRUE)
  }
  normalizePath(file, winslash = "/", mustWork = TRUE)
}

.gact_gwas_cleanup_download <- function(recipe) {
  if (!is.list(recipe)) return(NULL)
  file <- recipe[["file_name"]]
  url <- recipe[["input_url"]]
  md5 <- recipe[["input_md5"]]
  fields <- list(file, url, md5)
  if (!all(lengths(fields) == 1L) || !all(vapply(fields, is.character, logical(1))) ||
      anyNA(unlist(fields)) || !nzchar(file) ||
      !grepl("^https://", url, ignore.case = TRUE) || !grepl("^[a-f0-9]{32}$", md5)) {
    return(NULL)
  }
  receipt <- paste0(file, ".gact-download.json")
  if (!file.exists(receipt)) return(NULL)
  tryCatch({
    lock <- paste0(file, ".download-lock")
    if (!dir.create(lock, showWarnings = FALSE)) stop("the download cache is locked")
    on.exit(unlink(lock, recursive = TRUE, expand = FALSE), add = TRUE)
    if (!identical(jsonlite::read_json(receipt), list(url = url, md5 = md5))) {
      stop("the download receipt does not match")
    }
    if (!file.exists(file)) return(list(status = "already_absent"))
    if (dir.exists(file) || !identical(normalizePath(file, winslash = "/", mustWork = TRUE), file) ||
        !identical(unname(tools::md5sum(file)), md5)) {
      stop("the input path or checksum changed")
    }
    if (!suppressWarnings(file.remove(file))) stop("the file could not be removed")
    list(status = "removed")
  }, error = function(e) {
    list(status = "retained", warning = paste0("Download cleanup skipped: ", conditionMessage(e)))
  })
}
