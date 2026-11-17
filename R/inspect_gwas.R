.gact_gwas_metadata_fields <- c("source", "trait", "type", "category", "efo",
  "gender", "ancestry", "build", "n", "ncase", "ncontrol", "reference",
  "comments", "effect_scale", "p_scale", "strand",
  "study_n_total", "study_n_case", "study_n_control")

#' Inspect local GWAS files and prepare ingestion recipes
#'
#' Reads bounded leading previews, joins files to study metadata, and writes
#' recipes, a review table, reusable format observations, and an ingestion script.
#' The script is executed separately with the path to a GAlist RDS database.
#' Preview success describes only the sampled records. Alignment remains the
#' responsibility of the existing ingestion machinery.
#' All recipes must pass preview validation before the script starts writing.
#' After correcting a failed preview, inspect the file again to regenerate its
#' recipe. Omitted optional metadata uses the ingestion API defaults.
#' The preview and the script support plain text, gzip, bzip2, and XZ.
#' The script expands compressed inputs to temporary files for full-file reads.
#'
#' Format knowledge is explicit and persistent. Edit a generated formats.json
#' entry to correct its mapping and scales, add evidence, and set verified=true.
#' Supply that file as knowledge on a later run. Only verified entries with an
#' exact source and header match are reused. Source should identify a release
#' when conventions differ between releases. No statistical model is trained.
#'
#' @param data_dir Directory containing downloaded summary statistics.
#' @param metadata Path to a tab-separated metadata file or a data.frame.
#' @param out_dir New or empty output directory, outside data_dir.
#' @param file_column Metadata column containing file paths or basenames.
#'   If NULL, a unique filename, file_name, file, or path column is detected.
#'   Paths are absolute or relative to `data_dir`. A bare basename must
#'   uniquely identify a file among the inputs.
#' @param metadata_mapping Named character vector mapping gact metadata fields
#'   to metadata columns. Canonically named columns are used by default.
#' @param defaults Named list of study metadata supplied by the caller. Values
#'   in mapped metadata columns take precedence. Build is never inferred from
#'   a directory name.
#' @param knowledge Optional path to formats.json from an earlier inspection.
#' @param preview_rows Maximum number of leading data lines per file.
#' @param skip Number of leading lines before the header.
#' @param pattern Regular expression selecting supported text files.
#' @param recursive Search subdirectories.
#' @param reader Use the R preview reader or the C++ reader for gzip TSVs.
#' @return Invisibly, a list with review, recipes, formats, metadata_columns,
#'   unmatched_metadata, and script paths or values.
#' @export
inspectGWAS <- function(data_dir, metadata, out_dir, file_column = NULL,
                     metadata_mapping = NULL, defaults = list(), knowledge = NULL,
                     preview_rows = 1000L, skip = 0L,
                     pattern = "\\.(tsv|txt|csv|sumstats|imp)(\\.(gz|bgz|bz2|xz))?$",
                     recursive = TRUE, reader = c("R", "cpp")) {
  reader <- match.arg(reader)
  data_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)
  if (!dir.exists(data_dir)) stop("data_dir must be a directory")
  for (value in list(preview_rows, skip)) {
    if (length(value) != 1L || !is.finite(value) || value != floor(value)) {
      stop("preview_rows and skip must be finite integers")
    }
  }
  if (preview_rows < 1 || skip < 0) stop("Invalid preview_rows or skip")
  metadata_file <- character()
  if (is.character(metadata) && length(metadata) == 1L) {
    metadata_file <- normalizePath(metadata, winslash = "/", mustWork = TRUE)
    metadata <- data.table::fread(file = metadata, sep = "\t", data.table = FALSE)
  }
  if (!is.data.frame(metadata)) stop("metadata must be a TSV path or data.frame")
  if (is.null(file_column)) {
    candidates <- names(metadata)[tolower(names(metadata)) %in%
      c("filename", "file_name", "file", "path")]
    if (length(candidates) != 1L) {
      stop("Specify file_column. Metadata columns: ", paste(names(metadata), collapse = ", "))
    }
    file_column <- candidates
  }
  if (length(file_column) != 1L || !file_column %in% names(metadata)) {
    stop("file_column must name a metadata column")
  }
  fields <- .gact_gwas_metadata_fields
  if (is.null(metadata_mapping)) {
    shared <- intersect(fields, names(metadata))
    metadata_mapping <- stats::setNames(shared, shared)
  }
  if (length(metadata_mapping) && (is.null(names(metadata_mapping)) ||
      any(!names(metadata_mapping) %in% fields) ||
      any(!metadata_mapping %in% names(metadata)))) {
    stop("metadata_mapping must map supported metadata fields to existing columns")
  }
  if (length(defaults) && (is.null(names(defaults)) || any(!names(defaults) %in% fields))) {
    stop("defaults must be a named list of supported metadata fields")
  }
  profiles <- list()
  if (!is.null(knowledge)) profiles <- jsonlite::read_json(knowledge)$formats
  files <- sort(list.files(data_dir, pattern = pattern, full.names = TRUE,
                           recursive = recursive, ignore.case = TRUE))
  files <- files[!dir.exists(files)]
  if (!length(files)) stop("No matching summary-statistics files found")
  files <- normalizePath(files, winslash = "/", mustWork = TRUE)
  files <- files[!.gact_gwas_path_key(files) %in% .gact_gwas_path_key(metadata_file)]
  if (!length(files)) stop("No summary-statistics files beyond the metadata table")
  input_key <- .gact_gwas_path_key(data_dir)
  check_output_path <- function(path) {
    key <- .gact_gwas_path_key(path)
    if (key == input_key || startsWith(key, paste0(sub("/+$", "", input_key), "/"))) {
      stop("out_dir must be outside data_dir")
    }
  }
  out_dir <- .gact_gwas_output_path(path.expand(out_dir))
  check_output_path(out_dir)
  if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) {
    stop("out_dir must be new or empty; earlier recipes are preserved")
  }
  aliases <- read_aliases_yaml(must_exist = TRUE)
  if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create out_dir: ", out_dir)
  }
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  check_output_path(out_dir)
  dir.create(file.path(out_dir, "recipes"))

  metadata_paths <- gsub("\\\\", "/", as.character(metadata[[file_column]]))
  has_path <- !is.na(metadata_paths) & nzchar(metadata_paths)
  absolute <- grepl("^(/|[A-Za-z]:)", metadata_paths)
  bare <- has_path & !absolute & !grepl("/", metadata_paths, fixed = TRUE)
  relative <- has_path & !bare & !absolute
  metadata_paths[relative] <- file.path(data_dir, metadata_paths[relative])
  qualified <- has_path & !bare & file.exists(metadata_paths)
  metadata_paths[qualified] <- normalizePath(metadata_paths[qualified], winslash = "/", mustWork = TRUE)
  metadata_keys <- .gact_gwas_path_key(metadata_paths)
  file_keys <- .gact_gwas_path_key(files)
  basenames <- basename(file_keys)
  recipes <- formats <- rows <- list()
  used_metadata <- integer()
  present <- function(x) length(x) == 1L && !is.na(x) && nzchar(as.character(x))

  for (i in seq_along(files)) {
    file <- files[[i]]
    id <- sprintf("study_%05d", i)
    matches <- which(qualified & metadata_keys == file_keys[[i]])
    issues <- character()
    recipe <- defaults
    if (!length(matches)) {
      by_name <- which(bare & metadata_keys == basenames[[i]])
      if (length(by_name) && sum(basenames == basenames[[i]]) > 1L) {
        issues <- c(issues, "Ambiguous basename across data directories")
      } else {
        matches <- by_name
      }
    }
    if (length(matches) == 1L) {
      used_metadata <- c(used_metadata, matches)
      for (field in names(metadata_mapping)) {
        value <- metadata[[metadata_mapping[[field]]]][matches]
        if (present(value)) recipe[[field]] <- value
      }
    } else {
      issues <- c(issues, paste("Metadata matches:", length(matches)))
    }
    recipe$file_name <- file
    result <- tryCatch({
      preview <- if (reader == "cpp") .gact_cpp_gwas_preview(file, preview_rows, skip) else
        .gact_gwas_preview(file, preview_rows, skip)
      headers <- names(preview)
      det <- detectStatSchema(preview, sample_n = NULL, alias_dict = aliases)
      recipe$mapping <- as.list(det$mapping)
      recipe$effect_scale <- if (present(recipe$effect_scale)) recipe$effect_scale else "auto"
      recipe$p_scale <- if (present(recipe$p_scale)) recipe$p_scale else "auto"
      profile_matches <- Filter(function(p) {
        isTRUE(p$verified) && present(p$evidence) && present(recipe$source) &&
          identical(p$source, as.character(recipe$source)) &&
          identical(unlist(p$headers, use.names = FALSE), headers)
      }, profiles)
      if (length(profile_matches) > 1L) stop("Multiple verified format profiles match")
      reused <- length(profile_matches) == 1L
      if (reused) {
        profile <- profile_matches[[1L]]
        recipe$mapping <- lapply(profile$mapping, function(value) {
          if (!length(value)) return(NA_character_)
          if (!is.character(value) || length(value) != 1L) stop("Invalid profile mapping")
          value
        })
        for (field in c("effect_scale", "p_scale")) {
          if (present(profile[[field]])) recipe[[field]] <- profile[[field]]
        }
      }
      mapped <- unlist(recipe$mapping, use.names = FALSE)
      missing_columns <- setdiff(mapped[!is.na(mapped) & nzchar(mapped)], headers)
      if (length(missing_columns)) stop("Mapped columns absent: ", paste(missing_columns, collapse = ", "))
      warnings <- character()
      normalized <- withCallingHandlers(normalizeStatSchema(preview,
        schema = recipe$mapping, effect_scale = recipe$effect_scale,
        p_scale = recipe$p_scale), warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      validation <- validateStatSchema(normalized)
      recipe$validation <- validation
      recipe$inspection <- list(headers = headers, preview_rows = nrow(preview),
        skip = skip, size = unname(file.info(file)$size),
        mtime = as.numeric(file.info(file)$mtime), profile_reused = reused,
        detection = det, warnings = unique(c(warnings, validation$warnings)))
      if (!reused) {
        issues <- c(issues, "Review inferred mapping and source conventions", det$warnings)
      }
      issues <- c(issues, validation$errors, warnings)
      for (field in c("source", "trait", "type", "build", "ancestry")) {
        if (!present(recipe[[field]])) issues <- c(issues, paste("Missing", field))
      }
      profile <- list(source = recipe$source, headers = headers,
        mapping = recipe$mapping, effect_scale = recipe$effect_scale,
        p_scale = recipe$p_scale, verified = reused,
        evidence = if (reused) profile_matches[[1L]]$evidence else "",
        example_file = file)
      formats[[id]] <- profile
      list(recipe = recipe, ok = validation$ok, rows = nrow(preview))
    }, error = function(e) {
      issues <<- c(issues, conditionMessage(e))
      list(recipe = recipe, ok = FALSE, rows = NA_integer_)
    })
    recipe <- result$recipe
    recipe$inspection_issues <- unique(issues)
    recipe_path <- file.path(out_dir, "recipes", paste0(id, ".json"))
    jsonlite::write_json(recipe, recipe_path, pretty = TRUE, auto_unbox = TRUE, na = "null", digits = NA)
    recipes[[id]] <- recipe_path
    rows[[i]] <- data.frame(id = id, file = file, metadata_matches = length(matches),
      preview_rows = result$rows, preview_ok = result$ok, recipe = recipe_path,
      issues = paste(unique(issues), collapse = " | "), stringsAsFactors = FALSE)
  }
  review <- do.call(rbind, rows)
  signatures <- vapply(formats, function(p) jsonlite::toJSON(p[c("source", "headers")],
    auto_unbox = TRUE), character(1))
  formats <- formats[!duplicated(signatures)]
  jsonlite::write_json(list(formats = formats), file.path(out_dir, "formats.json"),
    pretty = TRUE, auto_unbox = TRUE, na = "null", digits = NA)
  data.table::fwrite(review, file.path(out_dir, "review.tsv"), sep = "\t")
  unmatched <- metadata[setdiff(seq_len(nrow(metadata)), used_metadata), , drop = FALSE]
  data.table::fwrite(unmatched, file.path(out_dir, "unmatched_metadata.tsv"), sep = "\t")
  script <- file.path(out_dir, "ingest.R")
  .gact_write_ingestion_script(script, unname(unlist(recipes)), skip)
  invisible(list(review = review, recipes = recipes, formats = formats,
    metadata_columns = names(metadata), unmatched_metadata = unmatched, script = script))
}

.gact_gwas_path_key <- function(path) {
  if (.Platform$OS.type == "windows") tolower(path) else path
}

.gact_gwas_output_path <- function(path) {
  if (file.exists(path)) return(normalizePath(path, winslash = "/", mustWork = TRUE))
  parent <- dirname(path)
  if (identical(parent, path)) stop("Cannot resolve out_dir: ", path)
  parent <- .gact_gwas_output_path(parent)
  leaf <- basename(path)
  switch(leaf, "." = parent, ".." = dirname(parent), file.path(parent, leaf))
}

.gact_gwas_connection <- function(path) {
  switch(tolower(tools::file_ext(path)), gz = gzfile, bgz = gzfile,
         bz2 = bzfile, xz = xzfile, base::file)
}

.gact_gwas_preview <- function(file, rows, skip) {
  connection <- .gact_gwas_connection(file)(file, "rt")
  on.exit(close(connection))
  if (skip) readLines(connection, n = skip, warn = FALSE)
  lines <- readLines(connection, n = rows + 1L, warn = FALSE)
  if (length(lines) < 2L) stop("Preview has no data rows")
  data.table::fread(text = paste(lines, collapse = "\n"), nrows = rows,
                    data.table = FALSE, check.names = FALSE)
}

.gact_gwas_read_source <- function(path, skip) {
  reader <- .gact_gwas_connection(path)
  if (!identical(reader, base::file)) {
    unpacked <- tempfile("gact-gwas-")
    on.exit(unlink(unpacked))
    R.utils::decompressFile(path, destname = unpacked, ext = tools::file_ext(path),
                            FUN = reader, remove = FALSE)
    path <- unpacked
  }
  data.table::fread(file = path, skip = skip, data.table = FALSE, check.names = FALSE)
}

.gact_write_ingestion_script <- function(path, recipes, skip) {
  paths <- paste(utils::capture.output(dput(recipes)), collapse = "\n")
  writeLines(c(
    "# Generated by inspectGWAS. Review and correct the recipe JSON files first.",
    "# Run with Rscript ingest.R /absolute/path/to/GAlist.rds",
    "library(gact)",
    "args <- commandArgs(trailingOnly = TRUE)",
    "if (length(args) != 1L) stop('Supply the GAlist RDS path')",
    paste0("recipe_paths <- ", paths),
    "recipes <- lapply(recipe_paths, jsonlite::read_json, simplifyVector = FALSE)",
    "# Check all recipes before the first database write.",
    "for (r in recipes) {",
    "  required <- c('trait', 'type', 'build', 'ancestry', 'source')",
    "  for (field in required) {",
    "    if (is.null(r[[field]]) || length(r[[field]]) != 1L ||",
    "        is.na(r[[field]]) || !nzchar(r[[field]])) stop('Missing ', field, ' in ', r$file_name)",
    "  }",
    "  if (is.null(r$inspection) || !isTRUE(r$validation$ok)) {",
    "    stop('File preview failed; inspect again: ', r$file_name,",
    "         '. ', paste(r$validation$errors, collapse = ' | '))",
    "  }",
    "  info <- file.info(r$file_name)",
    "  if (is.na(info$size) || info$size != r$inspection$size ||",
    "      abs(as.numeric(info$mtime) - r$inspection$mtime) > 1e-5) stop('File changed; inspect again: ', r$file_name)",
    "}",
    "GAlist <- readRDS(args[[1L]])",
    "for (r in recipes) {",
    paste0("  stat <- gact:::.gact_gwas_read_source(r$file_name, skip = ", skip, ")"),
    "  fields <- c('source', 'trait', 'type', 'category', 'efo', 'gender', 'ancestry',",
    "    'build', 'n', 'ncase', 'ncontrol', 'reference', 'comments', 'strand',",
    "    'effect_scale', 'p_scale')",
    "  supplied <- Filter(Negate(is.null), r[intersect(fields, names(r))])",
    "  result <- do.call(ingestStatDB, c(list(GAlist = GAlist, stat = stat, schema = r$mapping,",
    "    recipe = r), supplied))",
    "  GAlist <- result$GAlist",
    "  saveRDS(GAlist, args[[1L]], compress = FALSE)",
    "}"), path)
}
