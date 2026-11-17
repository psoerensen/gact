#' Scan Pan-UK Biobank studies for GWAS ingestion
#'
#' Match selected phenotypes to verified GRCh38 GWAS Catalog files. Use the
#' C++ reader on a bounded gzip prefix and prepare recipes and an import manifest.
#' The scanner does not open a database. The importer verifies the full file MD5.
#'
#' @param selection Data frame with Pan-UKB filename and analysis (EUR, hq, or
#'   all). An optional accession column supplies an explicit Catalog match.
#' @param out_dir New or empty directory for the recipes and scan reports.
#' @param phenotype_manifest Official compressed Pan-UKB manifest URL or local path.
#' @param preview_rows Number of complete leading rows to inspect.
#' @param max_preview_bytes Maximum compressed bytes requested per file.
#' @return Paths to the manifest and format profiles, and a per-study review table.
#' @export
scanGWAS <- function(selection, out_dir,
    phenotype_manifest = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz",
    preview_rows = 1000L, max_preview_bytes = 1024L * 1024L) {
  if (!requireNamespace("yaml", quietly = TRUE)) stop("Install yaml to read Catalog metadata")
  required <- c("filename", "analysis")
  if (!is.data.frame(selection) || !nrow(selection) || !all(required %in% names(selection)) ||
      anyNA(selection[required]) || anyDuplicated(selection[required]) ||
      !all(selection$analysis %in% c("EUR", "hq", "all"))) {
    stop("Supply a unique Pan-UKB filename and analysis per selection")
  }
  if (!"accession" %in% names(selection)) selection$accession <- NA_character_
  selection$accession <- as.character(selection$accession)
  explicit <- !is.na(selection$accession) & nzchar(selection$accession)
  if (!all(grepl("^GCST[0-9]{8}$", selection$accession[explicit]))) stop("Invalid Catalog accession")
  for (value in list(preview_rows, max_preview_bytes)) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value >= .Machine$integer.max || value != floor(value)) {
      stop("Preview limits must be positive integers")
    }
  }
  if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) {
    stop("out_dir must be new or empty")
  }
  phenotype <- .gact_read_panukb_manifest(phenotype_manifest)
  at <- match(selection$filename, phenotype$filename)
  catalog <- if (any(!explicit & !is.na(at))) {
    tryCatch(.gact_panukb_catalog_studies(), error = function(e) e)
  } else list()
  if (!dir.create(file.path(out_dir, "recipes"), recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create scan recipe directory")
  }
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  review <- data.frame(filename = selection$filename, analysis = selection$analysis,
    accession = selection$accession, status = "needs_review", preview_bytes = NA_integer_,
    preview_rows = NA_integer_, recipe = NA_character_, issues = NA_character_)
  rows <- formats <- vector("list", nrow(selection))
  for (i in seq_len(nrow(selection))) {
    result <- tryCatch({
      if (is.na(at[i])) stop("Filename is absent from the Pan-UKB phenotype manifest")
      row <- phenotype[at[i], , drop = FALSE]
      accession <- if (explicit[i]) selection$accession[i] else {
        if (inherits(catalog, "error")) stop(catalog)
        .gact_panukb_candidates(row, phenotype, catalog)
      }
      picked <- .gact_panukb_scan_match(row, accession, selection$analysis[i])
      review$accession[i] <- picked$accession
      prefix <- tempfile(fileext = ".tsv.gz")
      sample <- tryCatch({
        bytes <- .gact_fetch_gwas_prefix(picked$url, prefix, max_preview_bytes)
        review$preview_bytes[i] <- bytes
        .gact_cpp_gwas_preview(prefix, preview_rows, 0L, prefix = TRUE)
      }, finally = unlink(prefix))
      review$preview_rows[i] <- nrow(sample)
      if (nrow(sample) < preview_rows) stop("Byte limit reached before the requested number of complete rows")
      profile <- .gact_panukb_catalog_profile(picked$source)
      recipe <- .gact_gwas_scan_recipe(picked, sample, bytes, profile)
      path <- file.path(out_dir, "recipes", sprintf("study_%05d.json", i))
      jsonlite::write_json(recipe, path, auto_unbox = TRUE, pretty = TRUE, na = "null", digits = NA)
      review$recipe[i] <- path
      if (length(recipe$inspection_issues)) stop(paste(recipe$inspection_issues, collapse = " | "))
      list(manifest = picked, format = profile)
    }, error = function(e) e)
    if (inherits(result, "error")) {
      review$issues[i] <- conditionMessage(result)
    } else {
      review$status[i] <- "ready"
      rows[[i]] <- result$manifest
      formats[[i]] <- result$format
    }
  }
  accepted <- which(review$status == "ready")
  for (field in c("accession", "md5")) {
    values <- vapply(rows[accepted], function(row) row[[field]], character(1))
    duplicate <- accepted[duplicated(values) | duplicated(values, fromLast = TRUE)]
    review$status[duplicate] <- "needs_review"
    review$issues[duplicate] <- paste("Multiple selections resolved to the same", field)
  }
  accepted <- which(review$status == "ready")
  manifest <- file.path(out_dir, "manifest.tsv")
  knowledge <- file.path(out_dir, "formats.json")
  report <- file.path(out_dir, "scan_review.tsv")
  prepared <- if (length(accepted)) do.call(rbind, rows[accepted]) else
    as.data.frame(stats::setNames(rep(list(character()), 9L),
      c("accession", "filename", "url", "md5", "source", "trait", "type", "build", "ancestry")))
  data.table::fwrite(prepared, manifest, sep = "\t", na = "NA")
  data.table::fwrite(review, report, sep = "\t", na = "NA")
  jsonlite::write_json(list(formats = unique(formats[accepted])), knowledge,
    auto_unbox = TRUE, pretty = TRUE)
  invisible(list(manifest = manifest, knowledge = knowledge, review = review, report = report))
}

.gact_panukb_scan_match <- function(row, accession, analysis) {
  matches <- lapply(accession, function(id) {
    tryCatch(.gact_panukb_catalog_row(row, id, analysis, metadata_dir = NULL),
      gact_panukb_mismatch = function(e) e)
  })
  matched <- !vapply(matches, inherits, logical(1), "gact_panukb_mismatch")
  if (!any(matched)) {
    detail <- if (length(matches)) paste(unique(vapply(matches, conditionMessage, "")), collapse = " | ")
    stop("No verified Catalog match", if (length(detail)) paste0(": ", detail))
  }
  if (sum(matched) > 1L) stop("Multiple verified Catalog matches: ", paste(accession[matched], collapse = ", "))
  matches[[which(matched)]]
}

.gact_gwas_scan_recipe <- function(study, sample, bytes, profile) {
  recipe <- as.list(study)
  recipe$input_md5 <- study$md5
  recipe$input_url <- study$url
  reused <- identical(names(sample), profile$headers)
  recipe$inspection <- list(headers = names(sample), preview_rows = nrow(sample),
    preview_bytes = bytes, remote = TRUE, skip = 0L, profile_reused = reused)
  if (reused) {
    recipe$mapping <- profile$mapping
    recipe$effect_scale <- profile$effect_scale
    recipe$p_scale <- profile$p_scale
    recipe$validation <- validateStatSchema(normalizeStatSchema(sample, recipe$mapping,
      effect_scale = recipe$effect_scale, p_scale = recipe$p_scale))
    recipe$inspection_issues <- recipe$validation$errors
    if (!length(recipe$inspection_issues)) .gact_cpp_gwas_mapping(recipe)
  } else {
    guesses <- detectStatSchema(sample, sample_n = NULL)
    recipe$mapping <- as.list(guesses$mapping)
    recipe$inspection_issues <- c("Unexpected Catalog header. Review the inferred mapping", guesses$warnings)
  }
  recipe
}

.gact_panukb_catalog_studies <- function() {
  base <- "https://www.ebi.ac.uk/gwas/api/v2/publications/40968291/studies"
  pages <- list()
  repeat {
    page <- length(pages)
    response <- curl::curl_fetch_memory(paste0(base, "?size=500&page=", page),
      handle = curl::new_handle(timeout = 45))
    if (response$status_code != 200L) stop("Could not read the Pan-UKB Catalog study listing")
    data <- jsonlite::fromJSON(rawToChar(response$content), simplifyVector = FALSE)
    studies <- data[["_embedded"]]$studies
    if (!is.list(studies) || !length(studies) || !isTRUE(data$page$number == page) ||
        !isTRUE(data$page$totalPages > page)) stop("Unexpected Catalog pagination")
    pages[[page + 1L]] <- studies
    if (page + 1L == data$page$totalPages) break
  }
  unlist(pages, recursive = FALSE, use.names = FALSE)
}

.gact_panukb_candidates <- function(row, phenotype, catalog) {
  key <- .gact_panukb_trait_key(row$description)
  same_label <- .gact_panukb_trait_key(phenotype$description) == key &
    phenotype$pheno_sex == row$pheno_sex
  if (sum(same_label, na.rm = TRUE) > 1L) {
    stop("Several Pan-UKB files share this trait label and sex. Supply a reviewed Catalog accession")
  }
  hit <- Filter(function(study) {
    label <- study$reportedTrait
    if (!isTRUE(.gact_panukb_trait_key(label) == key) ||
        !isTRUE(grepl("^GCST[0-9]{8}$", study$accessionId)) ||
        !isTRUE(grepl("^https?://", study$summaryStatistics))) return(FALSE)
    if (!grepl("\\(phecode [^)]+\\)", label, ignore.case = TRUE)) return(TRUE)
    code <- sub(".*\\(phecode ([^)]+)\\).*", "\\1", label, ignore.case = TRUE)
    identical(row$trait_type, "phecode") && identical(code, row$phenocode)
  }, catalog)
  unique(vapply(hit, function(study) study$accessionId, character(1)))
}

.gact_fetch_gwas_prefix <- function(url, path, max_bytes) {
  connection <- file(path, "wb")
  on.exit(close(connection), add = TRUE)
  count <- 0
  magic <- raw()
  receive <- function(chunk) {
    count <<- count + length(chunk)
    if (count > max_bytes) stop("Catalog server exceeded the preview byte limit")
    if (length(magic) < 2L) magic <<- head(c(magic, chunk), 2L)
    writeBin(chunk, connection)
    invisible(NULL)
  }
  handle <- curl::new_handle(range = sprintf("0-%d", as.integer(max_bytes) - 1L), timeout = 45)
  response <- curl::curl_fetch_stream(url, receive, handle = handle)
  headers <- strsplit(rawToChar(response$headers), "\r?\n")[[1L]]
  range <- grep("^content-range:[[:space:]]*bytes 0-[0-9]+/[0-9]+$", trimws(headers),
                ignore.case = TRUE, value = TRUE)
  if (response$status_code != 206L || length(range) != 1L || !count ||
      as.numeric(sub("^.* 0-([0-9]+)/.*$", "\\1", range)) + 1 != count) {
    stop("Catalog server did not return a valid bounded byte range")
  }
  if (!identical(magic, as.raw(c(31L, 139L)))) stop("Catalog preview is not gzip data")
  count
}
