#' Prepare a Pan-UK Biobank manifest for harmonised Catalog files
#'
#' Resolve reviewed Pan-UKB phenotype filename and Catalog accession pairs
#' against their provider metadata. Keep cohort counts as study metadata.
#'
#' @param selection Data frame with filename, accession, and analysis (EUR, hq,
#'   or all). Each row identifies a reviewed Catalog match.
#' @param out_dir New or empty directory for the manifest and format profiles.
#' @param phenotype_manifest Official compressed Pan-UKB manifest URL or local path.
#' @param metadata_dir Optional directory with Catalog per-file YAML and
#'   accession-named .md5 lists. Otherwise fetch them from the Catalog.
#' @return Paths to the gact manifest and verified format profiles.
#' @export
panUKBManifest <- function(selection, out_dir,
    phenotype_manifest = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz",
    metadata_dir = NULL) {
  if (!requireNamespace("yaml", quietly = TRUE)) stop("Install yaml to read Catalog metadata")
  required <- c("filename", "accession", "analysis")
  if (!is.data.frame(selection) || !nrow(selection) ||
      !all(required %in% names(selection)) || anyNA(selection[required]) ||
      anyDuplicated(selection$accession) ||
      !all(grepl("^GCST[0-9]{8}$", selection$accession)) ||
      !all(selection$analysis %in% c("EUR", "hq", "all"))) stop("Supply a unique Catalog accession and analysis per selection")
  if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) {
    stop("out_dir must be new or empty")
  }
  phenotype <- .gact_read_panukb_manifest(phenotype_manifest)
  rows <- phenotype[match(selection$filename, phenotype$filename), , drop = FALSE]
  if (anyNA(rows$filename)) stop("A selected filename is absent from the Pan-UKB manifest")
  output <- lapply(seq_len(nrow(rows)), function(i) {
    .gact_panukb_catalog_row(rows[i, ], selection$accession[i], selection$analysis[i], metadata_dir)
  })
  result <- do.call(rbind, output)
  if (anyDuplicated(result$md5)) stop("Selected Catalog accessions share an input checksum")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  manifest <- file.path(out_dir, "manifest.tsv")
  knowledge <- file.path(out_dir, "formats.json")
  format <- .gact_panukb_catalog_profile(result$source[1L])
  data.table::fwrite(result, manifest, sep = "\t", na = "NA")
  jsonlite::write_json(list(formats = list(format)), knowledge, auto_unbox = TRUE, pretty = TRUE)
  invisible(list(manifest = manifest, knowledge = knowledge))
}

.gact_read_panukb_manifest <- function(phenotype_manifest) {
  manifest_file <- phenotype_manifest
  if (grepl("^https://", manifest_file)) {
    manifest_file <- tempfile(fileext = ".bgz")
    on.exit(unlink(manifest_file), add = TRUE)
    curl::curl_download(phenotype_manifest, manifest_file, quiet = TRUE)
  }
  connection <- gzfile(manifest_file, "rt")
  phenotype <- tryCatch(read.delim(connection, colClasses = "character", quote = "",
    check.names = FALSE, na.strings = "NA"), finally = close(connection))
  if (anyDuplicated(phenotype$filename)) stop("Pan-UKB manifest contains duplicate filenames")
  phenotype
}

.gact_panukb_catalog_profile <- function(source) {
  fields <- c(chr = "chromosome", pos = "base_pair_location", ea = "effect_allele",
    nea = "other_allele", b = "beta", seb = "standard_error",
    eaf = "effect_allele_frequency", p = "p_value", marker = "rsid")
  header <- c(unname(fields), "hm_coordinate_conversion", "hm_code", "variant_id")
  list(source = source, headers = header, mapping = as.list(fields),
    effect_scale = "beta", p_scale = "p", verified = TRUE,
    evidence = "GWAS-SSF v1.0 column definitions and verified harmonised GRCh38 Catalog metadata")
}

.gact_panukb_trait_key <- function(x) {
  x <- gsub(" \\(phecode [^)]+\\)| levels$", "", tolower(x))
  gsub("[^a-z0-9]", "", gsub("haemoglobin", "hemoglobin", x, fixed = TRUE))
}
.gact_panukb_pop_code <- function(sample) {
  categories <- tolower(unlist(sample$sample_ancestry_category))
  keys <- c(EUR = "^european$", CSA = "south asian|central asian",
    EAS = "east asian", MID = "middle eastern", AFR = "african unspecified|sub-saharan african",
    AMR = "hispanic|latin american|admixed american")
  match <- names(keys)[vapply(keys, function(x) any(grepl(x, categories)), logical(1))]
  if (length(match) != 1L) stop("Cannot map Catalog sample ancestry")
  match
}

.gact_panukb_catalog_row <- function(row, id, analysis, metadata_dir) {
  file <- paste0(id, ".h.tsv.gz")
  num <- as.integer(sub("^GCST", "", id))
  upper <- ceiling(num / 1000) * 1000
  range <- sprintf("GCST%08d-GCST%08d", upper - 999L, upper)
  base <- paste0("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/", range, "/", id, "/harmonised/")
  url <- paste0(base, file)
  fetch <- function(remote, local) {
    if (is.null(metadata_dir)) {
      response <- curl::curl_fetch_memory(remote, handle = curl::new_handle(timeout = 45))
      if (response$status_code != 200L) stop("Cannot fetch Catalog metadata: ", remote)
      rawToChar(response$content)
    } else paste(readLines(file.path(metadata_dir, local), warn = FALSE), collapse = "\n")
  }
  meta <- yaml::yaml.load(fetch(paste0(url, "-meta.yaml"), paste0(file, "-meta.yaml")), eval.expr = FALSE)
  if (!identical(meta$gwas_id, id) || !identical(meta$data_file_name, file) ||
      !identical(meta$genome_assembly, "GRCh38") || !isTRUE(meta$is_harmonised) ||
      !isTRUE(grepl("^GWAS-SSF", meta$file_type)) ||
      !any(grepl("pan.ukbb.broadinstitute.org", meta$author_notes, fixed = TRUE)) ||
      !identical(.gact_panukb_trait_key(unlist(meta$trait_description)), .gact_panukb_trait_key(row$description))) {
    .gact_panukb_mismatch("Catalog metadata does not match the selected Pan-UKB phenotype: ", id)
  }
  sexes <- c(both_sexes = "combined", females = "female", males = "male")
  if (!identical(meta$sex, unname(sexes[row$pheno_sex]))) .gact_panukb_mismatch("Catalog sex scope differs: ", id)
  pops <- switch(analysis, EUR = "EUR", hq = strsplit(row$pops_pass_qc, ",", fixed = TRUE)[[1L]],
                 all = strsplit(row$pops, ",", fixed = TRUE)[[1L]])
  catalog_pops <- vapply(meta$samples, .gact_panukb_pop_code, character(1))
  if (!length(pops) || anyNA(pops) || anyDuplicated(catalog_pops) ||
      !setequal(pops, catalog_pops)) .gact_panukb_mismatch("Catalog ancestry scope differs: ", id)
  binary <- !row$trait_type %in% c("continuous", "biomarkers")
  cases <- as.numeric(unlist(row[paste0("n_cases_", catalog_pops)]))
  controls <- if (binary) as.numeric(unlist(row[paste0("n_controls_", catalog_pops)])) else rep(0, length(cases))
  sizes <- vapply(meta$samples, function(sample) as.numeric(sample$sample_size), numeric(1))
  if (anyNA(c(cases, controls, sizes)) || !all(sizes == cases + controls) ||
      !all(vapply(meta$samples, function(sample) identical(sample$case_control_study, binary), logical(1)))) {
    .gact_panukb_mismatch("Catalog cohort counts differ from the selected Pan-UKB analysis: ", id)
  }
  listing <- strsplit(fetch(paste0(base, "md5sum.txt"), paste0(id, ".md5")), "\n", fixed = TRUE)[[1L]]
  entries <- strsplit(trimws(listing), "[[:space:]]+")
  match <- Filter(function(x) length(x) == 2L && sub("^\\*", "", x[2L]) == file, entries)
  checksums <- vapply(match, function(x) x[1L], character(1))
  checksum <- tolower(meta$data_file_md5sum)
  if (length(checksum) != 1L || !grepl("^[0-9a-f]{32}$", checksum) ||
      (length(checksums) && !identical(tolower(checksums), checksum))) {
    stop("Catalog checksum metadata is invalid or disagrees with the checksum list: ", id)
  }
  data.frame(accession = id, filename = file, url = url, md5 = checksum,
    source = "Pan-UK Biobank Catalog harmonised GWAS-SSF v1.0", trait = row$description,
    type = if (binary) "binary" else "quantitative", build = "GRCh38",
    ancestry = if (length(pops) == 1L) pops else "Mixed",
    gender = switch(row$pheno_sex, both_sexes = "both", females = "females", males = "males"),
    category = row$category, strand = "forward", study_n_total = sum(sizes),
    study_n_case = if (binary) sum(cases) else NA_real_,
    study_n_control = if (binary) sum(controls) else NA_real_,
    comments = paste0("Pan-UKB ", row$filename, "; analysis=", analysis, "; groups=", paste(pops, collapse = ",")),
    reference = .gact_gwas_catalog_url(id), stringsAsFactors = FALSE)
}

.gact_panukb_mismatch <- function(...) {
  stop(errorCondition(paste0(...), class = "gact_panukb_mismatch"))
}
