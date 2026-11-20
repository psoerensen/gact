# Study discovery and manifests --------

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

.gact_gwas_n_columns <- function(headers) {
  headers[gsub("[^a-z0-9]", "", tolower(headers)) %in%
    c("n", "nsamplesize", "samplesize", "ntotal", "totalsamplesize")]
}

.gact_catalog_preview_profile <- function(sample, source) {
  profile <- .gact_panukb_catalog_profile(source)
  mapping <- unlist(profile$mapping)
  legacy <- any(c("hm_chrom", "hm_pos", "hm_effect_allele", "hm_other_allele") %in% names(sample))
  if (legacy) {
    mapping[c("chr", "pos", "ea", "nea", "b", "eaf", "marker")] <-
      c("hm_chrom", "hm_pos", "hm_effect_allele", "hm_other_allele", "hm_beta",
        "hm_effect_allele_frequency", "hm_rsid")
  }
  if (!all(mapping[.gwas_required_fields] %in% names(sample))) return(profile)
  mapping <- mapping[mapping %in% names(sample)]
  n_columns <- .gact_gwas_n_columns(names(sample))
  if (length(n_columns) == 1L) mapping["n"] <- n_columns
  if (length(n_columns) > 1L) profile$issues <- "Ambiguous per-variant n columns; review the sample-size mapping"
  profile$mapping <- as.list(mapping)
  profile$headers <- names(sample)
  profile$evidence <- paste("Catalog harmonised GRCh38 metadata and explicit",
    if (legacy) "hm_* harmonised" else "standard", "column definitions; beta and SE on signed effect scale")
  # Extra ratio columns are acceptable only when they agree with the mapped beta.
  beta <- suppressWarnings(as.numeric(sample[[mapping[["b"]]]]))
  ratios <- if (legacy) "hm_odds_ratio" else c("odds_ratio", "hazard_ratio")
  for (column in intersect(ratios, names(sample))) {
    ratio <- suppressWarnings(as.numeric(sample[[column]]))
    valid <- is.finite(beta) & is.finite(ratio) & ratio > 0
    if (any(abs(beta[valid] - log(ratio[valid])) > 5e-4 * (1 + abs(beta[valid]))))
      profile$issues <- c(profile$issues, paste("Mapped beta disagrees with log", column))
  }
  pos <- suppressWarnings(as.numeric(sample[[mapping[["pos"]]]]))
  se <- suppressWarnings(as.numeric(sample[[mapping[["seb"]]]]))
  p <- suppressWarnings(as.numeric(sample[[mapping[["p"]]]]))
  usable <- is.finite(pos) & pos > 0 & is.finite(beta) & is.finite(se) & se > 0 &
    is.finite(p) & p >= 0 & p <= 1 &
    toupper(sample[[mapping[["ea"]]]]) %in% c("A", "C", "G", "T") &
    toupper(sample[[mapping[["nea"]]]]) %in% c("A", "C", "G", "T")
  if (!any(usable)) profile$issues <- c(profile$issues, "No usable harmonised variants in the bounded preview")
  profile
}

.gact_panukb_original_row <- function(row, ancestry) {
  if (!ancestry %in% .gact_panukb_ancestries ||
      !ancestry %in% trimws(strsplit(row$pops, ",", fixed = TRUE)[[1L]]))
    stop("Selected ancestry is absent from the Pan-UKB phenotype")
  required <- c("filename", "aws_path", "md5_hex", "pheno_sex", "trait_type")
  missing <- required[vapply(required, function(field) {
    value <- row[[field]]
    length(value) != 1L || is.na(value) || !nzchar(trimws(value))
  }, logical(1))]
  if (length(missing))
    stop("Missing original Pan-UKB metadata: ", paste(missing, collapse = ", "))
  if (!isTRUE(grepl("^[[:xdigit:]]{32}$", row$md5_hex)))
    stop("Invalid original Pan-UKB provider MD5")
  trait <- row$description
  fallback <- length(trait) != 1L || is.na(trait) || !nzchar(trimws(trait))
  if (fallback) {
    trait <- row$phenocode
    if (row$trait_type != "prescriptions" || length(trait) != 1L ||
        is.na(trait) || !nzchar(trimws(trait)))
      stop("Missing original Pan-UKB metadata: description")
  }
  filename <- trimws(row$filename)
  if (!identical(basename(filename), filename)) stop("Invalid Pan-UKB filename")
  url <- sub("^s3://pan-ukb-us-east-1/", "https://pan-ukb-us-east-1.s3.amazonaws.com/",
             row$aws_path)
  if (!startsWith(url, "https://pan-ukb-us-east-1.s3.amazonaws.com/") ||
      !identical(basename(url), filename)) stop("Unexpected original Pan-UKB download URL")
  url <- paste0(dirname(url), "/", utils::URLencode(filename, reserved = TRUE))
  quantitative <- row$trait_type %in% c("continuous", "biomarkers")
  if (!row$trait_type %in% c("continuous", "biomarkers", "prescriptions", "icd10", "phecode", "categorical"))
    stop("Unknown Pan-UKB trait type")
  cases <- suppressWarnings(as.numeric(row[[paste0("n_cases_", ancestry)]]))
  controls <- if (quantitative) 0 else
    suppressWarnings(as.numeric(row[[paste0("n_controls_", ancestry)]]))
  if (length(cases) != 1L || length(controls) != 1L ||
      !all(is.finite(c(cases, controls))) || cases <= 0 || controls < 0 ||
      any(c(cases, controls) %% 1 != 0) || (!quantitative && controls == 0))
    stop("Missing or invalid ancestry-specific sample counts")
  gender <- unname(c(both_sexes = "both", females = "females", males = "males")[row$pheno_sex])
  if (is.na(gender)) stop("Unknown Pan-UKB sex")
  id <- toupper(substr(as.character(openssl::sha256(filename)), 1L, 24L))
  cache_filename <- filename
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_. -]*\\.tsv\\.bgz$", filename))
    cache_filename <- paste0("panukb-", tolower(id), ".tsv.bgz")
  comments <- paste("Original Pan-UKB", filename, "ancestry", ancestry)
  if (fallback) comments <- paste(comments,
    "(trait label from provider phenocode, description missing)")
  data.frame(accession = paste0("PANUKB-", id, "-", ancestry), filename = cache_filename,
    url = url, md5 = tolower(row$md5_hex), source = paste("Pan-UK Biobank original", ancestry),
    trait = trait, type = if (quantitative) "quantitative" else "binary",
    build = "GRCh37", ancestry = ancestry, gender = gender, strand = "forward",
    category = row$category %||% NA_character_, study_n_total = cases + controls,
    study_n_case = if (quantitative) NA_real_ else cases,
    study_n_control = if (quantitative) NA_real_ else controls,
    cohort_evidence = "Pan-UKB ancestry-specific GWAS columns and phenotype manifest counts",
    reference = "https://pan.ukbb.broadinstitute.org/docs/per-phenotype-files/",
    comments = comments, coverage = "genome_wide",
    coverage_evidence = "Pan-UKB official genome-wide per-phenotype SNP release")
}

.gact_panukb_original_profile <- function(ancestry) {
  fields <- c(chr = "chr", pos = "pos", ea = "alt", nea = "ref",
    b = paste0("beta_", ancestry), seb = paste0("se_", ancestry),
    p = paste0("neglog10_pval_", ancestry))
  list(source = paste("Pan-UK Biobank original", ancestry), mapping = as.list(fields),
    effect_scale = "beta", p_scale = "neglog10p", verified = TRUE,
    evidence = "Pan-UKB original per-phenotype definitions, ALT effect and forward-strand REF")
}

.gact_panukb_code_pattern <- "\\((phecode|ukb data field|icd-?10(?: code)?)[ :]+([a-z0-9.]+)\\)"

.gact_panukb_catalog_code <- function(label) {
  label <- tolower(label)
  parts <- regmatches(label, regexec(.gact_panukb_code_pattern, label, perl = TRUE))
  vapply(parts, function(part) {
    if (!length(part)) return(NA_character_)
    namespace <- if (part[2L] == "ukb data field") "ukb" else
      if (part[2L] == "phecode") "phecode" else "icd10"
    code <- if (namespace == "icd10") gsub(".", "", part[3L], fixed = TRUE) else part[3L]
    paste0(namespace, ":", code)
  }, character(1))
}

.gact_panukb_provider_code <- function(type, code) {
  if (!length(type) || !length(code)) return(character())
  namespace <- unname(c(biomarkers = "ukb", continuous = "ukb", categorical = "ukb",
                        phecode = "phecode", icd10 = "icd10")[type])
  code <- tolower(code)
  icd <- which(type == "icd10")
  code[icd] <- gsub(".", "", code[icd], fixed = TRUE)
  result <- paste0(namespace, ":", code)
  result[is.na(namespace) | is.na(code) | !nzchar(code)] <- NA_character_
  result
}

.gact_panukb_trait_key <- function(x) {
  x <- trimws(gsub(.gact_panukb_code_pattern, "", tolower(x), perl = TRUE))
  x <- sub(" levels$", "", x)
  gsub("[^a-z0-9]", "", gsub("haemoglobin", "hemoglobin", x, fixed = TRUE))
}

.gact_panukb_label_keys <- function(row) {
  labels <- row$description
  if (!is.null(row$coding_description) && !is.na(row$coding_description) &&
      nzchar(row$coding_description)) {
    labels <- c(labels, paste(row$description, row$coding_description))
  }
  .gact_panukb_trait_key(labels)
}

.gact_panukb_keys <- function(row) {
  code <- .gact_panukb_provider_code(row$trait_type, row$phenocode)
  unique(c(paste0("label:", .gact_panukb_label_keys(row)), code[!is.na(code)]))
}

.gact_panukb_trait_score <- function(row, labels) {
  code <- .gact_panukb_provider_code(row$trait_type, row$phenocode)
  codes <- unique(na.omit(.gact_panukb_catalog_code(labels)))
  if (length(codes) && (length(code) != 1L || is.na(code) || !all(codes == code))) return(0L)
  if (any(.gact_panukb_label_keys(row) %in% .gact_panukb_trait_key(labels))) return(2L)
  if (length(codes)) 1L else 0L
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

.gact_panukb_ancestries <- c("AFR", "AMR", "CSA", "EAS", "EUR", "MID")

.gact_catalog_root <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"

.gact_catalog_index <- function() {
  url <- paste0(.gact_catalog_root, "harmonised_list.txt")
  response <- curl::curl_fetch_memory(url, handle = curl::new_handle(timeout = 120))
  if (response$status_code != 200L) stop("Cannot read the GWAS Catalog harmonised-file index")
  .gact_parse_catalog_index(strsplit(rawToChar(response$content), "\n", fixed = TRUE)[[1L]])
}

.gact_parse_catalog_index <- function(lines) {
  paths <- unique(sub("^\\./", "", trimws(lines)))
  pattern <- "^GCST[0-9]{6,8}-GCST[0-9]{6,8}/(GCST[0-9]{6,8})/harmonised/[A-Za-z0-9._-]+\\.h\\.tsv\\.gz$"
  paths <- paths[grepl(pattern, paths)]
  if (!length(paths)) stop("No harmonised gzip TSV files in the GWAS Catalog index")
  accession <- sub(pattern, "\\1", paths)
  ambiguous <- duplicated(accession) | duplicated(accession, fromLast = TRUE)
  keep <- !duplicated(accession)
  data.frame(filename = basename(paths[keep]), analysis = "catalog", accession = accession[keep],
             catalog_path = paths[keep], ambiguous = ambiguous[keep])
}

.gact_catalog_file <- function(id, metadata_dir = NULL, catalog_path = NULL) {
  file <- paste0(id, ".h.tsv.gz")
  num <- as.integer(sub("^GCST", "", id))
  upper <- ceiling(num / 1000) * 1000
  width <- nchar(id) - 4L
  range <- sprintf("GCST%0*d-GCST%0*d", width, upper - 999L, width, upper)
  base <- paste0(.gact_catalog_root, range, "/", id, "/harmonised/")
  if (!is.null(catalog_path)) {
    file <- basename(catalog_path)
    base <- paste0(.gact_catalog_root, dirname(catalog_path), "/")
  }
  url <- paste0(base, file)
  fetch <- function(remote, local) {
    if (is.null(metadata_dir)) {
      response <- curl::curl_fetch_memory(remote, handle = curl::new_handle(timeout = 45))
      if (response$status_code != 200L) stop("Cannot fetch Catalog metadata: ", remote)
      text <- rawToChar(response$content)
      Encoding(text) <- "UTF-8"
      text
    } else paste(readLines(file.path(metadata_dir, local), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  }
  meta <- yaml::yaml.load(fetch(paste0(url, "-meta.yaml"), paste0(file, "-meta.yaml")), eval.expr = FALSE)
  list(file = file, base = base, url = url, meta = meta, fetch = fetch)
}

.gact_catalog_md5 <- function(asset, id) {
  listing <- strsplit(asset$fetch(paste0(asset$base, "md5sum.txt"), paste0(id, ".md5")),
                      "\n", fixed = TRUE)[[1L]]
  entries <- strsplit(trimws(listing), "[[:space:]]+")
  match <- Filter(function(x) length(x) == 2L && sub("^\\*", "", x[2L]) == asset$file, entries)
  checksums <- vapply(match, function(x) x[1L], character(1))
  checksum <- tolower(asset$meta$data_file_md5sum)
  if (length(checksum) != 1L || !grepl("^[0-9a-f]{32}$", checksum) ||
      (length(checksums) && !identical(tolower(checksums), checksum))) {
    stop("Catalog checksum metadata is invalid or disagrees with the checksum list: ", id)
  }
  checksum
}

.gact_panukb_catalog_row <- function(row, id, analysis, metadata_dir,
                                    catalog_label = NULL, alternatives = NULL) {
  asset <- .gact_catalog_file(id, metadata_dir)
  meta <- asset$meta
  file <- asset$file
  if (!identical(meta$gwas_id, id) || !identical(meta$data_file_name, file) ||
      !identical(meta$genome_assembly, "GRCh38") ||
      !identical(meta$coordinate_system, "1-based") || !isTRUE(meta$is_harmonised) ||
      !isTRUE(grepl("^GWAS-SSF", meta$file_type)) ||
      !any(grepl("pan.ukbb.broadinstitute.org", meta$author_notes, fixed = TRUE))) {
    .gact_panukb_mismatch("Catalog metadata does not match the selected Pan-UKB phenotype: ", id)
  }
  labels <- c(unlist(meta$trait_description, use.names = FALSE), catalog_label)
  score <- .gact_panukb_trait_score(row, labels)
  if (!score) .gact_panukb_mismatch("Catalog phenotype code or label differs: ", id)
  sexes <- c(both_sexes = "combined", females = "female", males = "male")
  if (!identical(meta$sex, unname(sexes[row$pheno_sex]))) .gact_panukb_mismatch("Catalog sex scope differs: ", id)
  pops <- switch(analysis, hq = strsplit(row$pops_pass_qc, ",", fixed = TRUE)[[1L]],
                 all = strsplit(row$pops, ",", fixed = TRUE)[[1L]], analysis)
  catalog_pops <- vapply(meta$samples, .gact_panukb_pop_code, character(1))
  if (!length(pops) || anyNA(pops) || anyDuplicated(catalog_pops) ||
      !setequal(pops, catalog_pops)) .gact_panukb_mismatch("Catalog ancestry scope differs: ", id)
  sizes <- vapply(meta$samples, function(sample) as.numeric(sample$sample_size), numeric(1))
  flags <- lapply(meta$samples, `[[`, "case_control_study")
  counts <- .gact_panukb_cohort_counts(row, catalog_pops, sizes, flags)
  if (is.null(counts)) {
    .gact_panukb_mismatch("Catalog cohort counts differ from the selected Pan-UKB analysis: ", id)
  }
  coverage <- .gact_require_gwas_coverage(meta,
    "Pan-UKB official genome-wide per-phenotype SNP release", id)
  for (i in seq_len(nrow(alternatives) %||% 0L)) {
    other <- alternatives[i, , drop = FALSE]
    if (.gact_panukb_trait_score(other, labels) >= score &&
        !is.null(.gact_panukb_cohort_counts(other, catalog_pops, sizes, flags))) {
      .gact_panukb_mismatch("Catalog phenotype, sex and cohort also match ", other$filename, ": ", id)
    }
  }
  binary <- counts$binary
  data.frame(accession = id, filename = file, url = asset$url, md5 = .gact_catalog_md5(asset, id),
    source = "Pan-UK Biobank Catalog harmonised GWAS-SSF v1.0", trait = row$description,
    type = if (binary) "binary" else "quantitative", build = "GRCh38",
    ancestry = if (length(pops) == 1L) pops else "Mixed",
    gender = switch(row$pheno_sex, both_sexes = "both", females = "females", males = "males"),
    category = row$category, strand = "forward", study_n_total = sum(sizes),
    study_n_case = if (binary) sum(counts$cases) else NA_real_,
    study_n_control = if (binary) sum(counts$controls) else NA_real_,
    comments = paste0("Pan-UKB ", row$filename, "; analysis=", analysis, "; groups=", paste(pops, collapse = ",")),
    cohort_evidence = if (analysis %in% .gact_panukb_ancestries)
      "Pan-UKB separate ancestry GWAS and matching Catalog sample metadata" else NA_character_,
    reference = .gact_gwas_catalog_url(id), coverage = coverage$coverage,
    coverage_evidence = coverage$evidence, stringsAsFactors = FALSE)
}

.gact_panukb_cohort_counts <- function(row, pops, sizes, flags) {
  if (!is.null(row$pops) &&
      !all(pops %in% trimws(strsplit(row$pops, ",", fixed = TRUE)[[1L]]))) return(NULL)
  binary <- !row$trait_type %in% c("continuous", "biomarkers")
  cases <- as.numeric(unlist(row[paste0("n_cases_", pops)]))
  controls <- if (binary) as.numeric(unlist(row[paste0("n_controls_", pops)])) else
    rep(0, length(pops))
  if (length(cases) != length(pops) || length(controls) != length(pops) ||
      anyNA(c(cases, controls, sizes)) || !all(sizes == cases + controls) ||
      !all(vapply(flags, identical, logical(1), binary))) return(NULL)
  list(binary = binary, cases = cases, controls = controls)
}

.gact_gwas_coverage <- function(meta, evidence = NULL) {
  notes <- tolower(paste(unlist(meta[c("author_notes", "file_type")]), collapse = " "))
  technology <- tolower(paste(unlist(meta$genotyping_technology), collapse = " "))
  restricted <- c(
    cis_only = "\\bcis[-_ ]only\\b|\\bonly (cis|variants within)\\b|\\bcis[-_ ](pqtl|eqtl|sqtl|mqtl) summary",
    targeted = "\\b(targeted|candidate)[-_ ](region|gene|variant|snp)|\\bregional[-_ ]only\\b",
    exome = "\\b(exome[-_ ]only|exome[-_ ]wide|whole[-_ ]exome)\\b",
    gene_based = "\\bgene[-_ ](based|level)\\b|\\bburden (test|results|analysis)\\b",
    hits_only = "\\b(significant|lead|top)[-_ ](variants|snps|hits)[-_ ]only\\b"
  )
  matched <- names(restricted)[vapply(restricted, function(pattern)
    grepl(pattern, notes, perl = TRUE), logical(1L))]
  if (!length(matched) && grepl("\\bexome\\b|exomechip", technology, perl = TRUE) &&
      !grepl("genome[-_ ]wide|whole[-_ ]genome", technology, perl = TRUE))
    matched <- "exome"
  if (length(matched)) return(list(coverage = "restricted",
    evidence = paste("Metadata indicates", paste(matched, collapse = ", "))))
  explicit <- grepl("\\bgenome[-_ ]wide[ -](association|summary|snp|variant|pqtl|mqtl|scan|analysis|results)\\b",
                    notes, perl = TRUE)
  negated <- grepl("\\b(not(?: a)?|non)[-_ ]genome[-_ ]wide\\b", notes, perl = TRUE)
  if (negated) return(list(coverage = "restricted", evidence = "Metadata indicates non-genome-wide results"))
  if (explicit) return(list(coverage = "genome_wide", evidence = paste(
    "Provider metadata", paste(unlist(meta[c("author_notes", "file_type")]), collapse = " "))))
  if (!is.null(evidence)) return(list(coverage = "genome_wide", evidence = evidence))
  list(coverage = "unclear", evidence = "Metadata does not establish genome-wide SNP coverage")
}

.gact_require_gwas_coverage <- function(meta, evidence, id) {
  coverage <- .gact_gwas_coverage(meta, evidence)
  if (coverage$coverage != "genome_wide")
    stop(errorCondition(paste(coverage$evidence, id), class = "gact_coverage_review",
                        coverage = coverage$coverage, coverage_evidence = coverage$evidence))
  coverage
}

.gact_catalog_coverage <- function(asset, evidence, id) {
  coverage <- .gact_gwas_coverage(asset$meta, evidence)
  if (coverage$coverage == "unclear") {
    technology <- paste(unlist(asset$meta$genotyping_technology), collapse = " ")
    # The small provider log supplies file-level evidence, not just an assay label.
    if (grepl("genome[-_ ]wide|whole[-_ ]genome", technology, ignore.case = TRUE)) {
      url <- paste0(asset$base, id, ".running.log")
      log <- tryCatch(asset$fetch(url, paste0(id, ".running.log")), error = function(e) "")
      lines <- strsplit(paste(log, collapse = "\n"), "\n", fixed = TRUE)[[1L]]
      line <- grep("sites successfully harmoni[sz]ed", lines, value = TRUE)
      pattern <- "\\([[:space:]]*([0-9,]+)[[:space:]]+of[[:space:]]+[0-9,]+[[:space:]]*\\)"
      matches <- regmatches(line, regexec(pattern, line))
      counts <- vapply(matches, function(x) if (length(x) > 1L)
        as.numeric(gsub(",", "", x[2L], fixed = TRUE)) else NA_real_, numeric(1))
      if (length(counts) == 1L && is.finite(counts) && counts >= 1000000) {
        coverage <- list(coverage = "genome_wide", evidence = paste0(
          "Genome-wide assay and provider harmonisation log: ", counts,
          " successfully harmonised variants; ", url,
          " (coverage screening, not a full-file coverage audit)"))
      }
    }
  }
  if (coverage$coverage != "genome_wide")
    stop(errorCondition(paste(coverage$evidence, id), class = "gact_coverage_review",
                        coverage = coverage$coverage, coverage_evidence = coverage$evidence))
  coverage
}

.gact_catalog_study_scope <- function(meta, id) {
  design <- paste(unlist(meta[c("analysis_software", "author_notes", "trait_description")]),
                  collapse = " ")
  if (grepl(paste0("\\b(METAL|METASOFT|GWAMA|MR[-_ ]MEGA|MTAG|METACARPA|",
                   "meta[-_ ]?analys[ie]s|meta[-_ ]?analy[sz](ed|ing))\\b"),
            design, ignore.case = TRUE, perl = TRUE))
    stop(errorCondition(paste("Catalog metadata indicates a meta-analysis:", id),
                        class = "gact_scan_excluded"))
  samples <- meta$samples
  ancestry <- unique(trimws(unlist(lapply(samples, `[[`, "sample_ancestry_category"),
                                  use.names = FALSE)))
  ancestry <- unique(trimws(unlist(strsplit(ancestry, ",|;| and | or ", perl = TRUE))))
  if (!length(ancestry) || anyNA(ancestry) || any(!nzchar(ancestry)) ||
      any(grepl("^(unknown|other|not reported|not available)( ancestry| ancestries)?$",
                ancestry, ignore.case = TRUE)))
    stop("Catalog metadata does not establish a single ancestry: ", id)
  if (any(grepl("^(mixed|multiple|multi[- ]ancestry|trans[- ]ancestry)( ancestry| ancestries)?$",
                ancestry, ignore.case = TRUE)))
    stop(errorCondition(paste("Catalog file contains multiple ancestries:", id),
                        class = "gact_scan_excluded"))
  codes <- vapply(ancestry, function(x) tryCatch(
    .gact_panukb_pop_code(list(sample_ancestry_category = x)), error = function(e)
      if (grepl("African American|Afro-Caribbean", x, ignore.case = TRUE)) "AFR" else x), "")
  if (length(unique(codes)) != 1L)
    stop(errorCondition(paste("Catalog file contains multiple ancestries:", id),
                        class = "gact_scan_excluded"))
  if (length(samples) != 1L) stop("Multiple Catalog sample groups need study-design review: ", id)
  list(ancestry = unname(codes[1L]), design = design)
}

.gact_catalog_recipe_row <- function(id, evidence, catalog_path = NULL, coverage_evidence = NULL) {
  asset <- .gact_catalog_file(id, catalog_path = catalog_path)
  meta <- asset$meta
  scope <- .gact_catalog_study_scope(meta, id)
  if (!identical(meta$gwas_id, id) || !identical(meta$data_file_name, asset$file) ||
      !identical(meta$genome_assembly, "GRCh38") || !identical(meta$coordinate_system, "1-based") ||
      !isTRUE(meta$is_harmonised)) {
    stop("Expected a harmonised, 1-based GRCh38 file: ", id)
  }
  coverage <- .gact_catalog_coverage(asset, coverage_evidence, id)
  trait <- unlist(meta$trait_description, use.names = FALSE)
  samples <- meta$samples
  if (is.null(evidence) && grepl("\\bsingle[- ]cohort\\b", scope$design, ignore.case = TRUE))
    evidence <- scope$design
  flags <- lapply(samples, `[[`, "case_control_study")
  if (!length(flags) || !all(vapply(flags, function(x) identical(x, TRUE) || identical(x, FALSE),
                                   logical(1))) || length(unique(flags)) != 1L) {
    stop("Catalog metadata does not establish the study type: ", id)
  }
  ancestry <- unlist(samples[[1L]]$sample_ancestry_category, use.names = FALSE)
  if (length(trait) != 1L || anyNA(trait) || !nzchar(trait))
    stop("Catalog metadata does not establish a single trait: ", id)
  size <- if (length(samples) == 1L) samples[[1L]]$sample_size else NULL
  counts <- stats::setNames(rep(NA_real_, 2L), c("case_count", "control_count"))
  if (flags[[1L]]) {
    for (field in names(counts)) {
      value <- samples[[1L]][[field]]
      if (is.null(value)) next
      if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
          value <= 0 || value %% 1 != 0 ||
          (is.numeric(size) && length(size) == 1L && is.finite(size) && value > size))
        stop("Invalid Catalog ", field, ": ", id)
      counts[field] <- value
    }
    if (all(is.finite(counts)) && is.numeric(size) && length(size) == 1L &&
        is.finite(size) && sum(counts) != size)
      stop("Catalog case/control counts disagree with sample_size: ", id)
  }
  data.frame(accession = id, filename = asset$file, url = asset$url,
    md5 = .gact_catalog_md5(asset, id), source = paste("GWAS Catalog harmonised", meta$file_type %||% "legacy"),
    trait = trait, type = if (flags[[1L]]) "binary" else "quantitative", build = "GRCh38",
    ancestry = scope$ancestry, strand = "forward", cohort_evidence = evidence %||% NA_character_,
    study_design = "not_reported_as_meta_analysis",
    gender = unname(c(combined = "both", female = "females", male = "males")[meta$sex %||% ""]),
    study_n_total = if (is.numeric(size) && length(size) == 1L && is.finite(size) && size > 0)
      size else NA_real_,
    study_n_case = unname(counts["case_count"]),
    study_n_control = unname(counts["control_count"]),
    comments = paste("Catalog ancestry:", paste(ancestry, collapse = ", ")),
    reference = .gact_gwas_catalog_url(id), coverage = coverage$coverage,
    coverage_evidence = coverage$evidence, stringsAsFactors = FALSE)
}

.gact_panukb_mismatch <- function(...) {
  stop(errorCondition(paste0(...), class = "gact_panukb_mismatch"))
}

#' Scan Pan-UK Biobank and GWAS Catalog studies for ingestion
#'
#' Match phenotypes or Catalog accessions to verified GRCh38 files. Use the
#' C++ reader on a bounded gzip prefix and prepare recipes and an import manifest.
#' Verified recipes are keyed by accession and input MD5. Use the same
#' `recipe_dir` to reuse them on later scans; metadata is checked again.
#' Sample counts are compared numerically so JSON integer and double
#' representations of the same count reuse the same recipe.
#' Pan-UKB discovery uses phenotype codes and labels, then checks sex,
#' ancestry and cohort counts against the Catalog metadata. Ambiguous
#' phenotype matches and shared checksums require review.
#' Unresolved studies appear in the review report. The scanner only writes
#' recipes and reports. The importer verifies the full file MD5.
#' On completion, recipes.tsv in recipe_dir is refreshed from all top-level
#' verified JSON recipes, including those outside the current selection. It
#' contains one metadata row per recipe. JSONs are not modified. Unverified
#' drafts are omitted; unreadable or malformed JSONs are skipped with a warning.
#' The TSV is replaced only after the new file has been fully written. It is
#' an inventory of cached recipes; the scan manifest records current readiness.
#' Explicit Catalog case_count/control_count fields are retained as study
#' metadata and checked against sample_size. Missing breakdowns are not inferred.
#' Original Pan-UKB scans use the provider manifest's download paths, MD5s,
#' sex and ancestry-specific counts. Recipes select one ancestry's signed beta,
#' SE and negative log10 p-value columns. Cross-ancestry meta-analysis columns
#' remain unmapped. Binary case/control frequencies remain unmapped.
#' Prescription traits with no description use the provider phenocode as the
#' label and record that choice. Reserved filename characters are URL-encoded
#' and use a stable safe local cache filename.
#' Catalog scans accept standard and legacy harmonised columns after preview
#' validation; a GWAS-SSF label is not required. Legacy hm_* coordinates and
#' alleles take precedence over source columns. Conflicting beta/ratio effects
#' and inferred or missing mappings require review.
#' Coverage is screened before previews and recipe reuse. Explicit restricted
#' scope is excluded. If scope wording is absent, a genome-wide assay plus at
#' least one million successfully harmonised variants in the provider's running
#' log supplies screening evidence. Neither an assay label nor a leading preview
#' alone suffices. This is not a full-file coverage audit. Only one ancestry is
#' accepted; reported meta-analyses (including MTAG) are excluded. Absence of a
#' meta-analysis statement is recorded as screening, not proof of a single cohort.
#' A unique standard sample-size column (n, N, N_sample_size, sample_size,
#' n_total or total_sample_size) is mapped to per-variant n. Missing n stays
#' missing; observed values must be finite and positive. Multiple candidates
#' require review. Study-level sample counts are not assigned to variants.
#'
#' @param selection Pan-UKB data frame with filename, an ancestry in the
#'   analysis column, and optional accession, or a vector of Catalog accessions.
#'   Omit it to select all available phenotype and ancestry pairs.
#' @param out_dir New or empty scan report directory. Defaults to a new
#'   directory under `recipe_dir/reviews`.
#' @param phenotype_manifest Official compressed Pan-UKB manifest URL or local path.
#' @param preview_rows Number of complete leading rows to inspect.
#' @param max_preview_bytes Maximum compressed bytes requested per file.
#' @param recipe_dir Persistent recipe directory. Defaults to `out_dir/recipes`.
#' @param source `"panukb"` or `"catalog"`.
#' @param analysis `"per_ancestry"` (default) or one Pan-UKB ancestry code
#'   for scans of the full phenotype manifest.
#' @param cohort_evidence Optional named character vector of reviewed
#'   study-design evidence keyed by Catalog accession, retained in the recipe.
#'   Literal single-cohort wording is not required. Evidence cannot override
#'   reported meta-analyses or multiple ancestries.
#' @param coverage_evidence Optional named character vector of reviewed evidence
#'   establishing genome-wide SNP coverage, keyed by Catalog accession.
#'   Explicitly restricted metadata takes precedence. Unclear coverage needs review.
#'   Proteins and metabolites are included when their files meet the same coverage rules.
#' @param GAlist Existing GRCh38 database descriptor. Its recipes directory is
#'   the default when supplied.
#' @param repository Provider name: `"pan-ukb"` or `"gwas-catalog"`.
#'   Catalog scans with no selection use the published harmonised-file index.
#' @param workers Number of concurrent recipe scans. Defaults to 1.
#' @param files For Pan-UKB, use `"catalog"` (default) for harmonised Catalog
#'   matches or `"original"` for the provider's GRCh37 files. Original files
#'   require GAlist and its installed liftover map. Each recipe selects one
#'   ancestry and records the preview build check and liftover index checksum.
#' @return Paths to the manifest and format profiles, and a per-study review table.
#'   Coverage is `genome_wide`, `restricted`, or `unclear`. Restricted files
#'   have status `excluded` and unclear files have status `needs_review`.
#'   The import manifest contains only ready studies.
#'   `recipes_tsv` is the path to the refreshed inventory in the recipe directory.
#' @export
scanGWAS <- function(selection = NULL, out_dir = NULL,
    phenotype_manifest = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz",
    preview_rows = 1000L, max_preview_bytes = 1024L * 1024L, recipe_dir = NULL,
    source = c("panukb", "catalog"), analysis = "per_ancestry", cohort_evidence = NULL,
    GAlist = NULL, repository = NULL, workers = 1L,
    files = c("catalog", "original"), coverage_evidence = NULL) {
  if (!is.null(repository)) {
    provider <- c("pan-ukb" = "panukb", "gwas-catalog" = "catalog")
    selected <- unname(provider[match.arg(repository, names(provider))])
    if (!missing(source) && !identical(match.arg(source), selected))
      stop("source and repository specify different providers")
    source <- selected
  }
  source <- match.arg(source)
  files <- match.arg(files)
  original <- files == "original"
  if (original && (source != "panukb" || is.null(GAlist)))
    stop("Original Pan-UKB scans require repository = 'pan-ukb' and GAlist")
  map_dir <- map_info <- NULL
  if (!is.numeric(workers) || length(workers) != 1L || !is.finite(workers) ||
      workers < 1L || workers > .Machine$integer.max || workers != floor(workers))
    stop("workers must be a positive integer")
  if (!is.null(GAlist)) {
    if (!identical(.gact_catalog_build(GAlist), "GRCh38"))
      stop("scanGWAS() requires a GRCh38 database when GAlist is supplied")
    if (is.null(recipe_dir)) {
      if (!"recipes" %in% names(GAlist$dirs)) stop("GAlist has no recipes directory")
      recipe_dir <- GAlist$dirs[["recipes"]]
    }
    if (original) {
      map_dir <- .gact_liftover_dir(GAlist)
      map_info <- .gact_liftover_index(map_dir, GAlist)
    }
  }
  if (!original && !requireNamespace("yaml", quietly = TRUE)) stop("Install yaml to read Catalog metadata")
  if (is.null(out_dir) && is.null(recipe_dir)) stop("Supply GAlist, out_dir, or recipe_dir")
  if (is.null(out_dir)) out_dir <- tempfile("scan-", tmpdir = file.path(recipe_dir, "reviews"))
  if (is.null(recipe_dir)) recipe_dir <- file.path(out_dir, "recipes")
  if (!is.null(cohort_evidence) && (!is.character(cohort_evidence) ||
      is.null(names(cohort_evidence)) || anyNA(cohort_evidence) ||
      any(!nzchar(trimws(cohort_evidence))) || anyDuplicated(names(cohort_evidence)) ||
      !all(grepl("^GCST[0-9]{6,8}$", names(cohort_evidence))))) {
    stop("cohort_evidence must contain reviewed evidence named by Catalog accession")
  }
  cohort_evidence <- as.list(cohort_evidence)
  if (!is.null(coverage_evidence) && (!is.character(coverage_evidence) ||
      is.null(names(coverage_evidence)) || anyNA(coverage_evidence) ||
      any(!nzchar(trimws(coverage_evidence))) || anyDuplicated(names(coverage_evidence)) ||
      !all(grepl("^GCST[0-9]{6,8}$", names(coverage_evidence)))))
    stop("coverage_evidence must contain reviewed evidence named by Catalog accession")
  coverage_evidence <- as.list(coverage_evidence)
  if (source == "panukb") {
    phenotype <- .gact_read_panukb_manifest(phenotype_manifest)
    if (is.null(selection)) {
      if (length(analysis) != 1L || !analysis %in% c("per_ancestry", .gact_panukb_ancestries))
        stop("Choose per_ancestry or one Pan-UKB ancestry code")
      populations <- lapply(strsplit(phenotype$pops, ",", fixed = TRUE), function(pops) {
        eligible <- if (analysis == "per_ancestry") .gact_panukb_ancestries else analysis
        intersect(trimws(pops), eligible)
      })
      selection <- data.frame(filename = rep(phenotype$filename, lengths(populations)),
                              analysis = unlist(populations, use.names = FALSE))
    }
  } else {
    if (is.null(selection)) {
      selection <- .gact_catalog_index()
    } else {
      if (!is.character(selection) || !length(selection) || anyNA(selection) ||
          anyDuplicated(selection) || !all(grepl("^GCST[0-9]{6,8}$", selection))) {
        stop("Supply unique GWAS Catalog accessions")
      }
      selection <- data.frame(filename = paste0(selection, ".h.tsv.gz"),
                              analysis = "catalog", accession = selection)
    }
  }
  required <- c("filename", "analysis")
  if (!is.data.frame(selection) || !nrow(selection) || !all(required %in% names(selection)) ||
      anyNA(selection[required]) || (source == "panukb" &&
      (anyDuplicated(selection[required]) || !all(selection$analysis %in% .gact_panukb_ancestries)))) {
    stop("Supply unique Pan-UKB filename and single-ancestry analysis pairs")
  }
  if (!"accession" %in% names(selection)) selection$accession <- NA_character_
  selection$accession <- as.character(selection$accession)
  explicit <- !is.na(selection$accession) & nzchar(selection$accession)
  if (original && any(explicit)) stop("Original Pan-UKB scans use provider IDs. Omit Catalog accessions")
  if (!all(grepl("^GCST[0-9]{6,8}$", selection$accession[explicit]))) stop("Invalid Catalog accession")
  for (value in list(preview_rows, max_preview_bytes)) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value >= .Machine$integer.max || value != floor(value)) {
      stop("Preview limits must be positive integers")
    }
  }
  if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) {
    stop("out_dir must be new or empty")
  }
  at <- if (source == "panukb") match(selection$filename, phenotype$filename)
  catalog <- if (source == "panukb" && !original && any(!explicit & !is.na(at))) {
    tryCatch(.gact_panukb_index_catalog(.gact_panukb_catalog_studies()),
             error = function(e) e)
  } else list()
  phenotype_index <- if (source == "panukb") .gact_panukb_index_phenotypes(phenotype)
  for (path in c(out_dir, recipe_dir)) {
    if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE))
      stop("Cannot create scan directory: ", path)
  }
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  recipe_dir <- normalizePath(recipe_dir, winslash = "/", mustWork = TRUE)
  review <- data.frame(filename = selection$filename, analysis = selection$analysis,
    accession = selection$accession, status = "needs_review", preview_bytes = NA_integer_,
    preview_rows = NA_integer_, reused = FALSE, recipe = NA_character_, issues = NA_character_,
    coverage = NA_character_, coverage_evidence = NA_character_)
  rows <- formats <- vector("list", nrow(selection))
  context <- list(
    selection = selection, review = review, source = source, at = at, explicit = explicit,
    phenotype = if (source == "panukb") phenotype else NULL, catalog = catalog,
    phenotype_index = phenotype_index, cohort_evidence = cohort_evidence,
    coverage_evidence = coverage_evidence,
    out_dir = out_dir, recipe_dir = recipe_dir, preview_rows = preview_rows,
    max_preview_bytes = max_preview_bytes, original = original,
    map_dir = map_dir, map_info = map_info
  )
  count <- min(as.integer(workers), nrow(selection))
  if (count > 1L) {
    cluster <- parallel::makePSOCKcluster(count)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, .libPaths, .libPaths())
    .gact_scan_context <- context
    parallel::clusterExport(cluster, ".gact_scan_context", envir = environment())
  }
  for (start in seq.int(1L, nrow(selection), by = 100L)) {
    indices <- seq.int(start, min(start + 99L, nrow(selection)))
    message("Scanning GWAS recipes ", start, "–", tail(indices, 1L), "/", nrow(selection),
            " (", count, " workers)")
    results <- if (count > 1L) {
      parallel::parLapplyLB(cluster, indices, .gact_scan_recipe_worker, chunk.size = 1L)
    } else {
      lapply(indices, .gact_scan_recipe_one, context = context)
    }
    for (j in seq_along(indices)) {
      i <- indices[j]
      result <- results[[j]]
      review[i, ] <- result$review
      if (identical(result$review$status, "ready")) {
        rows[[i]] <- result$manifest
        formats[[i]] <- result$format
      }
    }
  }
  accepted <- which(review$status == "ready")
  for (field in c("accession", "md5")) {
    values <- vapply(rows[accepted], function(row) row[[field]], character(1))
    duplicate <- accepted[duplicated(values) | duplicated(values, fromLast = TRUE)]
    if (original && field == "md5" && length(duplicate)) {
      shared <- do.call(rbind, rows[duplicate])
      same_file <- vapply(split(shared$filename, shared$md5),
                          function(x) length(unique(x)) == 1L, logical(1L))
      allowed <- .gact_shared_panukb_inputs(shared) & same_file[shared$md5]
      duplicate <- duplicate[!allowed]
    }
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
  recipes_tsv <- .gact_write_recipe_inventory(recipe_dir)
  invisible(list(manifest = manifest, knowledge = knowledge, review = review, report = report,
                 recipes_tsv = recipes_tsv))
}

.gact_write_recipe_inventory <- function(recipe_dir) {
  fields <- c("accession", "filename", "source", "trait", "type", "ancestry",
    "build", "gender", "study_n_total", "study_n_case", "study_n_control",
    "reference", "input_url", "input_md5")
  numeric_fields <- c("study_n_total", "study_n_case", "study_n_control")
  skipped <- character()
  paths <- list.files(recipe_dir, pattern = "\\.json$", full.names = TRUE)
  rows <- lapply(paths, function(path) {
    tryCatch({
      recipe <- jsonlite::read_json(path, simplifyVector = TRUE)
      if (!is.list(recipe) || !isTRUE(recipe$verified)) return(NULL)
      row <- stats::setNames(lapply(fields, function(field) {
        value <- recipe[[field]]
        if (is.null(value) || !length(value)) return(if (field %in% numeric_fields) NA_real_ else NA_character_)
        if (!is.atomic(value) || length(value) != 1L) stop("Nonscalar recipe metadata")
        value
      }), fields)
      if (anyNA(unlist(row[c("accession", "input_md5")])) ||
          any(!nzchar(unlist(row[c("accession", "input_md5")]))))
        stop("Missing recipe identity")
      row$recipe_file <- basename(path)
      row
    }, error = function(e) {
      skipped <<- c(skipped, basename(path))
      NULL
    })
  })
  # Keep a stable header even when no verified recipes exist.
  empty <- stats::setNames(lapply(c(fields, "recipe_file"), function(field)
    if (field %in% numeric_fields) numeric() else character()), c(fields, "recipe_file"))
  table <- data.table::rbindlist(c(list(empty), Filter(Negate(is.null), rows)), fill = TRUE)
  if (nrow(table)) table <- table[order(table$accession, table$recipe_file, method = "radix")]
  output <- file.path(recipe_dir, "recipes.tsv")
  staged <- tempfile(".recipes-", tmpdir = recipe_dir, fileext = ".tsv")
  on.exit(unlink(staged), add = TRUE)
  data.table::fwrite(table, staged, sep = "\t", na = "NA")
  .publish_gwas_import_file(staged, output, replace = TRUE)
  if (length(skipped)) warning("Skipped ", length(skipped),
    " unreadable or malformed recipe JSON file(s) while updating recipes.tsv: ",
    paste(head(skipped, 5L), collapse = ", "), call. = FALSE)
  output
}

.gact_scan_recipe_worker <- function(i) {
  .gact_scan_recipe_one(i, get(".gact_scan_context", envir = .GlobalEnv))
}

.gact_scan_recipe_one <- function(i, context) {
  with(context, {
    entry <- review[i, , drop = FALSE]
    result <- tryCatch({
      picked <- if (original) {
        if (is.na(at[i])) stop("Filename is absent from the Pan-UKB phenotype manifest")
        .gact_panukb_original_row(phenotype[at[i], , drop = FALSE], selection$analysis[i])
      } else if (source == "catalog") {
        id <- selection$accession[i]
        if (isTRUE(selection$ambiguous[i])) stop("Multiple harmonised files listed for ", id)
        .gact_catalog_recipe_row(id, cohort_evidence[[id]], selection$catalog_path[i],
                                 coverage_evidence[[id]])
      } else {
        if (is.na(at[i])) stop("Filename is absent from the Pan-UKB phenotype manifest")
        row <- phenotype[at[i], , drop = FALSE]
        candidates <- if (explicit[i]) {
          data.frame(accession = selection$accession[i], label = NA_character_)
        } else {
          if (inherits(catalog, "error")) stop(catalog)
          .gact_panukb_candidates(row, catalog, selection$analysis[i])
        }
        alternatives <- if (!explicit[i]) {
          siblings <- unlist(phenotype_index[paste(.gact_panukb_keys(row), row$pheno_sex)],
                             use.names = FALSE)
          phenotype[setdiff(siblings, at[i]), , drop = FALSE]
        }
        .gact_panukb_scan_match(row, candidates, selection$analysis[i], alternatives)
      }
      if (source == "panukb" && is.null(picked$coverage)) {
        picked$coverage <- "genome_wide"
        picked$coverage_evidence <- "Pan-UKB documented genome-wide per-phenotype SNP release"
      }
      entry$accession <- picked$accession
      entry$coverage <- picked$coverage
      entry$coverage_evidence <- picked$coverage_evidence
      path <- .gact_verified_recipe_path(recipe_dir, picked$accession, picked$md5)
      stored <- .gact_read_verified_recipe(path, picked)
      profile <- .gact_panukb_catalog_profile(picked$source)
      if (original) profile <- .gact_panukb_original_profile(picked$ancestry)
      if (is.null(stored)) {
        prefix <- tempfile(fileext = ".tsv.gz")
        sample <- tryCatch({
          bytes <- .gact_fetch_gwas_prefix(picked$url, prefix, max_preview_bytes)
          entry$preview_bytes <- bytes
          .gact_cpp_gwas_preview(prefix, preview_rows, 0L, prefix = TRUE)
        }, finally = unlink(prefix))
        entry$preview_rows <- nrow(sample)
        if (nrow(sample) < preview_rows)
          stop("Byte limit reached before the requested number of complete rows")
        if (source == "catalog") profile <- .gact_catalog_preview_profile(sample, picked$source)
        standard <- all(unlist(profile$mapping) %in% names(sample)) &&
          (source == "catalog" ||
           !any(c("odds_ratio", "hazard_ratio", "neg_log_10_p_value") %in% names(sample)))
        known_panukb <- all(names(sample) %in% c(profile$headers, "low_confidence"))
        if (original) {
          eaf <- paste0("af_", picked$ancestry)
          if (picked$type == "quantitative" && eaf %in% names(sample))
            profile$mapping$eaf <- eaf
        }
        if (original || (source != "catalog" && standard && known_panukb)) {
          profile$headers <- names(sample)
        }
        recipe <- .gact_gwas_scan_recipe(picked, sample, bytes, profile)
        if (original) {
          recipe$build_check <- verifyGWASBuild(sample, map_dir)
          recipe$build_check$index_md5 <- map_info$md5
          recipe$liftover <- list(source_build = "GRCh37", build = "GRCh38",
                                  index_md5 = map_info$md5)
          if (!identical(recipe$build_check$build, "GRCh37"))
            recipe$inspection_issues <- c(recipe$inspection_issues,
              "Preview does not verify the original Pan-UKB GRCh37 build")
        }
        if (length(recipe$inspection_issues)) {
          draft <- file.path(out_dir, sprintf("study_%05d.json", i))
          jsonlite::write_json(recipe, draft, auto_unbox = TRUE, pretty = TRUE,
                               na = "null", digits = NA)
          entry$recipe <- draft
          stop(paste(recipe$inspection_issues, collapse = " | "))
        }
        recipe$verified <- TRUE
        .gact_save_verified_recipe(recipe, path, picked)
      } else {
        if (original && (!identical(stored$build_check$build, "GRCh37") ||
                         !identical(stored$build_check$index_md5, map_info$md5)))
          stop("Saved recipe build check differs from the installed liftover index. Review the recipe")
        entry$reused <- TRUE
        entry$preview_bytes <- 0L
        entry$preview_rows <- stored$inspection$preview_rows
        profile$headers <- unlist(stored$inspection$headers, use.names = FALSE)
        profile$mapping <- stored$mapping
        profile$effect_scale <- stored$effect_scale
        profile$p_scale <- stored$p_scale
      }
      entry$recipe <- path
      list(manifest = picked, format = profile)
    }, error = function(e) e)
    if (inherits(result, "error")) {
      if (inherits(result, "gact_coverage_review")) {
        entry$coverage <- result$coverage
        entry$coverage_evidence <- result$coverage_evidence
        if (result$coverage == "restricted") entry$status <- "excluded"
      }
      if (inherits(result, "gact_scan_excluded")) entry$status <- "excluded"
      entry$issues <- conditionMessage(result)
      return(list(review = entry))
    }
    entry$status <- "ready"
    c(list(review = entry), result)
  })
}

.gact_verified_recipe_path <- function(directory, accession, md5) {
  md5 <- tolower(md5)
  if (!isTRUE(grepl("^[A-Z0-9][A-Z0-9._-]*$", accession)) ||
      !isTRUE(grepl("^[a-f0-9]{32}$", md5))) {
    stop("Supply a valid accession and the provider's MD5 checksum")
  }
  file.path(directory, paste0(accession, "-", md5, ".json"))
}

.gact_shared_panukb_inputs <- function(manifest) {
  result <- rep(FALSE, nrow(manifest))
  groups <- split(seq_len(nrow(manifest)), .gact_gwas_path_key(manifest$filename))
  for (rows in groups[lengths(groups) > 1L]) {
    selected <- manifest[rows, , drop = FALSE]
    if (!all(c("source", "ancestry") %in% names(selected))) next
    valid <- !anyNA(selected[c("url", "md5", "source", "ancestry")]) &&
      length(unique(tolower(selected$md5))) == 1L && length(unique(selected$url)) == 1L &&
      all(selected$ancestry %in% .gact_panukb_ancestries) && !anyDuplicated(selected$ancestry) &&
      all(selected$source == paste("Pan-UK Biobank original", selected$ancestry))
    if (valid) result[rows] <- TRUE
  }
  result
}

.gact_read_verified_recipe <- function(path, study) {
  if (!file.exists(path)) return(NULL)
  recipe <- jsonlite::read_json(path)
  if (!isTRUE(recipe$verified)) return(NULL)
  if (!identical(recipe$accession, toupper(study$accession)) ||
      !identical(recipe$input_md5, tolower(study$md5)) ||
      !identical(recipe$input_url, study$url) ||
      !isTRUE(recipe$inspection$profile_reused) || !isTRUE(recipe$validation$ok) ||
      length(recipe$inspection_issues)) {
    stop("Stored verified recipe does not match the input or validation: ", basename(path))
  }
  fields <- intersect(c("filename", "cohort_evidence", .gact_gwas_metadata_fields), names(study))
  same <- vapply(fields, function(field) {
    value <- study[[field]]
    if (is.na(value)) return(TRUE)
    if (field %in% c("n", "ncase", "ncontrol", "study_n_total", "study_n_case", "study_n_control")) {
      saved <- suppressWarnings(as.numeric(recipe[[field]]))
      expected <- suppressWarnings(as.numeric(value))
      return(length(saved) == 1L && is.finite(saved) && is.finite(expected) && saved == expected)
    }
    identical(as.character(recipe[[field]]), as.character(value))
  }, logical(1))
  if (!all(same)) stop("Stored verified recipe metadata changed: ", paste(fields[!same], collapse = ", "))
  if (!all(unlist(recipe$mapping) %in% unlist(recipe$inspection$headers)))
    stop("Stored recipe mapping differs from its preview header")
  .gact_cpp_gwas_mapping(recipe)
  recipe
}

.gact_save_verified_recipe <- function(recipe, path, study) {
  key <- substr(as.character(openssl::sha256(basename(path))), 1L, 16L)
  lock <- file.path(dirname(path), paste0(".recipe-", key))
  if (!dir.create(lock, showWarnings = FALSE)) stop("Cannot lock recipe: ", basename(path))
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (file.exists(path)) {
    saved <- .gact_read_verified_recipe(path, study)
    if (is.null(saved) || !identical(saved$mapping, recipe$mapping) ||
        !identical(unlist(saved$inspection$headers), recipe$inspection$headers)) {
      stop("Existing recipe differs from the verified preview: ", basename(path))
    }
    return(invisible(path))
  }
  staged <- file.path(lock, "recipe.json")
  jsonlite::write_json(recipe, staged, auto_unbox = TRUE, pretty = TRUE, na = "null", digits = NA)
  .publish_gwas_import_file(staged, path)
  invisible(path)
}

.gact_panukb_scan_match <- function(row, candidates, analysis, alternatives) {
  matches <- lapply(seq_len(nrow(candidates)), function(i) {
    tryCatch(.gact_panukb_catalog_row(row, candidates$accession[i], analysis, metadata_dir = NULL,
                                    catalog_label = candidates$label[i], alternatives = alternatives),
      gact_panukb_mismatch = function(e) e)
  })
  matched <- !vapply(matches, inherits, logical(1), "gact_panukb_mismatch")
  if (!any(matched)) {
    detail <- if (length(matches)) paste(unique(vapply(matches, conditionMessage, "")), collapse = " | ")
    stop("No verified Catalog match", if (length(detail)) paste0(": ", detail))
  }
  if (sum(matched) > 1L)
    stop("Multiple verified Catalog matches: ", paste(candidates$accession[matched], collapse = ", "))
  matches[[which(matched)]]
}

.gact_gwas_scan_recipe <- function(study, sample, bytes, profile) {
  recipe <- as.list(study)
  recipe$input_md5 <- study$md5
  recipe$input_url <- study$url
  reused <- identical(names(sample), profile$headers)
  recipe$inspection <- list(headers = names(sample), preview_rows = nrow(sample),
    preview_bytes = bytes, remote = TRUE, skip = 0L, profile_reused = reused,
    evidence = if (reused) profile$evidence else NULL,
    unmapped_columns = if (reused) setdiff(names(sample), unlist(profile$mapping)) else NULL)
  if (reused) {
    recipe$mapping <- profile$mapping
    recipe$effect_scale <- profile$effect_scale
    recipe$p_scale <- profile$p_scale
    recipe$validation <- validateStatSchema(normalizeStatSchema(sample, recipe$mapping,
      effect_scale = recipe$effect_scale, p_scale = recipe$p_scale))
    recipe$inspection_issues <- c(profile$issues, recipe$validation$errors)
    if (!length(recipe$inspection_issues)) {
      recipe$inspection_issues <- tryCatch({
        .gact_cpp_gwas_mapping(recipe)
        character()
      }, error = conditionMessage)
    }
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

.gact_panukb_index_catalog <- function(studies) {
  studies <- Filter(function(study) isTRUE(grepl("^GCST[0-9]{8}$", study$accessionId)) &&
    isTRUE(grepl("^https?://", study$summaryStatistics)), studies)
  records <- data.frame(
    accession = vapply(studies, `[[`, character(1), "accessionId"),
    label = vapply(studies, function(study) study$reportedTrait %||% "", character(1)),
    ancestry = vapply(studies, function(study) {
      groups <- unlist(study$discoverySampleAncestry, use.names = FALSE)
      if (!length(groups) || any(!grepl("^[0-9,]+ ", groups))) return(NA_character_)
      groups <- sub("^[0-9,]+ ", "", groups)
      pops <- vapply(groups, function(group) tryCatch(
        .gact_panukb_pop_code(list(sample_ancestry_category = group)),
        error = function(e) NA_character_), character(1))
      if (anyNA(pops)) return(NA_character_)
      if (length(unique(pops)) == 1L) pops[1L] else "Mixed"
    }, character(1))
  )
  keys <- c(paste0("label:", .gact_panukb_trait_key(records$label)),
            .gact_panukb_catalog_code(records$label))
  list(records = records, index = split(rep(seq_len(nrow(records)), 2L), keys))
}

.gact_panukb_index_phenotypes <- function(phenotype) {
  labels <- paste0("label:", .gact_panukb_trait_key(phenotype$description))
  codes <- .gact_panukb_provider_code(phenotype$trait_type %||% rep(NA_character_, nrow(phenotype)),
                                     phenotype$phenocode %||% rep(NA_character_, nrow(phenotype)))
  keys <- c(paste(labels, phenotype$pheno_sex),
            ifelse(is.na(codes), NA_character_, paste(codes, phenotype$pheno_sex)))
  split(rep(seq_len(nrow(phenotype)), 2L), keys)
}

.gact_panukb_candidates <- function(row, catalog, analysis) {
  indices <- unique(unlist(catalog$index[.gact_panukb_keys(row)], use.names = FALSE))
  records <- catalog$records[indices, , drop = FALSE]
  keep <- vapply(records$label, function(label)
    .gact_panukb_trait_score(row, label) > 0L, logical(1))
  keep <- keep & (is.na(records$ancestry) | records$ancestry == analysis)
  records <- records[keep, , drop = FALSE]
  records[!duplicated(records$accession), , drop = FALSE]
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
  handle <- curl::new_handle(range = sprintf("0-%d", as.integer(max_bytes) - 1L),
                             connecttimeout = 20, timeout = 45)
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

# GWAS coordinate conversion --------

.gact_read_liftover_chain <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- trimws(readLines(con))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  fields <- strsplit(lines, "[[:space:]]+")
  header <- startsWith(lines, "chain ") | startsWith(lines, "chain\t")
  if (!length(fields) || !header[1L] || any(lengths(fields[header]) != 13L))
    stop("Invalid UCSC chain header.")
  chain <- cumsum(header)[!header]
  blocks <- fields[!header]
  last <- !duplicated(chain, fromLast = TRUE)
  if (length(unique(chain)) != sum(header) ||
      !identical(lengths(blocks), ifelse(last, 1L, 3L)))
    stop("Invalid UCSC chain blocks.")
  headers <- do.call(rbind, fields[header])
  coordinates <- suppressWarnings(matrix(as.integer(headers[, c(4, 6, 7, 9, 11, 12)]),
                                          ncol = 6L))
  values <- suppressWarnings(matrix(as.integer(unlist(lapply(blocks, function(x)
    c(x, rep("0", 3L - length(x)))))), byrow = TRUE, ncol = 3L))
  if (anyNA(coordinates) || anyNA(values) || any(values < 0L) || any(values[, 1L] == 0L) ||
      any(headers[, 5L] != "+") || !all(headers[, 10L] %in% c("+", "-")) ||
      any(coordinates[, c(2, 5)] < 0L) || any(coordinates[, c(3, 6)] > coordinates[, c(1, 4)]))
    stop("Invalid UCSC chain coordinates.")
  blocks <- data.table::data.table(chain = chain, size = values[, 1L],
                                   dt = values[, 2L], dq = values[, 3L])
  blocks[, c("source_start", "query_start") := list(
    coordinates[chain[1L], 2L] + cumsum(size + dt) - size - dt,
    coordinates[chain[1L], 5L] + cumsum(size + dq) - size - dq), by = chain]
  if (any(blocks$source_start[last] + blocks$size[last] != coordinates[, 3L]) ||
      any(blocks$query_start[last] + blocks$size[last] != coordinates[, 6L]))
    stop("UCSC chain blocks disagree with their header spans.")
  blocks[, `:=`(
    chr = sub("^chr", "", headers[chain, 3L]), target_chr = sub("^chr", "", headers[chain, 8L]),
    lo = source_start + 1L, hi = source_start + size,
    reversed = headers[chain, 10L] == "-")]
  blocks[, offset := ifelse(reversed, coordinates[chain, 4L] - query_start + lo,
                            query_start - source_start)]
  blocks[, list(chr, target_chr, lo, hi, reversed, offset)]
}

.gact_liftover_map_chromosome <- function(variants, blocks, markers, chr) {
  bases <- c(A = "T", C = "G", G = "C", T = "A")
  snp <- variants$ref %in% names(bases) & variants$alt %in% names(bases) &
    variants$ref != variants$alt & !is.na(variants$pos) & variants$pos > 0
  snps <- variants[snp]
  duplicate <- duplicated(snps, by = c("pos", "ref", "alt")) |
    duplicated(snps, by = c("pos", "ref", "alt"), fromLast = TRUE)
  sites <- snps[!duplicate]
  hits <- unique(blocks[data.table::data.table(pos = unique(sites$pos)),
    on = list(lo <= pos, hi >= pos), allow.cartesian = TRUE, nomatch = 0L,
    list(pos = i.pos, target_chr, reversed,
         target_pos = as.integer(ifelse(reversed, offset - i.pos, offset + i.pos)))])
  mapped <- sites$pos %in% hits$pos
  ambiguous <- hits$pos[duplicated(hits$pos) | duplicated(hits$pos, fromLast = TRUE)]
  unique_hits <- hits[!pos %in% ambiguous]
  positions <- unique_hits[sites[mapped & !pos %in% ambiguous], on = "pos", nomatch = 0L]
  same_chr <- positions$target_chr == chr
  lifted <- positions[same_chr]
  lifted[, `:=`(nea = as.character(ifelse(reversed, bases[ref], ref)),
                 ea = as.character(ifelse(reversed, bases[alt], alt)))]
  if (anyDuplicated(markers, by = c("pos", "nea", "ea")))
    stop("Duplicate allele keys in the GRCh38 marker catalog on chromosome ", chr)
  at <- markers[lifted, on = list(pos = target_pos, nea, ea), which = TRUE]
  matched <- !is.na(at)
  result <- lifted[matched, list(chr, pos, ref, alt, target_pos, nea, ea, reversed)]
  result[, rsids := markers$rsids[at[matched]]]
  duplicate_target <- duplicated(result$rsids) | duplicated(result$rsids, fromLast = TRUE)
  map <- result[!duplicate_target]
  counts <- c(input = nrow(variants), unsupported = sum(!snp), duplicate_source = sum(duplicate),
              unmapped = sum(!mapped), ambiguous = sum(sites$pos %in% ambiguous),
              different_chromosome = sum(!same_chr),
              outside_backbone = sum(!lifted$target_pos %in% markers$pos),
              allele_mismatch = sum(!matched & lifted$target_pos %in% markers$pos),
              duplicate_target = sum(duplicate_target), retained = nrow(map))
  list(map = map, counts = as.list(counts))
}

#' Build a reusable GRCh37 to GRCh38 SNP map
#'
#' Lift a GRCh37 source variant list once with a UCSC hg19-to-hg38
#' chain. Process autosomal SNPs by chromosome and keep unique mappings
#' whose REF and ALT agree with the release marker catalog. Saves one
#' compressed table per chromosome, a reusable C++ index, QC counts, and checksums.
#' Tables use R gzip compression and work with data.table builds without gzip support.
#' @param GAlist GRCh38 database descriptor from the IGSR marker build.
#' @param variant_manifest Local GRCh37 variant list with chr (or chrom), pos, ref, and alt columns.
#'   Defaults to the Pan-UKB variant list cached in the release reference directory.
#' @param chain Local UCSC hg19-to-hg38 chain file, plain or gzip compressed.
#'   Defaults to a checksum-verified download cached in the release reference directory.
#' @param out_dir New output directory. Defaults to liftover-grch37-grch38 under the release marker directory.
#' @return Invisibly, the map directory. Pass it to [liftGWASStats()].
#' @export
buildLiftoverMap <- function(GAlist, variant_manifest = NULL, chain = NULL, out_dir = NULL) {
  if (!identical(.gact_catalog_build(GAlist), "GRCh38"))
    stop("A GRCh38 marker backbone is required.")
  if (is.null(out_dir)) out_dir <- .gact_liftover_dir(GAlist)
  parent <- normalizePath(dirname(out_dir), winslash = "/", mustWork = TRUE)
  out_dir <- file.path(parent, basename(out_dir))
  if (file.exists(out_dir)) {
    if (!is.null(variant_manifest) || !is.null(chain))
      stop("Liftover output already exists: ", out_dir)
    .gact_verify_liftover_map(GAlist, out_dir)
    message("Using verified liftover map at ", out_dir)
    return(invisible(out_dir))
  }
  defaults <- .gact_liftover_sources()
  sources <- list(if (is.null(variant_manifest)) defaults$variant_manifest,
                  if (is.null(chain)) defaults$chain, NULL)
  if (is.null(chain)) chain <- .gact_cache_liftover_source(GAlist, defaults$chain)
  if (is.null(variant_manifest))
    variant_manifest <- .gact_cache_liftover_source(GAlist, defaults$variant_manifest)
  catalog <- file.path(GAlist$dirs[["marker"]], "markers.txt.gz")
  provenance <- jsonlite::read_json(file.path(dirname(catalog), "marker_manifest.json"))
  if (!identical(provenance$build, "GRCh38") || !identical(provenance$effect_allele, "ALT"))
    stop("Expected a GRCh38 marker catalog with ALT as its effect allele.")
  inputs <- c(variant_manifest = variant_manifest, chain = chain, marker_catalog = catalog)
  hashes <- tools::md5sum(inputs)
  if (anyNA(hashes) || !identical(unname(hashes[3L]), provenance$checksums[["markers.txt.gz"]]))
    stop("Missing liftover input or marker catalog checksum mismatch.")
  connection <- gzfile(variant_manifest, "rt")
  header <- tryCatch(strsplit(readLines(connection, n = 1L), "\t", fixed = TRUE)[[1L]],
                     finally = close(connection))
  chr_column <- intersect(c("chr", "chrom"), header)
  required <- c("pos", "ref", "alt")
  if (length(chr_column) != 1L || anyDuplicated(header) || !all(required %in% header))
    stop("Variant list must have exactly one chromosome column (chr or chrom), plus pos, ref and alt")
  variants <- data.table::fread(variant_manifest, select = c(chr_column, required),
    colClasses = list(character = c(chr_column, "ref", "alt"), integer = "pos"))
  data.table::setnames(variants, chr_column, "chr")
  variants[, chr := sub("^chr", "", chr)]
  data.table::setkeyv(variants, "chr")
  markers <- data.table::fread(catalog, select = c("rsids", "chr", "pos", "nea", "ea"),
                                colClasses = list(character = "chr", integer = "pos"))
  data.table::setkeyv(markers, "chr")
  blocks <- .gact_read_liftover_chain(chain)
  data.table::setkeyv(blocks, "chr")
  stage <- tempfile(".gact-liftover-", tmpdir = parent)
  if (!dir.create(stage)) stop("Cannot stage liftover map.")
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  autosomes <- as.character(seq_len(22L))
  receipts <- lapply(autosomes, function(chromosome) {
    result <- .gact_liftover_map_chromosome(variants[chromosome, nomatch = 0L],
      blocks[chromosome, nomatch = 0L], markers[chromosome, nomatch = 0L], chromosome)
    chr <- chromosome
    file <- paste0("chr", chr, ".tsv.gz")
    path <- file.path(stage, file)
    plain <- sub("\\.gz$", "", path)
    data.table::fwrite(result$map, plain, sep = "\t", compress = "none")
    R.utils::gzip(plain, destname = path, remove = TRUE)
    message("Liftover chr", chr, ": ", result$counts$retained,
            " of ", result$counts$input, " variants retained")
    list(chromosome = chr, file = file, md5 = unname(tools::md5sum(path)), counts = result$counts)
  })
  populated <- Filter(function(x) x$counts$retained > 0L, receipts)
  if (!length(populated))
    stop("No variants matched the GRCh38 marker backbone.")
  index <- file.path(stage, "index.bin")
  shards <- vapply(populated, function(x) x$file, character(1L))
  .Call(C_gact_gwas_liftover_index, file.path(stage, shards), index)
  jsonlite::write_json(list(source_build = "GRCh37", build = "GRCh38",
    index = list(file = basename(index), md5 = unname(tools::md5sum(index))),
    inputs = Map(function(file, md5, source) {
      input <- list(file = basename(file), md5 = md5)
      if (!is.null(source)) input$url <- source$url
      input
    }, as.list(inputs), as.list(unname(hashes)), sources),
    non_autosomal = sum(!variants$chr %in% autosomes), chromosomes = receipts),
    file.path(stage, "manifest.json"), pretty = TRUE, auto_unbox = TRUE)
  if (file.exists(out_dir) || !file.rename(stage, out_dir)) stop("Cannot publish liftover map.")
  invisible(out_dir)
}

.gact_liftover_stat_keys <- c("chr", "pos", "ref", "alt")

.gact_liftover_sites <- function(stat) {
  if (!is.data.frame(stat) || !all(.gact_liftover_stat_keys %in% names(stat)))
    stop("Supply chr, pos, ref and alt columns.")
  pos <- stat$pos
  if (!is.numeric(pos)) pos <- suppressWarnings(as.numeric(as.character(pos)))
  data.table::data.table(chr = sub("^chr", "", as.character(stat$chr), ignore.case = TRUE),
    pos = pos,
    ref = toupper(as.character(stat$ref)), alt = toupper(as.character(stat$alt)))
}

.gact_liftover_map_manifest <- function(map_dir) {
  manifest <- jsonlite::read_json(file.path(map_dir, "manifest.json"))
  if (!identical(manifest$source_build, "GRCh37") || !identical(manifest$build, "GRCh38"))
    stop("Expected a GRCh37 to GRCh38 map.")
  manifest
}

.gact_read_liftover_shard <- function(map_dir, receipt) {
  path <- file.path(map_dir, receipt$file)
  if (!identical(unname(tools::md5sum(path)), receipt$md5))
    stop("Liftover map checksum mismatch: ", receipt$file)
  data.table::fread(path, colClasses = list(character = c("chr", "ref", "alt", "nea", "ea", "rsids"),
                                            integer = c("pos", "target_pos")))
}

#' Verify the genome build of a GWAS preview
#'
#' Compare distinct coordinate and allele-pair keys with both sides of the map.
#' A key supports a build when it matches that build only. Verification
#' requires at least `min_informative` such keys and 90 percent agreement.
#' An unresolved preview has an `NA` build and requires review.
#' @param preview Data frame with chr, pos, ref and alt columns.
#' @param map_dir Output directory from [buildLiftoverMap()].
#' @param min_informative Minimum number of distinct build-informative keys.
#' @return List with the verified build and counts for GRCh37-only, GRCh38-only,
#'   both, and neither.
#' @export
verifyGWASBuild <- function(preview, map_dir, min_informative = 10L) {
  if (!is.numeric(min_informative) || length(min_informative) != 1L ||
      !is.finite(min_informative) || min_informative < 1 || min_informative %% 1 != 0)
    stop("min_informative must be a positive integer.")
  sites <- .gact_liftover_sites(preview)
  sites[, c("ref", "alt") := list(pmin(ref, alt), pmax(ref, alt))]
  sites <- unique(sites)
  matches37 <- matches38 <- rep(FALSE, nrow(sites))
  for (receipt in .gact_liftover_map_manifest(map_dir)$chromosomes) {
    rows <- which(sites$chr == receipt$chromosome)
    if (!length(rows) || !receipt$counts$retained) next
    map <- .gact_read_liftover_shard(map_dir, receipt)
    if (!nrow(map)) next
    keys <- sites[rows]
    matches37[rows] <- !is.na(map[keys, on = .gact_liftover_stat_keys, which = TRUE])
    matches38[rows] <- !is.na(map[keys, on = list(chr, target_pos = pos, nea = ref, ea = alt),
                                   which = TRUE, mult = "first"])
    keys[, c("ref", "alt") := list(alt, ref)]
    matches37[rows] <- matches37[rows] |
      !is.na(map[keys, on = .gact_liftover_stat_keys, which = TRUE])
    matches38[rows] <- matches38[rows] |
      !is.na(map[keys, on = list(chr, target_pos = pos, nea = ref, ea = alt),
                 which = TRUE, mult = "first"])
  }
  support <- c(GRCh37 = sum(matches37 & !matches38), GRCh38 = sum(matches38 & !matches37))
  informative <- sum(support)
  build <- if (informative >= min_informative && max(support) / informative >= 0.9)
    names(support)[which.max(support)] else NA_character_
  list(build = build, grch37_only = unname(support[1L]), grch38_only = unname(support[2L]),
       both = sum(matches37 & matches38), neither = sum(!matches37 & !matches38))
}

#' Lift GRCh37 SNP summary statistics to GRCh38
#'
#' Join the original chr, pos, ref and alt to a saved liftover map.
#' Replaces coordinates and allele labels in the retained rows, preserving
#' order and association values, and adds the catalog rsID as marker.
#' Reports the retained and dropped row counts.
#' @param stat GRCh37 rows with chr, pos, ref and alt columns.
#' @param map_dir Output directory from [buildLiftoverMap()].
#' @return Data frame of lifted rows, with liftover counts in the `liftover_qc` attribute.
#' @export
liftGWASStats <- function(stat, map_dir) {
  keys <- .gact_liftover_stat_keys
  sites <- .gact_liftover_sites(stat)
  chromosome_rows <- split(seq_len(nrow(sites)), sites$chr)
  retained <- lapply(.gact_liftover_map_manifest(map_dir)$chromosomes, function(receipt) {
    rows <- chromosome_rows[[receipt$chromosome]]
    if (!length(rows) || !receipt$counts$retained) return(NULL)
    map <- .gact_read_liftover_shard(map_dir, receipt)
    if (!nrow(map)) return(NULL)
    at <- map[sites[rows], on = keys, which = TRUE]
    keep <- !is.na(at)
    list(rows = rows[keep], chr = map$chr[at[keep]], pos = map$target_pos[at[keep]],
         ref = map$nea[at[keep]], alt = map$ea[at[keep]], marker = map$rsids[at[keep]])
  })
  mapped <- data.table::rbindlist(retained)
  if (nrow(mapped)) {
    data.table::setorderv(mapped, "rows")
    result <- as.data.frame(stat)[mapped$rows, , drop = FALSE]
    fields <- c(keys, "marker")
    result[fields] <- as.data.frame(mapped[, fields, with = FALSE])
  } else {
    result <- as.data.frame(stat)[FALSE, , drop = FALSE]
    result$marker <- character()
  }
  attr(result, "liftover_qc") <- list(n_input = nrow(stat), n_lifted = nrow(result),
                                     n_dropped = nrow(stat) - nrow(result))
  result
}

# Study inspection and recipe generation --------

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

# Schemas, recipes and single-study ingestion --------

.gwas_required_fields <- c("chr", "pos", "ea", "nea", "b", "seb", "p")

.gwas_neglog10p_column <- function(name) {
 gsub("[^a-z0-9]", "", tolower(name)) %in%
  c("lp", "mlogp", "neglog10p", "minuslog10p", "neglog10pvalue", "minuslog10pvalue")
}

scoreTemplateMatch <- function(colnames, template_mapping) {

 clean <- function(x) tolower(gsub("[^a-z0-9]", "", x, ignore.case = TRUE))
 cols <- clean(colnames)

 score <- 0
 total <- 0

 for (f in names(template_mapping)) {

  aliases <- template_mapping[[f]]

  if (is.null(aliases) || all(is.na(aliases)) || length(aliases) == 0) next

  aliases <- clean(aliases)

  total <- total + 1

  if (any(aliases %in% cols)) {
   score <- score + 1
  }
 }

 if (total == 0) return(0)

 score / total
}

detectTemplate <- function(stat,
                           templates = c("UKB","FinnGen","SAIGE"),
                           threshold = 0.6) {

 coln <- colnames(stat)

 scores <- sapply(templates, function(tpl) {
  mapping <- getTemplateMapping(tpl)
  scoreTemplateMatch(coln, mapping)
 })

 if (all(scores == 0)) {
  return(list(template = NA, score = 0, all_scores = scores))
 }

 best <- names(which.max(scores))
 best_score <- max(scores)

 list(
  template = if (best_score >= threshold) best else NA,
  score = best_score,
  all_scores = scores
 )
}

resolveTemplateMapping <- function(stat, template) {

 tpl <- getTemplateMapping(template)
 cols <- colnames(stat)

 clean <- function(x) tolower(gsub("[^a-z0-9]", "", x, ignore.case = TRUE))
 clean_cols <- clean(cols)

 mapping <- setNames(vector("list", length(tpl)), names(tpl))

 for (f in names(tpl)) {

  aliases <- tpl[[f]]

  if (is.null(aliases) || all(is.na(aliases))) {
   mapping[[f]] <- NA_character_
   next
  }

  aliases_clean <- clean(aliases)

  match_idx <- which(clean_cols %in% aliases_clean)

  mapping[[f]] <- if (length(match_idx)) cols[match_idx[1]] else NA_character_
 }

 mapping
}

autoDetectSchema <- function(stat,
                             template_threshold = 0.6,
                             verbose = TRUE) {

 tpl <- detectTemplate(stat, threshold = template_threshold)

 # -------------------------
 # TEMPLATE MODE
 # -------------------------
 if (!is.na(tpl$template)) {

  if (verbose) {
   message("Template detected: ", tpl$template,
           " (score=", round(tpl$score,2), ")")
  }

  mapping <- resolveTemplateMapping(stat, tpl$template)

  return(list(
   mode = "template",
   template = tpl$template,
   mapping = mapping,
   confidence = setNames(rep(tpl$score, length(mapping)), names(mapping)),
   schema = NULL
  ))
 }

 # -------------------------
 # AUTO MODE
 # -------------------------
 if (verbose) message("Falling back to schema detection")

 schema <- detectStatSchema(stat)

 return(list(
  mode = "auto",
  template = NA,
  mapping = schema$mapping,
  confidence = schema$confidence,
  schema = schema
 ))
}


ingestWithAuto <- function(GAlist, stat, ..., confidence_threshold = 0.15) {

 det <- autoDetectSchema(stat)

 mapping <- det$mapping
 confidence <- det$confidence

 # -------------------------
 # Confidence check
 # -------------------------
 low_conf <- names(mapping)[
  !is.na(mapping) & confidence < confidence_threshold
 ]

 if (length(low_conf)) {
  stop("Low confidence fields: ",
       paste(low_conf, collapse = ", "),
       ". Please review mapping.")
 }

 # -------------------------
 # Ingest
 # -------------------------
 ingestStatDB(
  GAlist = GAlist,
  stat = stat,
  schema = mapping,
  ...
 )
}


cleanRecipe <- function(recipe) {

 `%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || is.na(a) || identical(a, "")) b else a
 }

 recipe$source    <- recipe$source    %||% "unknown"
 recipe$ancestry  <- recipe$ancestry  %||% "unknown"
 recipe$build     <- recipe$build     %||% "GRCh37"
 recipe$gender    <- recipe$gender    %||% "both"
 recipe$comments  <- recipe$comments  %||% ""

 # Ensure mapping is named list
 if (!is.null(recipe$mapping)) {
  recipe$mapping <- as.list(recipe$mapping)
 }

 recipe
}



validateRecipe <- function(recipe, level = c("minimal","standard","strict")) {

 level <- match.arg(level)

 errors <- character()
 warnings <- character()
 info <- character()

 is_empty <- function(x) {
  is.null(x) || length(x) == 0 || all(is.na(x)) || identical(x, "")
 }

 # -----------------------------
 # REQUIRED
 # -----------------------------
 if (is_empty(recipe$trait)) {
  errors <- c(errors, "Missing trait")
 }

 if (is_empty(recipe$type)) {
  errors <- c(errors, "Missing type (quantitative/binary)")
 }

 # -----------------------------
 # Mapping checks
 # -----------------------------
 if (is.null(recipe$mapping)) {

  errors <- c(errors, "Missing mapping")

 } else {

  required_fields <- .gwas_required_fields

  for (f in required_fields) {

   val <- recipe$mapping[[f]]

   if (is_empty(val)) {
    errors <- c(errors, paste("Mapping missing for:", f))
   } else if (!is.character(val)) {
    errors <- c(errors, paste("Mapping not character for:", f))
   }
  }

  # duplicates
  mapped_cols <- unlist(recipe$mapping)
  mapped_cols <- mapped_cols[!is.na(mapped_cols)]

  dup <- mapped_cols[duplicated(mapped_cols)]

  if (length(dup)) {
   warnings <- c(warnings,
                 paste("Duplicate column mapping:", paste(unique(dup), collapse=", ")))
  }

 }

 # -----------------------------
 # STANDARD
 # -----------------------------
 if (level %in% c("standard","strict")) {

  if (is_empty(recipe$source)) {
   warnings <- c(warnings, "Missing source → fallback used")
  }

  if (is_empty(recipe$ancestry)) {
   warnings <- c(warnings, "Missing ancestry → set to 'unknown'")
  }

  if (is_empty(recipe$build)) {
   warnings <- c(warnings, "Missing genome build")
  }

  if (is_empty(recipe$n) && is_empty(recipe$ncase)) {
   warnings <- c(warnings, "Missing sample size (n or ncase/ncontrol)")
  }
 }

 # -----------------------------
 # STRICT
 # -----------------------------
 if (level == "strict") {

  if (is_empty(recipe$category)) {
   warnings <- c(warnings, "Missing category")
  }

  if (is_empty(recipe$efo)) {
   warnings <- c(warnings, "Missing EFO term")
  }

  if (is_empty(recipe$reference)) {
   warnings <- c(warnings, "Missing reference (PMID)")
  }
 }

 # -----------------------------
 # Validation block
 # -----------------------------
 if (!is.null(recipe$validation)) {

  if (!isTRUE(recipe$validation$ok)) {
   errors <- c(errors, "Underlying data validation failed")
  }

  if (!is.null(recipe$validation$metrics$n_rows)) {
   if (recipe$validation$metrics$n_rows < 10) {
    warnings <- c(warnings, "Very small dataset preview (<10 rows)")
   }
  }
 }

 list(
  ok = length(errors) == 0,
  level = level,
  errors = unique(errors),
  warnings = unique(warnings),
  info = info
 )
}


# recipes <- list.files("recipes/", full.names = TRUE)
#
# for (f in recipes) {
#  r <- jsonlite::read_json(f)
#
#  v <- validateRecipe(r, "standard")
#
#  if (!v$ok) {
#   message("Skipping:", f)
#   next
#  }
#
#  runRecipe(r)
# }

runRecipe <- function(recipe, GAlist, db_path,
                      validate_level = "standard",
                      verbose = TRUE) {

 # -----------------------------
 # 1. Validate recipe
 # -----------------------------
 v <- validateRecipe(recipe, level = validate_level)

 if (!v$ok) {
  stop("Recipe failed validation:\n", paste(v$errors, collapse = "\n"))
 }

 if (length(v$warnings) > 0 && verbose) {
  message("Warnings:\n", paste(v$warnings, collapse = "\n"))
 }

 # -----------------------------
 # 2. Clean recipe
 # -----------------------------
 recipe <- cleanRecipe(recipe)

 # -----------------------------
 # 3. Load GWAS
 # -----------------------------
 has_file <- !is.null(recipe$file_name) &&
  nzchar(recipe$file_name) &&
  file.exists(recipe$file_name)

 has_url  <- !is.null(recipe$url) && nzchar(recipe$url)

 if (verbose) {
  message("Reading GWAS: ",
          if (has_file) recipe$file_name else recipe$url)
 }

 if (!has_file && has_url) {

  tmp <- tempfile(fileext = ".gz")

  tryCatch({
   download.file(recipe$url, tmp, mode = "wb")
  }, error = function(e) {
   stop("Failed to download GWAS from URL")
  })

  stat <- data.table::fread(tmp)

  # checksum
  recipe$file_md5 <- unname(tools::md5sum(tmp))

 } else if (has_file) {

  stat <- data.table::fread(recipe$file_name)
  recipe$file_md5 <- unname(tools::md5sum(recipe$file_name))

 } else {
  stop("No valid file or URL provided in recipe")
 }

 # -----------------------------
 # 4. Normalize
 # -----------------------------
 if (verbose) message("Normalizing...")

 stat_norm <- normalizeStatSchema(
  stat = stat,
  schema = recipe$mapping,
  effect_scale = if (is.null(recipe$effect_scale)) "auto" else recipe$effect_scale,
  p_scale = if (is.null(recipe$p_scale)) "auto" else recipe$p_scale
 )

 # -----------------------------
 # 5. Validate stats
 # -----------------------------
 val <- validateStatSchema(stat_norm)

 if (!val$ok) {
  stop("Stat validation failed:\n", paste(val$errors, collapse = "\n"))
 }

 if (verbose) {
  message("Rows:", nrow(stat_norm))
 }

 # -----------------------------
 # 6. Ingest
 # -----------------------------
 if (verbose) message("Updating database...")

 G_updated <- updateStatDB(
  GAlist = GAlist,
  stat = stat_norm,
  source = recipe$source,
  trait = recipe$trait,
  type = recipe$type,
  category = recipe$category,
  efo = recipe$efo,
  gender = recipe$gender,
  ancestry = recipe$ancestry,
  build = recipe$build,
  n = recipe$n,
  ncase = recipe$ncase,
  ncontrol = recipe$ncontrol,
  reference = recipe$reference,
  comments = recipe$comments,
  recipe = recipe,
  writeStatDB = TRUE,
  strand = recipe$strand
 )

 # -----------------------------
 # 7. Save DB
 # -----------------------------
 saveRDS(G_updated, db_path, compress = FALSE)

 if (verbose) message("Done:", recipe$trait)

 return(G_updated)
}

createRecipesFromManifest <- function(meta, recipes_dir) {

 dir.create(recipes_dir, showWarnings = FALSE, recursive = TRUE)

 for (i in seq_len(nrow(meta))) {

  row <- meta[i, ]

  recipe <- list(
   file_name = row$filename,        # or NA if remote only
   url = row$aws_link,             # ⭐ important

   trait = row$Description,
   type = ifelse(row$trait_type == "biomarkers", "quantitative", "binary"),
   category = row$Category,
   efo = row$trait_efo_terms,

   gender = row$pheno_sex,
   ancestry = "EUR",               # or from column if available
   build = "GRCh37",

   n = row$n_cases_full_cohort_both_sexes,
   ncase = row$n_cases_EUR,
   ncontrol = row$n_controls_EUR,

   source = "UKB",
   reference = row$trait_efo_terms,
   comments = paste("phenocode:", row$phenocode),

   # ⚠️ mapping can be default or template-based
   mapping = list(
    marker = "rsid",
    chr = "#CHROM",
    pos = "POS",
    ea = "ALT",
    nea = "REF",
    eaf = "AF",
    b = "BETA",
    seb = "SE",
    p = "P",
    n = NA,
    ncase = NA,
    ncontrol = NA,
    info = NA
   ),

   created_at = as.character(Sys.time())
  )

  file_name <- file.path(
   recipes_dir,
   paste0(row$phenocode, ".json")
  )

  jsonlite::write_json(recipe, file_name, pretty = TRUE, auto_unbox = TRUE)
 }
}

validateAllRecipes <- function(recipes_dir) {

 files <- list.files(recipes_dir, "*.json", full.names = TRUE)

 for (f in files) {

  r <- jsonlite::read_json(f)
  v <- validateRecipe(r, "standard")

  if (!v$ok) {
   message("❌ Invalid:", f)
   print(v$errors)
  } else {
   message("✅ OK:", f)
  }
 }
}

runRecipes <- function(recipes_dir, GAlist, db_path) {

 files <- list.files(recipes_dir, "*.json", full.names = TRUE)

 for (f in files) {

  message("Processing:", f)

  recipe <- jsonlite::read_json(f)

  tryCatch({
   GAlist <- runRecipe(recipe, GAlist, db_path)
  }, error = function(e) {
   message("FAILED:", f)
  })
 }

 return(GAlist)
}

getTemplateMapping <- function(template = c("UKB","FinnGen","SAIGE")) {

 template <- match.arg(template)

 # canonical fields
 canon <- c(
  "marker","chr","pos","ea","nea","eaf",
  "b","seb","p","n","ncase","ncontrol","info"
 )

 mapping <- switch(template,

                   # -------------------------
                   # UK Biobank (Pan-UKB / BOLT / Neale)
                   # -------------------------
                   UKB = list(
                    marker   = c("rsid","variant","markername"),
                    chr      = c("#chrom","chr"),
                    pos      = c("pos","position"),
                    ea       = c("alt","allele1","effect_allele"),
                    nea      = c("ref","allele0","other_allele"),
                    eaf      = c("af","alt_af","eaf"),
                    b        = c("beta","b","effect"),
                    seb      = c("se","standard_error"),
                    p        = c("pval","p","p_bolt_lmm","p.value"),
                    n        = c("n","n_complete_samples"),
                    ncase    = c("ncase"),
                    ncontrol = c("ncontrol"),
                    info     = c("info")
                   ),

                   # -------------------------
                   # FinnGen
                   # -------------------------
                   FinnGen = list(
                    marker   = c("rsid"),
                    chr      = c("#chrom","chr"),
                    pos      = c("pos"),
                    ea       = c("alt"),
                    nea      = c("ref"),
                    eaf      = c("af_alt","eaf"),
                    b        = c("beta"),
                    seb      = c("sebeta","se"),
                    p        = c("pval","p"),
                    n        = c("n_total"),
                    ncase    = c("n_cases"),
                    ncontrol = c("n_controls"),
                    info     = c("info")
                   ),

                   # -------------------------
                   # SAIGE
                   # -------------------------
                   SAIGE = list(
                    marker   = c("markerid","snpid"),
                    chr      = c("chr"),
                    pos      = c("pos"),
                    ea       = c("allele2","alt"),
                    nea      = c("allele1","ref"),
                    eaf      = c("af","eaf"),
                    b        = c("beta"),
                    seb      = c("se"),
                    p        = c("p.value","pval"),
                    n        = c("n"),
                    ncase    = c("ncase"),
                    ncontrol = c("ncontrol"),
                    info     = c("info")
                   )
 )

 # ensure all canonical fields exist
 missing <- setdiff(canon, names(mapping))
 if (length(missing) > 0) {
  for (m in missing) mapping[[m]] <- NA
 }

 return(mapping)
}


#' Read GWAS alias dictionary from YAML
#'
#' @param path Optional path to YAML file. If `NULL`, uses
#'   `inst/extdata/gwas_aliases.yml` from the package.
#' @param must_exist Logical; if `TRUE`, throws an error when file is missing.
#' @return A normalized alias dictionary list.
#' @keywords internal
read_aliases_yaml <- function(path = NULL, must_exist = TRUE) {
  if (is.null(path)) {
    path <- system.file("extdata", "gwas_aliases.yml", package = "gact")
  }

  if (!nzchar(path) || !file.exists(path)) {
    if (must_exist) stop("Alias YAML not found: ", path)
    return(NULL)
  }

  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read alias dictionaries.")
  }

  x <- yaml::read_yaml(path)
  if (is.null(x$fields) || !is.list(x$fields)) {
    stop("Alias YAML must contain a top-level 'fields:' list.")
  }

  clean_name <- function(z) tolower(gsub("[^a-z0-9]", "", as.character(z), ignore.case = TRUE))
  `%||%` <- function(a, b) if (is.null(a)) b else a

  out <- list(version = x$version %||% NA, fields = list())
  for (fld in names(x$fields)) {
    node <- x$fields[[fld]]
    aliases <- unique(clean_name(node$aliases %||% character(0)))
    aliases <- aliases[nzchar(aliases)]
    regex <- unique(as.character(node$regex %||% character(0)))
    out$fields[[fld]] <- list(aliases = aliases, regex = regex)
  }
  out
}

#' Apply alias dictionary to column names and return alias score matrix
#'
#' @param col_names Character vector of column names from incoming data.
#' @param alias_dict Alias dictionary from [read_aliases_yaml()].
#' @param fields Canonical fields to score.
#' @return A list with `score` matrix (`field x column`) and `hits`.
#' @keywords internal
apply_alias_dictionary <- function(col_names, alias_dict, fields) {
  clean_name <- function(z) tolower(gsub("[^a-z0-9]", "", as.character(z), ignore.case = TRUE))
  clean_cols <- clean_name(col_names)

  score <- matrix(0, nrow = length(fields), ncol = length(col_names),
                  dimnames = list(fields, col_names))
  hits <- setNames(vector("list", length(fields)), fields)

  if (is.null(alias_dict) || is.null(alias_dict$fields)) {
    return(list(score = score, hits = hits))
  }

  for (f in fields) {
    fnode <- alias_dict$fields[[f]]
    if (is.null(fnode)) next

    alias_hits <- rep(FALSE, length(clean_cols))
    regex_hits <- rep(FALSE, length(clean_cols))

    if (length(fnode$aliases)) alias_hits <- clean_cols %in% fnode$aliases

    if (length(fnode$regex)) {
      regex_hits <- vapply(col_names, function(nm) {
        any(vapply(fnode$regex, function(rx) {
          #grepl(rx, nm, ignore.case = TRUE, perl = TRUE)
         clean_nm <- tolower(gsub("[^a-z0-9]", "", nm, ignore.case = TRUE))
         grepl(rx, clean_nm, perl = TRUE)
        }, logical(1)))
      }, logical(1))
    }

    sc <- numeric(length(col_names))
    sc[alias_hits] <- pmax(sc[alias_hits], 1.0)
    sc[regex_hits] <- pmax(sc[regex_hits], 0.7)

    score[f, ] <- sc
    hits[[f]] <- col_names[alias_hits | regex_hits]
  }

  list(score = score, hits = hits)
}

#' Detect GWAS summary-statistics schema using column names and values
#'
#' Scores candidate columns for canonical GWAS fields and returns a mapping
#' with confidence scores, alternatives, and warnings.
#'
#' Canonical fields are: `marker`, `chr`, `pos`, `ea`, `nea`, `eaf`, `b`,
#' `seb`, `p`, `n`, `ncase`, `ncontrol`, and `info`.
#'
#' @param stat A `data.frame` containing summary statistics.
#' @param sample_n Optional integer; number of rows to sample for detection.
#' @param alias_path Optional path to alias YAML. Defaults to package extdata file.
#' @param alias_dict Optional pre-loaded alias dictionary from [read_aliases_yaml()].
#' @param w_name Weight for name-pattern score.
#' @param w_value Weight for value-distribution score.
#' @param w_alias Weight for alias-dictionary score.
#' @return A list with elements `mapping`, `confidence`, `candidates`,
#'   `warnings`, `details`, and `alias_hits`.
#' @export

detectStatSchema <- function(stat, sample_n = 100000,
                             alias_path = NULL, alias_dict = NULL,
                             w_name = 0.45, w_value = 0.30, w_alias = 0.25) {
 if (is.null(stat)) stop("`stat` is NULL")
 if (!is.data.frame(stat)) stat <- as.data.frame(stat)
 if (!nrow(stat)) stop("`stat` has zero rows")

 if (!is.null(sample_n) && is.finite(sample_n) && nrow(stat) > sample_n) {
  idx <- sort(sample.int(nrow(stat), sample_n))
  stat <- stat[idx, , drop = FALSE]
 }

 col_names <- colnames(stat)
 clean_names <- tolower(gsub("[^a-z0-9]", "", col_names, ignore.case = TRUE))

 # ---- FIX 1: load alias_dict BEFORE use ----
 if (is.null(alias_dict)) {
  alias_dict <- tryCatch(
   read_aliases_yaml(path = alias_path, must_exist = FALSE),
   error = function(e) NULL
  )
 }

 fields <- c("marker", "chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n", "ncase", "ncontrol", "info")

 # ---- FIX 2: safe name_patterns construction ----
 name_patterns <- setNames(vector("list", length(fields)), fields)

 if (!is.null(alias_dict) && !is.null(alias_dict$fields)) {
  for (f in intersect(names(alias_dict$fields), fields)) {
   fnode <- alias_dict$fields[[f]]
   if (length(fnode$aliases)) {
    name_patterns[[f]] <- paste0("^", fnode$aliases, "$")
   }
  }
 }

 # score_name <- function(field) {
 #  pats <- name_patterns[[field]]
 #  if (length(pats) == 0) return(rep(0, length(clean_names)))
 #  sapply(clean_names, function(nm) {
 #   if (any(grepl(paste(pats, collapse = "|"), nm, perl = TRUE))) 1 else 0
 #  })
 # }
 score_name <- function(field) {
  pats <- name_patterns[[field]]
  fnode <- alias_dict$fields[[field]]

  sapply(seq_along(col_names), function(i) {
   nm <- col_names[i]
   clean_nm <- clean_names[i]
   if (clean_nm == field || (field == "p" && .gwas_neglog10p_column(nm))) return(1)

   alias_hit <- FALSE
   if (length(pats)) {
    alias_hit <- any(grepl(paste(pats, collapse = "|"), clean_nm, perl = TRUE))
   }

   regex_hit <- FALSE
   if (!is.null(fnode$regex) && length(fnode$regex)) {
    regex_hit <- any(vapply(fnode$regex, function(rx) {
     grepl(rx, nm, ignore.case = TRUE, perl = TRUE) ||
      grepl(rx, clean_nm, perl = TRUE)
    }, logical(1)))
   }

   if (alias_hit) return(1.0)
   if (regex_hit) return(0.8)

   return(0)
  })
 }

 get_col <- function(j) stat[[j]]
 is_numeric_col <- function(x) is.numeric(x) || is.integer(x)

 value_scores <- lapply(seq_along(col_names), function(j) {
  x <- get_col(j)
  x_char <- toupper(trimws(as.character(x)))
  x_num <- suppressWarnings(as.numeric(x))

  non_na_num <- x_num[!is.na(x_num)]
  non_na_char <- x_char[!is.na(x_char) & nzchar(x_char)]

  frac_numeric <- if (!length(x)) 0 else mean(!is.na(x_num))
  uniq <- length(unique(non_na_char))

  chr_score <- if (length(non_na_num)) mean(non_na_num %in% c(1:26, 23, 24, 25, 26)) else 0
  pos_score <- if (length(non_na_num)) mean(non_na_num > 1000 & non_na_num < 4e8) else 0

  allele_score <- if (length(non_na_char)) mean(grepl("^[ACGT]+$", non_na_char)) else 0

  p_score <- 0
  if (length(non_na_num)) {
   if (.gwas_neglog10p_column(col_names[j])) {
    p_score <- mean(is.finite(non_na_num) & non_na_num >= 0)
   } else {
    p_score <- mean(non_na_num >= 0 & non_na_num <= 1)
   }
  }

  eaf_score <- 0
  if (length(non_na_num)) {
   in01 <- mean(non_na_num >= 0 & non_na_num <= 1)
   med <- stats::median(non_na_num)
   shape <- as.numeric(med > 0.01 && med < 0.99)
   eaf_score <- in01 + 0.2 * shape
  }

  seb_score <- if (length(non_na_num)) mean(non_na_num > 0 & non_na_num < 1) else 0

  b_score <- if (length(non_na_num)) mean(is.finite(non_na_num)) else 0

  n_score <- if (length(non_na_num)) mean(non_na_num >= 100) else 0

  info_score <- 0
  if (length(non_na_num)) {
   info_score <- mean(non_na_num >= 0 & non_na_num <= 1) +
    0.3 * as.numeric(stats::median(non_na_num) > 0.7)
  }

  # ---- FIX 5: marker_score restored ----
  marker_score <- as.numeric(!is_numeric_col(x)) * as.numeric(uniq > 1000)
  rs_frac <- if (length(non_na_char)) mean(grepl("^rs[0-9]+$", non_na_char)) else 0
  marker_score <- max(marker_score, rs_frac)

  list(
   chr = chr_score,
   pos = pos_score,
   ea = allele_score,
   nea = allele_score,
   p = p_score,
   eaf = eaf_score,
   seb = seb_score,
   b = b_score,
   n = n_score,
   ncase = n_score,
   ncontrol = n_score,
   info = info_score,
   marker = marker_score,
   frac_numeric = frac_numeric,
   uniq = uniq
  )
 })

 score_tbl <- matrix(0, nrow = length(fields), ncol = length(col_names),
                     dimnames = list(fields, col_names))

 alias_res <- apply_alias_dictionary(col_names = col_names,
                                     alias_dict = alias_dict,
                                     fields = fields)
 alias_score_tbl <- alias_res$score

 for (f in fields) {
  nscore <- score_name(f)
  vscore <- vapply(value_scores, function(v) v[[f]], numeric(1))
  score_tbl[f, ] <- w_name * nscore +
   w_value * pmin(vscore, 1.25) +
   w_alias * alias_score_tbl[f, ]
  if (f %in% c(.gwas_required_fields, "eaf")) {
   score_tbl[f, nscore == 0 & alias_score_tbl[f, ] == 0] <- 0
  }
 }

 warnings <- character(0)

 pick_best <- function(field, exclude = character()) {
  scores <- score_tbl[field, ]
  if (length(exclude)) scores[names(scores) %in% exclude] <- -Inf

  ord <- order(scores, decreasing = TRUE)
  top <- scores[ord[1]]
  second <- if (length(ord) > 1) scores[ord[2]] else NA_real_

  if (!is.finite(top)) {
   return(list(best = NA_character_, confidence = 0,
               candidates = character(0), ambiguous = FALSE))
  }

  # ---- FIX 6: tie-breaking ----
  if (!is.na(second) && top == second) {
   alias_bonus <- alias_score_tbl[field, ]
   ord <- order(scores, alias_bonus, decreasing = TRUE)
  }

  best <- names(scores)[ord[1]]
  confidence <- if (is.na(second) || is.infinite(second)) top else top - second
  ambiguous <- !is.na(second) && abs(top - second) < 0.05

  list(best = if (is.finite(top) && top > 0.25) best else NA_character_,
       confidence = max(0, confidence),
       candidates = names(scores)[ord[seq_len(min(3, length(ord)))]],
       ambiguous = ambiguous)
 }

 mapping <- setNames(rep(NA_character_, length(fields)), fields)
 confidence <- setNames(rep(0, length(fields)), fields)
 candidates <- vector("list", length(fields)); names(candidates) <- fields

 used <- character(0)
 for (f in c("marker", "chr", "pos", "p", "b", "seb", "eaf", "n", "ncase", "ncontrol", "info")) {
  p <- pick_best(f, exclude = used)
  mapping[f] <- p$best
  confidence[f] <- p$confidence
  candidates[[f]] <- p$candidates
  if (p$ambiguous) warnings <- c(warnings, paste("Ambiguous mapping for", f))
  if (!is.na(p$best)) used <- c(used, p$best)
 }

 ea_pick <- pick_best("ea")

 nea_scores <- score_tbl["nea", ]
 if (!is.na(ea_pick$best)) nea_scores[ea_pick$best] <- -Inf

 nea_ord <- order(nea_scores, decreasing = TRUE)
 nea_best <- names(nea_scores)[nea_ord[1]]

 mapping["ea"] <- if (!is.na(ea_pick$best) && score_tbl["ea", ea_pick$best] > 0.25) ea_pick$best else NA_character_
 mapping["nea"] <- if (is.finite(nea_scores[nea_best]) && nea_scores[nea_best] > 0.25) nea_best else NA_character_

 confidence["ea"] <- ea_pick$confidence
 confidence["nea"] <- if (length(nea_ord) > 1) nea_scores[nea_ord[1]] - nea_scores[nea_ord[2]] else nea_scores[nea_ord[1]]

 candidates[["ea"]] <- ea_pick$candidates
 candidates[["nea"]] <- names(nea_scores)[nea_ord[seq_len(min(3, length(nea_ord)))]]

 required <- .gwas_required_fields
 missing_required <- required[is.na(mapping[required])]
 if (length(missing_required)) {
  warnings <- c(warnings, paste("Missing required fields:", paste(missing_required, collapse = ", ")))
 }

 low_conf <- names(confidence)[confidence < 0.1 & !is.na(mapping)]
 if (length(low_conf)) {
  warnings <- c(warnings, paste("Low-confidence mappings:", paste(low_conf, collapse = ", ")))
 }

 list(
  mapping = mapping,
  confidence = confidence,
  candidates = candidates,
  alias_hits = alias_res$hits,
  warnings = unique(warnings),
  details = as.data.frame(t(score_tbl), stringsAsFactors = FALSE)
 )
}

# detectStatSchema <- function(stat, sample_n = 100000,
#                              alias_path = NULL, alias_dict = NULL,
#                              w_name = 0.45, w_value = 0.30, w_alias = 0.25) {
#   if (is.null(stat)) stop("`stat` is NULL")
#   if (!is.data.frame(stat)) stat <- as.data.frame(stat)
#   if (!nrow(stat)) stop("`stat` has zero rows")
#
#   if (!is.null(sample_n) && is.finite(sample_n) && nrow(stat) > sample_n) {
#     idx <- sort(sample.int(nrow(stat), sample_n))
#     stat <- stat[idx, , drop = FALSE]
#   }
#
#   col_names <- colnames(stat)
#   clean_names <- tolower(gsub("[^a-z0-9]", "", col_names))
#
#   name_patterns <- lapply(alias_dict$fields, function(fnode) {
#    if (length(fnode$aliases)) {
#     paste0("^", fnode$aliases, "$")
#    } else character(0)
#   })
#
#   score_name <- function(field) {
#     pats <- name_patterns[[field]]
#     sapply(clean_names, function(nm) {
#       if (any(grepl(paste(pats, collapse = "|"), nm, perl = TRUE))) 1 else 0
#     })
#   }
#
#   get_col <- function(j) stat[[j]]
#
#   is_numeric_col <- function(x) is.numeric(x) || is.integer(x)
#
#   value_scores <- lapply(seq_along(col_names), function(j) {
#     x <- get_col(j)
#     x_char <- toupper(trimws(as.character(x)))
#     x_num <- suppressWarnings(as.numeric(x))
#
#     non_na_num <- x_num[!is.na(x_num)]
#     non_na_char <- x_char[!is.na(x_char) & nzchar(x_char)]
#
#     frac_numeric <- if (!length(x)) 0 else mean(!is.na(x_num))
#     uniq <- length(unique(non_na_char))
#
#     chr_score <- 0
#     if (length(non_na_num)) {
#       chr_score <- mean(non_na_num %in% c(1:26, 23, 24, 25, 26))
#     }
#
#     pos_score <- 0
#     if (length(non_na_num)) {
#       pos_score <- mean(non_na_num > 1000 & non_na_num < 4e8)
#     }
#
#     allele_score <- 0
#     if (length(non_na_char)) {
#       #allele_score <- mean(non_na_char %in% c("A", "C", "G", "T"))
#       allele_score <- mean(grepl("^[ACGT]+$", non_na_char))
#     }
#
#     p_score <- 0
#     #if (length(non_na_num)) {
#     #  p_score <- mean(non_na_num >= 0 & non_na_num <= 1)
#     #  tiny <- mean(non_na_num < 1e-6)
#     #  p_score <- p_score + 0.2 * tiny
#     #}
#     if (length(non_na_num)) {
#      if (max(non_na_num, na.rm = TRUE) > 50) {
#       # likely -log10(p)
#       p_score <- 0.3
#      } else {
#       p_score <- mean(non_na_num >= 0 & non_na_num <= 1)
#      }
#     }
#
#     eaf_score <- 0
#     if (length(non_na_num)) {
#       in01 <- mean(non_na_num >= 0 & non_na_num <= 1)
#       med <- stats::median(non_na_num)
#       shape <- as.numeric(med > 0.01 && med < 0.99)
#       eaf_score <- in01 + 0.2 * shape
#     }
#
#     seb_score <- 0
#     if (length(non_na_num)) {
#       seb_score <- mean(non_na_num > 0 & non_na_num < 1)
#     }
#
#     b_score <- 0
#     #if (length(non_na_num)) {
#     #  b_score <- as.numeric(any(non_na_num < 0)) + 0.5 * mean(abs(non_na_num) < 5)
#     #}
#     if (length(non_na_num)) {
#      if (all(non_na_num > 0, na.rm = TRUE) && median(non_na_num, na.rm = TRUE) > 0.5) {
#       # likely OR
#       b_score <- b_score + 0.3
#      } else if (any(non_na_num < 0, na.rm = TRUE)) {
#       # likely beta
#       b_score <- b_score + 0.5
#      }
#     }
#
#     n_score <- 0
#     if (length(non_na_num)) {
#       n_score <- mean(non_na_num >= 100)
#     }
#
#     info_score <- 0
#     if (length(non_na_num)) {
#       info_score <- mean(non_na_num >= 0 & non_na_num <= 1) +
#         0.3 * as.numeric(stats::median(non_na_num) > 0.7)
#     }
#
#     #marker_score <- as.numeric(!is_numeric_col(x)) * as.numeric(uniq > 1000)
#     rs_frac <- mean(grepl("^rs[0-9]+$", non_na_char))
#     marker_score <- max(marker_score, rs_frac)
#
#     list(
#       chr = chr_score,
#       pos = pos_score,
#       ea = allele_score,
#       nea = allele_score,
#       p = p_score,
#       eaf = eaf_score,
#       seb = seb_score,
#       b = b_score,
#       n = n_score,
#       ncase = n_score,
#       ncontrol = n_score,
#       info = info_score,
#       marker = marker_score,
#       frac_numeric = frac_numeric,
#       uniq = uniq
#     )
#   })
#
#   fields <- c("marker", "chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n", "ncase", "ncontrol", "info")
#   score_tbl <- matrix(0, nrow = length(fields), ncol = length(col_names),
#                       dimnames = list(fields, col_names))
#
#   if (is.null(alias_dict)) {
#     alias_dict <- tryCatch(
#       read_aliases_yaml(path = alias_path, must_exist = FALSE),
#       error = function(e) NULL
#     )
#   }
#   alias_res <- apply_alias_dictionary(col_names = col_names,
#                                       alias_dict = alias_dict,
#                                       fields = fields)
#   alias_score_tbl <- alias_res$score
#
#   for (f in fields) {
#     nscore <- score_name(f)
#     vscore <- vapply(value_scores, function(v) v[[f]], numeric(1))
#     score_tbl[f, ] <- w_name * nscore +
#       w_value * pmin(vscore, 1.25) +
#       w_alias * alias_score_tbl[f, ]
#   }
#
#   pick_best <- function(field, exclude = character()) {
#     scores <- score_tbl[field, ]
#     if (length(exclude)) scores[names(scores) %in% exclude] <- -Inf
#     ord <- order(scores, decreasing = TRUE)
#     best <- names(scores)[ord[1]]
#     top <- scores[ord[1]]
#     second <- if (length(ord) > 1) scores[ord[2]] else NA_real_
#     if (top == second) {
#      # prefer alias-supported column
#      alias_bonus <- alias_score_tbl[field, ord]
#      ord <- ord[order(alias_bonus[ord], decreasing = TRUE)]
#     }
#     confidence <- if (is.na(second) || is.infinite(second)) top else top - second
#     if (!is.na(second) && abs(top - second) < 0.05) {
#      warnings <- c(warnings, paste("Ambiguous mapping for", field))
#     }
#     list(best = if (is.finite(top) && top > 0.25) best else NA_character_,
#          confidence = max(0, confidence),
#          candidates = names(scores)[ord[seq_len(min(3, length(ord)))]] )
#   }
#
#   mapping <- setNames(rep(NA_character_, length(fields)), fields)
#   confidence <- setNames(rep(0, length(fields)), fields)
#   candidates <- vector("list", length(fields)); names(candidates) <- fields
#
#   used <- character(0)
#   for (f in c("marker", "chr", "pos", "p", "b", "seb", "eaf", "n", "ncase", "ncontrol", "info")) {
#     p <- pick_best(f, exclude = used)
#     mapping[f] <- p$best
#     confidence[f] <- p$confidence
#     candidates[[f]] <- p$candidates
#     if (!is.na(p$best)) used <- c(used, p$best)
#   }
#
#   ea_pick <- pick_best("ea")
#   nea_scores <- score_tbl["nea", ]
#   #if (!is.na(ea_pick$best)) nea_scores[ea_pick$best] <- -Inf
#   if (!is.na(ea_pick$best) && ea_pick$best == nea_best) {
#    nea_scores[nea_best] <- -Inf
#    nea_best <- names(sort(nea_scores, decreasing = TRUE))[1]
#   }
#   nea_ord <- order(nea_scores, decreasing = TRUE)
#   nea_best <- names(nea_scores)[nea_ord[1]]
#
#   mapping["ea"] <- if (!is.na(ea_pick$best) && score_tbl["ea", ea_pick$best] > 0.25) ea_pick$best else NA_character_
#   mapping["nea"] <- if (is.finite(nea_scores[nea_best]) && nea_scores[nea_best] > 0.25) nea_best else NA_character_
#   confidence["ea"] <- ea_pick$confidence
#   confidence["nea"] <- if (length(nea_ord) > 1) nea_scores[nea_ord[1]] - nea_scores[nea_ord[2]] else nea_scores[nea_ord[1]]
#   candidates[["ea"]] <- ea_pick$candidates
#   candidates[["nea"]] <- names(nea_scores)[nea_ord[seq_len(min(3, length(nea_ord)))]]
#
#   warnings <- character(0)
#   required <- c("chr", "pos", "ea", "nea", "p")
#   missing_required <- required[is.na(mapping[required])]
#   if (length(missing_required)) warnings <- c(warnings, paste("Missing required fields:", paste(missing_required, collapse = ", ")))
#   low_conf <- names(confidence)[confidence < 0.1 & !is.na(mapping)]
#   if (length(low_conf)) warnings <- c(warnings, paste("Low-confidence mappings:", paste(low_conf, collapse = ", ")))
#
#   list(
#     mapping = mapping,
#     confidence = confidence,
#     candidates = candidates,
#     alias_hits = alias_res$hits,
#     warnings = unique(warnings),
#     details = as.data.frame(t(score_tbl), stringsAsFactors = FALSE)
#   )
# }

#' Normalize GWAS summary statistics to gact canonical column names
#'
#' @param stat A `data.frame` with raw summary statistics.
#' @param schema Output from [detectStatSchema()] or a named character vector.
#' @param drop_unmapped Logical; if `TRUE`, return only canonical + passthrough columns.
#' @param require_fields Canonical fields required to be present after normalization.
#'   Defaults to `chr`, `pos`, `ea`, `nea`, `b`, `seb`, and `p`.
#' @param effect_scale Scale of the mapped effect column. `"auto"` converts
#'   columns named `OR` or `odds_ratio` (ignoring case and punctuation).
#'   Use `"beta"` for beta or log-odds values, or `"or"` to convert an odds
#'   ratio from any column name. Supply `seb` on the beta or log-odds scale.
#' @param p_scale Scale of the mapped p-value column. `"auto"` converts
#'   `LP`, `mlogp`, `neglog10p`, and `minus_log10_p` names (ignoring case and
#'   punctuation). Use `"p"` for ordinary p-values or `"neglog10p"` for
#'   negative base-10 log p-values from any column name.
#' @return A normalized `data.frame` with canonical columns where possible.
#' @export
normalizeStatSchema <- function(stat,
                                schema,
                                drop_unmapped = FALSE,
                                require_fields = .gwas_required_fields,
                                effect_scale = c("auto", "beta", "or"),
                                p_scale = c("auto", "p", "neglog10p")) {
 if (is.null(stat)) stop("`stat` is NULL")
 if (!is.data.frame(stat)) stat <- as.data.frame(stat)

 mapping <- if (is.list(schema) && !is.null(schema$mapping)) schema$mapping else schema
 if (is.null(names(mapping))) stop("`schema` mapping must be a named vector")
 effect_scale <- match.arg(effect_scale)
 p_scale <- match.arg(p_scale)

 out <- stat
 canon <- c("marker", "chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n", "ncase", "ncontrol", "info")

 for (nm in canon) {
  src <- if (nm %in% names(mapping)) mapping[[nm]] else NULL
  if (!is.null(src) && !is.na(src) && src %in% colnames(stat)) {
   out[[nm]] <- stat[[src]]
  }
 }

 # ---- Standardize basic fields ----
 if (!is.null(out$chr)) {
  chr <- as.character(out$chr)
  chr <- gsub("^chr", "", chr, ignore.case = TRUE)
  chr[chr %in% c("X", "x")] <- "23"
  chr[chr %in% c("Y", "y")] <- "24"
  chr[chr %in% c("MT", "M", "mt", "m")] <- "25"
  out$chr <- suppressWarnings(as.integer(chr))
 }

 if (!is.null(out$pos)) out$pos <- suppressWarnings(as.integer(out$pos))
 if (!is.null(out$ea)) out$ea <- toupper(as.character(out$ea))
 if (!is.null(out$nea)) out$nea <- toupper(as.character(out$nea))

 for (nm in c("eaf", "b", "seb", "p", "n", "ncase", "ncontrol", "info")) {
  values <- out[[nm]]
  if (is.null(values)) next
  if (is.factor(values)) values <- as.character(values)
  numeric_values <- suppressWarnings(as.numeric(values))
  if (nm %in% c("b", "seb", "eaf", "p", "n") && is.character(values)) {
   invalid <- is.na(numeric_values) & !is.nan(numeric_values) &
    !is.na(values) & !(trimws(values) %in% c("", "NA"))
   if (any(invalid)) stop("Non-numeric values in ", nm, ": ", sum(invalid), ".")
  }
  out[[nm]] <- numeric_values
 }

 # ---- AUTO TRANSFORMATIONS ----
 transformations <- character(0)

 if (p_scale == "auto") {
  p_col <- if ("p" %in% names(mapping)) mapping[["p"]] else "p"
  p_scale <- if (isTRUE(.gwas_neglog10p_column(p_col))) "neglog10p" else "p"
 }
 if (p_scale == "neglog10p" && "p" %in% names(out)) {
  p_values <- out[["p"]]
  observed <- p_values[!is.na(p_values) | is.nan(p_values)]
  if (any(!is.finite(observed) | observed < 0)) {
   stop("-log10(p) values must be finite and nonnegative.")
  }
  out$p <- 10^(-p_values)
  transformations <- c(transformations, "-log10(p) -> p")
 }

 effect_col <- if ("b" %in% names(mapping)) mapping[["b"]] else NA_character_
 if (effect_scale == "auto") {
  effect_name <- gsub("[^a-z0-9]", "", tolower(effect_col))
  effect_scale <- if (effect_name %in% c("or", "oddsratio")) "or" else "beta"
 }
 if (effect_scale == "or" && "b" %in% names(out)) {
  if (any(!is.na(out$b) & (!is.finite(out$b) | out$b <= 0))) {
   stop("Odds ratios must be finite and strictly positive.")
  }
  out$b <- log(out$b)
  transformations <- c(transformations, "OR -> log(OR)")
 }

 # ---- Clamp extreme p-values ----
 if ("p" %in% names(out)) {
  out$p[out$p == 0] <- .Machine$double.xmin
 }

 # ---- Optional warning: SE mismatch ----
 if ("OR -> log(OR)" %in% transformations && !is.null(out$seb)) {
  warning("SE may not correspond to log(OR); verify input scale.")
 }

 # ---- Store transformation metadata ----
 if (length(transformations)) {
  attr(out, "transformations") <- transformations
 }

 # ---- Fallback marker handling ----
 if (is.null(out$marker) && "rsids" %in% colnames(out)) out$marker <- out$rsids
 if (is.null(out$marker) && "SNP" %in% colnames(out)) out$marker <- out$SNP

 # ---- Validate required fields ----
 missing <- require_fields[
  !require_fields %in% colnames(out) |
   sapply(require_fields, function(x) all(is.na(out[[x]])))
 ]
 if (length(missing)) {
  stop("Missing required canonical fields after normalization: ",
       paste(missing, collapse = ", "))
 }

 # ---- Optional column filtering ----
 if (drop_unmapped) {
  keep <- unique(c(canon, setdiff(colnames(out), names(mapping))))
  out <- out[, keep[keep %in% colnames(out)], drop = FALSE]
 }

 out
}

#' Validate normalized GWAS summary statistics
#'
#' The ingestion schema requires `chr`, `pos`, `ea`, `nea`, `b`, `seb`, and `p`.
#' Observed beta values must be finite, standard errors must be finite and
#' strictly positive, and p-values and effect-allele frequencies must be in [0, 1].
#' Optional per-variant sample sizes (`n`) must be finite and strictly positive
#' when observed. Missing sample sizes are allowed.
#' When `eaf` is absent, the writer uses effect-allele frequencies from the
#' marker catalog after alignment.
#'
#' @param stat A normalized `data.frame` (typically from [normalizeStatSchema()]).
#' @return A list with `ok`, `errors`, `warnings`, and `metrics`.
#' @export
validateStatSchema <- function(stat) {
  if (is.null(stat)) stop("`stat` is NULL")
  if (!is.data.frame(stat)) stat <- as.data.frame(stat)

  errors <- character(0)
  warnings <- character(0)

  required <- .gwas_required_fields
  miss <- required[!required %in% colnames(stat)]
  if (length(miss)) errors <- c(errors, paste("Missing required columns:", paste(miss, collapse = ", ")))
  empty <- Filter(function(field) all(is.na(stat[[field]])), intersect(required, names(stat)))
  if (length(empty)) errors <- c(errors, paste("Empty required columns:", paste(empty, collapse = ", ")))

  numeric_fields <- intersect(c("b", "seb", "p", "eaf", "n"), names(stat))
  non_numeric <- Filter(function(field) !is.numeric(stat[[field]]), numeric_fields)
  if (length(non_numeric)) errors <- c(errors, paste("Columns must be numeric:", paste(non_numeric, collapse = ", ")))

  count_invalid <- function(field, invalid) {
    values <- stat[[field]]
    if (!is.numeric(values)) return(0L)
    observed <- values[!is.na(values) | is.nan(values)]
    sum(!is.finite(observed) | invalid(observed))
  }
  bad_beta <- count_invalid("b", function(x) FALSE)
  bad_seb <- count_invalid("seb", function(x) x <= 0)
  bad_p <- count_invalid("p", function(x) x < 0 | x > 1)
  bad_eaf <- count_invalid("eaf", function(x) x < 0 | x > 1)
  bad_n <- count_invalid("n", function(x) x <= 0)
  if (bad_beta > 0) errors <- c(errors, sprintf("Non-finite beta values (b): %s", bad_beta))
  if (bad_seb > 0) errors <- c(errors, sprintf("Standard errors (seb) must be finite and strictly positive: %s", bad_seb))
  if (bad_p > 0) errors <- c(errors, sprintf("Invalid p-values outside [0,1]: %s", bad_p))
  if (bad_eaf > 0) errors <- c(errors, sprintf("Invalid effect-allele frequencies (eaf) outside [0,1]: %s", bad_eaf))
  if (bad_n > 0) errors <- c(errors, sprintf("Per-variant sample sizes (n) must be finite and strictly positive: %s", bad_n))

  if (all(c("ea", "nea") %in% colnames(stat))) {
    valid <- c("A", "C", "G", "T")
    bad_ea <- sum(!toupper(stat$ea) %in% valid, na.rm = TRUE)
    bad_nea <- sum(!toupper(stat$nea) %in% valid, na.rm = TRUE)
    if (bad_ea > 0) warnings <- c(warnings, sprintf("EA non-ACGT rows: %s", bad_ea))
    if (bad_nea > 0) warnings <- c(warnings, sprintf("NEA non-ACGT rows: %s", bad_nea))
  }

  if (all(c("chr", "pos", "ea", "nea") %in% colnames(stat))) {
    cpra <- paste(stat$chr, stat$pos, stat$ea, stat$nea, sep = "_")
    dup <- sum(duplicated(cpra))
    if (dup > 0) warnings <- c(warnings, sprintf("Duplicate chr-pos-allele rows: %s", dup))
  }

  metrics <- list(
    n_rows = nrow(stat),
    n_cols = ncol(stat),
    p_missing = if ("p" %in% colnames(stat)) mean(is.na(stat$p)) else NA_real_
  )

  list(ok = !length(errors), errors = unique(errors), warnings = unique(warnings), metrics = metrics)
}

#' Ingest GWAS summary statistics with AI-ready schema detection pipeline
#'
#' Runs schema detection, normalization, validation, and then delegates
#' storage/QC to [updateStatDB()].
#' A gact-managed download is removed after the study and its descriptor
#' are saved, using the URL and checksum in the ingestion recipe.
#'
#' @param GAlist A gact database object.
#' @param stat A raw summary-statistics `data.frame`.
#' @param schema Optional schema output from [detectStatSchema()] or named mapping.
#' @param min_confidence Minimum per-field confidence used for auto mode.
#' @param strict If `TRUE`, stop when low confidence exists for required fields.
#' @param write_to_db Passed to [updateStatDB()] as `writeStatDB`.
#' @param ... Arguments passed to [updateStatDB()], including study metadata
#'   and `strand` for documented forward-strand data.
#' @param effect_scale Passed to [normalizeStatSchema()].
#' @param p_scale Passed to [normalizeStatSchema()].
#' @return A list with updated `GAlist`, `schema`, `validation`, normalized
#'   input, transformations, the import `qc` summary, and download `cleanup`
#'   status when applicable.
#' @export
ingestStatDB <- function(GAlist,
                         stat,
                         schema = NULL,
                         min_confidence = 0.1,
                         strict = TRUE,
                         write_to_db = TRUE,
                         ...,
                         effect_scale = c("auto", "beta", "or"),
                         p_scale = c("auto", "p", "neglog10p")) {
  effect_scale <- match.arg(effect_scale)
  p_scale <- match.arg(p_scale)
  if (is.null(schema)) schema <- detectStatSchema(stat)

  mapping <- if (is.list(schema) && !is.null(schema$mapping)) schema$mapping else schema
  confidence <- if (is.list(schema) && !is.null(schema$confidence)) schema$confidence else setNames(rep(1, length(mapping)), names(mapping))

  mapped_required <- intersect(.gwas_required_fields, names(mapping))
  scores <- confidence[mapped_required]
  low_conf <- mapped_required[!is.na(unlist(mapping[mapped_required])) &
                               !is.na(scores) & scores < min_confidence]
  if (strict && length(low_conf)) {
    stop("Low-confidence required mappings: ", paste(low_conf, collapse = ", "),
         ". Provide `schema=` explicitly or set strict = FALSE.")
  }

  stat_norm <- normalizeStatSchema(stat, schema = mapping, effect_scale = effect_scale,
                                  p_scale = p_scale)
  validation <- validateStatSchema(stat_norm)
  if (!validation$ok) stop(paste(validation$errors, collapse = " | "))

  updated <- updateStatDB(
    GAlist = GAlist,
    stat = stat_norm,
    writeStatDB = write_to_db,
    ...
  )

  cleanup <- NULL
  if (isTRUE(write_to_db)) {
    args <- list(...)
    if (!is.null(args[["descriptor"]])) cleanup <- .gact_gwas_cleanup_download(args[["recipe"]])
  }
  list(
   GAlist = updated,
   schema = schema,
   validation = validation,
   qc = attr(updated, "qc_metrics"),
   normalized_stat = stat_norm,
   transformations = attr(stat_norm, "transformations"),
   cleanup = cleanup
  )
}

#' Align a summary-statistics file to the gact marker catalog
#'
#' Match by rsid, chromosome, position and allele pair. Swap the allele labels,
#' reverse beta and Z signs, and complement effect-allele frequency when needed.
#' Repair older Z signs that equal the negative of beta divided by SE.
#' @param file An existing gact-format CSV or TSV file, plain or gzipped.
#' @param GAlist A gact descriptor pointing to the marker catalog.
#' @param out_file Output path. Defaults to replacing the source after validation.
#' @return Invisibly, the output path, recode counts and the file's MD5.
#' @export
recodeStat <- function(file, GAlist, out_file = file) {
  one_path <- function(x) is.character(x) && length(x) == 1L &&
    !is.na(x) && nzchar(x)
  if (!one_path(file) || !one_path(out_file))
    stop("file and out_file must each name one file.", call. = FALSE)
  file <- normalizePath(file, winslash = "/", mustWork = TRUE)
  out_file <- path.expand(out_file)
  parent <- normalizePath(dirname(out_file), winslash = "/", mustWork = TRUE)
  out_file <- file.path(parent, basename(out_file))
  if (file.exists(out_file))
    out_file <- normalizePath(out_file, winslash = "/", mustWork = TRUE)
  if (dir.exists(file) || dir.exists(out_file) ||
      (file != out_file && file.exists(out_file)))
    stop("Provide a study file and either the same file or a new output path.", call. = FALSE)
  marker_file <- file.path(GAlist$dirs[["marker"]], "markers.txt.gz")
  if (length(marker_file) != 1L || !file.exists(marker_file))
    stop("GAlist must point to an installed markers.txt.gz.", call. = FALSE)
  marker_file <- normalizePath(marker_file, winslash = "/", mustWork = TRUE)
  if (identical(out_file, marker_file))
    stop("Cannot overwrite the marker catalog.", call. = FALSE)
  stage <- tempfile(".gact-recode-", tmpdir = parent)
  if (!dir.create(stage)) stop("Cannot create a staging directory.", call. = FALSE)
  on.exit(unlink(stage, recursive = TRUE))
  counts <- .Call(C_gact_recode_sumstat, enc2utf8(file), enc2utf8(marker_file),
    enc2utf8(file.path(stage, "output")), enc2utf8(out_file),
    grepl("\\.gz$", out_file, ignore.case = TRUE))
  invisible(list(file = out_file, reversed = counts[1L], repaired_z = counts[2L],
    checksum = unname(tools::md5sum(out_file))))
}

#' @rdname recodeStat
#' @export
recodeSumstat <- recodeStat

# Batch ingestion and C++ worker integration --------

#' Ingest explicitly selected verified GWAS recipes
#'
#' Reads recipes.tsv and uses the selected JSON recipes as the source of study
#' metadata and column mappings. No selection means no import. Existing imports
#' are checked against GWAS_information.csv before downloading inputs. Matching
#' accession, ancestry and input checksum are skipped; changed checksums or
#' recipes require review. GWAS1, GWAS2, ... are assigned locally on publication.
#' The CSV and saved descriptor are updated after successful imports. recipes.tsv
#' and JSON recipes are not used as import-status records and are not modified.
#' The native catalog reader accepts gact CSV/TSV catalogs, quoted fields and
#' integer-valued decimal/scientific coordinates. CPRA keys are checked by parsed
#' coordinates and alleles; invalid rows identify their catalog row number.
#'
#' @param GAlist Database descriptor list, or path to its saved RDS.
#' @param selection Required accessions or recipe filenames, or selected rows
#'   from recipes.tsv. A character selection must match exactly one recipe.
#'   Use recipe_file, ancestry and/or input_md5 to disambiguate data-frame rows.
#'   data.table selections are supported.
#'   Metadata edits to selected TSV rows do not override the JSON recipe.
#' @param workers Number of C++ staging workers. Windows uses one worker.
#' @param import Whether to import. FALSE checks selected recipes without
#'   downloading full inputs or changing the database.
#' @param descriptor Saved GAlist RDS path. By default find the single
#'   GAlist*.rds in the database root. The current saved descriptor is reloaded.
#' @param data_dir Optional input cache directory; defaults to download/.
#' @param out_dir Optional new report directory; defaults to recipes/reviews/.
#' @param scratch_dir Optional parent of temporary native working files.
#' @return With `import = TRUE`, invisibly returns the updated GAlist. Its
#'   `ingest_report` contains the per-study status, local IDs and QC counts for
#'   this call. With `import = FALSE`, returns the readiness report directly.
#' @export
ingestGWAS <- function(GAlist, selection, workers = 1L, import = TRUE,
                       descriptor = NULL, data_dir = NULL, out_dir = NULL,
                       scratch_dir = NULL) {
  if (missing(selection) || is.null(selection) || !length(selection) ||
      (is.data.frame(selection) && !nrow(selection)))
    stop("Select studies explicitly using accessions or rows from recipes.tsv")
  descriptor_loaded <- is.character(GAlist) && length(GAlist) == 1L
  if (descriptor_loaded) {
    if (!is.null(descriptor)) stop("Supply the descriptor path only once")
    descriptor <- GAlist
    GAlist <- readRDS(descriptor)
  }
  if (!is.list(GAlist) || length(GAlist$dirs[["gstat"]]) != 1L)
    stop("Supply a GAlist database descriptor")
  database <- dirname(normalizePath(GAlist$dirs[["gstat"]], winslash = "/", mustWork = TRUE))
  if (is.null(descriptor)) {
    candidates <- list.files(database, pattern = "^GAlist.*\\.rds$", full.names = TRUE)
    if (length(candidates) != 1L) stop("Supply descriptor: expected one saved GAlist RDS in the database root")
    descriptor <- candidates
  }
  descriptor <- normalizePath(descriptor, winslash = "/", mustWork = TRUE)
  current <- if (descriptor_loaded) GAlist else readRDS(descriptor)
  if (!identical(normalizePath(current$dirs[["gstat"]], winslash = "/", mustWork = TRUE),
                 normalizePath(GAlist$dirs[["gstat"]], winslash = "/", mustWork = TRUE)))
    stop("descriptor belongs to a different database")
  recipes <- file.path(database, "recipes")
  tsv <- file.path(recipes, "recipes.tsv")
  if (!file.exists(tsv)) stop("Run scanGWAS() to create recipes/recipes.tsv first")
  inventory <- data.table::fread(tsv, sep = "\t", colClasses = "character", data.table = FALSE)
  keys <- c("accession", "ancestry", "input_md5", "recipe_file")
  if (!all(keys %in% names(inventory))) stop("recipes.tsv is missing recipe identity columns")
  if (is.character(selection)) {
    if (anyNA(selection) || any(!nzchar(selection)) || anyDuplicated(selection))
      stop("Supply unique, nonempty study accessions or recipe filenames")
    hits <- lapply(selection, function(key)
      which(inventory$accession == key | inventory$recipe_file == key))
  } else if (is.data.frame(selection) && !anyDuplicated(names(selection))) {
    selection <- as.data.frame(selection)
    fields <- intersect(keys, names(selection))
    if (!any(c("accession", "recipe_file") %in% fields) || anyNA(selection[fields]))
      stop("Select rows using accession or recipe_file, with nonmissing identity fields")
    hits <- lapply(seq_len(nrow(selection)), function(i) {
      match <- rep(TRUE, nrow(inventory))
      for (field in fields) match <- match & inventory[[field]] == as.character(selection[[field]][i])
      which(match)
    })
  } else stop("selection must be accessions, recipe filenames or selected recipes.tsv rows")
  if (any(lengths(hits) != 1L))
    stop("Each selection must match exactly one recipe; use recipe_file and ancestry/input_md5 to disambiguate")
  indices <- unlist(hits, use.names = FALSE)
  if (anyDuplicated(indices)) stop("The same recipe was selected more than once")
  rows <- lapply(indices, function(i) {
    key <- inventory[i, , drop = FALSE]
    name <- key$recipe_file
    if (is.na(name) || !grepl("^[A-Za-z0-9][A-Za-z0-9._-]*\\.json$", name))
      stop("Invalid recipe_file in recipes.tsv")
    path <- file.path(recipes, name)
    recipe <- jsonlite::read_json(path)
    for (field in c("accession", "ancestry", "input_md5")) {
      if (!identical(recipe[[field]], key[[field]])) stop("Recipe identity differs from recipes.tsv: ", name)
    }
    if (!identical(name, basename(.gact_verified_recipe_path(recipes, recipe$accession, recipe$input_md5))))
      stop("Recipe filename does not match its accession/checksum: ", name)
    fields <- intersect(c("accession", "filename", "cohort_evidence", .gact_gwas_metadata_fields), names(recipe))
    row <- lapply(recipe[fields], function(value) {
      if (is.null(value) || !length(value)) return(NA)
      if (!is.atomic(value) || length(value) != 1L) stop("Recipe study metadata must be scalar: ", name)
      value
    })
    row$url <- recipe$input_url
    row$md5 <- recipe$input_md5
    row <- as.data.frame(row, stringsAsFactors = FALSE)
    if (is.null(.gact_read_verified_recipe(path, row))) stop("Recipe is not verified: ", name)
    if (tolower(row$ancestry) %in% c("mixed", "multiple", "unknown", "not reported"))
      stop("Select single-ancestry recipes")
    row$recipe_file <- name
    row
  })
  manifest <- as.data.frame(data.table::rbindlist(rows, fill = TRUE))
  result <- .gact_run_gwas_manifest(manifest, descriptor, GAlist = current,
    data_dir = data_dir, out_dir = out_dir,
    import = import, engine = "cpp", workers = workers, scratch_dir = scratch_dir)
  if (!import) return(invisible(result$report))
  result$GAlist$ingest_report <- result$report
  invisible(result$GAlist)
}

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
#'   Original Pan-UKB ancestry recipes may share a filename when their URL and
#'   checksum agree and the ancestries differ. The download cache reuses that file.
#'   Shared study metadata can be supplied through `defaults`. Other metadata
#'   fields from [inspectGWAS()] are optional.
#' @param descriptor Path to the existing database's GAlist RDS. Successful
#'   imports update this file. Previously imported studies are detected by
#'   accession and ancestry in GWAS_information.csv, including legacy Catalog
#'   study URLs in `reference`. Input checksum conflicts require review.
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
#'   to check readiness without importing studies. A matching verified
#'   [scanGWAS()] recipe in the database's `recipes/` directory avoids a full
#'   download in this mode. Importing still verifies the full input checksum.
#' @param preview_rows Maximum number of rows per preview.
#' @param engine Select the R importer or the C++ batch engine for supported
#'   verified gzip TSV profiles. C++ supports optional mapped per-variant `n`.
#'   Missing counts are retained and never filled from study-level metadata.
#'   Duplicate rows with different `n` conflict. The study QC reports retained
#'   rows with and without sample sizes as `n_with_sample_size` and
#'   `n_missing_sample_size`.
#' @param workers Number of C++ staging workers on Unix. Windows uses one worker.
#' @param scratch_dir Parent directory for temporary C++ working files.
#'   Defaults to SLURM_TMPDIR when available, or R's temporary directory.
#'   Each run removes its own working directory when it exits.
#' @param liftover_map Directory from [buildLiftoverMap()]. The C++ engine
#'   automatically uses the installed map for GRCh37 inputs targeting GRCh38.
#'   Supply this argument to override the installed map. With a map,
#'   inputs with a verified forward-strand profile are checked against both
#'   builds using the mapped chromosome, position and allele columns.
#'   The declared build must agree. GRCh37
#'   variants are lifted during alignment to the GRCh38 catalog using the saved
#'   index. The study QC includes lifted rows and unmapped liftover drops.
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
                            scratch_dir = NULL, liftover_map = NULL) {
  result <- .gact_run_gwas_manifest(manifest, descriptor, data_dir = data_dir,
    out_dir = out_dir, knowledge = knowledge, defaults = defaults, import = import,
    preview_rows = preview_rows, engine = match.arg(engine), workers = workers,
    scratch_dir = scratch_dir, liftover_map = liftover_map)
  invisible(result$report)
}

.gact_run_gwas_manifest <- function(manifest, descriptor, GAlist = readRDS(descriptor),
                                    data_dir = NULL, out_dir = NULL, knowledge = NULL,
                                    defaults = list(), import = TRUE, preview_rows = 1000L,
                                    engine, workers = 1L, scratch_dir = NULL, liftover_map = NULL) {
  if (!is.null(liftover_map) && engine != "cpp") stop("liftover_map requires engine = 'cpp'")
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
  build <- .gact_catalog_build(GAlist)
  map_info <- NULL
  gstat <- GAlist$dirs[["gstat"]]
  if (length(gstat) != 1L || is.na(gstat) || !dir.exists(gstat)) {
    stop("GAlist must point to an existing gstat directory")
  }
  inventory <- file.path(gstat, "GWAS_information.csv")
  imported <- data.frame()
  inventory_ids <- if (file.exists(inventory)) {
    imported <- read.csv2(inventory, colClasses = "character")
    ids <- imported$id
    if (is.null(ids)) stop("GWAS_information.csv is missing the id column")
    ids
  } else character()
  study_ids <- as.character(GAlist$study$id)
  if (anyDuplicated(study_ids) || anyDuplicated(inventory_ids) ||
      !setequal(study_ids, as.character(inventory_ids))) {
    stop("GAlist study IDs differ from GWAS_information.csv. Load the current descriptor.")
  }
  for (field in intersect(c("accession", "reference", "ancestry", "input_md5", "recipe_checksum"),
                          intersect(names(imported), names(GAlist$study)))) {
    expected <- as.character(GAlist$study[[field]][match(imported$id, study_ids)])
    actual <- imported[[field]]
    same <- (is.na(expected) & is.na(actual)) | (!is.na(expected) & !is.na(actual) & expected == actual)
    if (!all(same)) stop("GAlist import identity differs from GWAS_information.csv: ", field)
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
  repeated_file <- duplicated(filenames) | duplicated(filenames, fromLast = TRUE)
  identity <- paste(accessions, toupper(manifest$ancestry))
  duplicate <- duplicated(identity) | duplicated(identity, fromLast = TRUE) |
    (repeated_file & !.gact_shared_panukb_inputs(manifest))
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
      hits <- .gact_gwas_inventory_match(imported, row, strict = "recipe_file" %in% names(row))
      if (length(hits) > 1L) stop("Accession matches multiple imported studies")
      if (length(hits)) {
        study_id <- imported$id[hits]
        study_file <- unname(GAlist$studyfiles[study_id])
        if (length(study_file) != 1L || !file.exists(study_file)) {
          stop("Stored summary statistics are missing for ", study_id)
        }
        old_md5 <- imported$input_md5[hits]
        if (length(old_md5) && !is.na(old_md5) &&
            !identical(tolower(old_md5), tolower(row$md5))) {
          stop("Accession was imported with a different input checksum")
        }
        old_recipe <- imported$recipe_checksum[hits]
        if (length(old_recipe) && !is.na(old_recipe) && "recipe_file" %in% names(row)) {
          actual <- unname(tools::md5sum(file.path(recipes, basename(row$recipe_file))))
          if (is.na(actual) || !identical(actual, old_recipe)) stop("Imported study has a different verified recipe; review before retrying")
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
      needs_lift <- !identical(build, row$build)
      map_dir <- liftover_map
      if (is.null(map_dir) && needs_lift && engine == "cpp" && identical(build, "GRCh38")) {
        installed <- .gact_liftover_dir(GAlist)
        if (dir.exists(installed)) map_dir <- installed
      }
      if (needs_lift && (is.null(map_dir) || !identical(build, "GRCh38") ||
                         !identical(row$build, "GRCh37"))) {
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
      path <- .gact_verified_recipe_path(recipes, accession, row$md5)
      recipe <- .gact_read_verified_recipe(path, row)
      if (is.null(recipe) && "recipe_file" %in% names(row))
        stop("Selected verified recipe is missing or unverified; rescan before ingestion")
      inspected <- is.null(recipe)
      if (inspected) {
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
        issues <- recipe$inspection_issues
        if (!isTRUE(entry$preview_ok) || !isTRUE(recipe$inspection$profile_reused) || length(issues)) {
          stop(paste(c("Review the preview or verified format", issues), collapse = " | "))
        }
        recipe$accession <- accession
        recipe$input_md5 <- tolower(row$md5)
        recipe$input_url <- row$url
        if (engine == "cpp") .gact_cpp_gwas_mapping(recipe)
      } else {
        report$recipe[i] <- path
        if (import) {
          file <- .gact_gwas_manifest_file(row, data_dir)
          recipe$file_name <- file
        }
      }
      forward <- identical(recipe$strand, "forward")
      if (!is.null(map_dir) && (needs_lift || forward)) {
        if (!forward || !identical(build, "GRCh38"))
          stop("Liftover requires a GRCh38 catalog and a verified forward-strand profile")
        if (is.null(map_info)) map_info <- .gact_liftover_index(map_dir, GAlist)
        check <- recipe$build_check
        if (is.null(check)) {
          if (!inspected && !import) file <- .gact_gwas_manifest_file(row, data_dir)
          preview <- .gact_cpp_gwas_preview(file, preview_rows, 0L)
          mapping <- unlist(recipe$mapping)
          preview <- stats::setNames(preview[unname(mapping[c("chr", "pos", "nea", "ea")])],
                                    .gact_liftover_stat_keys)
          check <- verifyGWASBuild(preview, map_dir)
          check$index_md5 <- map_info$md5
          recipe$build_check <- check
          if (inspected) jsonlite::write_json(recipe, entry$recipe, pretty = TRUE, auto_unbox = TRUE)
        }
        if (!identical(check$index_md5, map_info$md5) || !identical(check$build, row$build)) {
          stop("Preview does not verify the declared build ", row$build, " with this map (",
               check$grch37_only, " GRCh37-only, ", check$grch38_only, " GRCh38-only)")
        }
        recipe$filename <- row$filename
        recipe$verified <- TRUE
      }
      report$warnings[i] <- paste(unlist(recipe$validation$warnings), collapse = " | ")
      if (!import) {
        report$status[i] <- "ready"
        next
      }
      if (engine == "cpp") {
        jobs[[i]] <- list(recipe = recipe, file = file,
                          lift_index = if (needs_lift) map_info$path else NULL)
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
    batch <- .gact_cpp_gwas_batch(GAlist, descriptor, jobs, report, recipes, work_dir, as.integer(workers))
    GAlist <- batch$GAlist
    report <- batch$report
  }
  data.table::fwrite(report, file.path(out_dir, "manifest_review.tsv"), sep = "\t", na = "NA")
  attr(report, "out_dir") <- out_dir
  list(GAlist = GAlist, report = report)
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

.gact_gwas_inventory_match <- function(study, row, strict = FALSE) {
  hits <- .gact_gwas_manifest_matches(study, toupper(row$accession))
  if (!length(hits)) return(hits)
  ancestry <- toupper(trimws(as.character(study$ancestry[hits])))
  if (length(ancestry) != length(hits)) ancestry <- rep(NA_character_, length(hits))
  unknown <- is.na(ancestry) | !nzchar(ancestry) | ancestry %in% c("UNKNOWN", "NOT REPORTED")
  if (strict && any(unknown)) stop("Stored ancestry is missing; review the existing import")
  hits <- hits[unknown | ancestry == toupper(trimws(row$ancestry))]
  if (strict && length(hits)) {
    checksum <- study$input_md5[hits]
    if (length(checksum) != length(hits) || anyNA(checksum) || any(!nzchar(checksum)))
      stop("Stored input checksum is missing; review the existing import")
  }
  hits
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

.gact_cpp_gwas_preview <- function(file, rows, skip, prefix = FALSE) {
  if (!grepl("\\.tsv\\.(gz|bgz)$", file, ignore.case = TRUE) || skip != 0L ||
      rows >= .Machine$integer.max) {
    stop("The C++ reader expects a gzip TSV with a header on the first line")
  }
  text <- .Call(C_gact_gwas_preview, file, as.integer(rows), prefix)
  if (length(strsplit(text, "\n", fixed = TRUE)[[1L]]) < 2L) stop("Preview has no data rows")
  data.table::fread(text = text, sep = "\t", nrows = rows, data.table = FALSE, check.names = FALSE)
}

.gact_cpp_gwas_mapping <- function(recipe) {
  fields <- c("chr", "pos", "ea", "nea", "b", "seb", "p", "eaf", "n")
  mapping <- unlist(recipe$mapping)
  mapping <- mapping[!is.na(mapping)]
  header <- unlist(recipe$inspection$headers, use.names = FALSE)
  if (!all(fields[1:7] %in% names(mapping)) ||
      !all(names(mapping) %in% c("marker", fields)) ||
      anyDuplicated(header) || any(!nzchar(header)) ||
      any(grepl('[,"\t\r\n]', header)) ||
      any(c("ncase", "ncontrol", "info", "z", "z_score") %in% header) ||
      !identical(recipe$effect_scale, "beta") ||
      !isTRUE(recipe$p_scale %in% c("p", "neglog10p")) ||
      !is.null(recipe[["n"]])) {
    stop("Verified recipe is outside the C++ profile. Use the R engine for this study")
  }
  if (length(setdiff(.gact_gwas_n_columns(header), mapping["n"])) ||
      ("n" %in% names(mapping) && !mapping[["n"]] %in% header))
    stop("Per-variant n columns require an explicit, unambiguous sample-size mapping")
  unname(as.character(mapping[fields]))
}

.gact_liftover_index <- function(map_dir, ga) {
  manifest <- .gact_liftover_map_manifest(map_dir)
  index <- manifest$index
  if (!identical(index$file, "index.bin") || !isTRUE(grepl("^[a-f0-9]{32}$", index$md5)))
    stop("Use buildLiftoverMap() to prepare the reusable liftover index")
  catalog <- file.path(ga$dirs[["marker"]], "markers.txt.gz")
  if (!identical(unname(tools::md5sum(catalog)), manifest$inputs$marker_catalog$md5))
    stop("Liftover map was built for a different GRCh38 marker catalog")
  path <- normalizePath(file.path(map_dir, index$file), winslash = "/", mustWork = TRUE)
  if (!identical(unname(tools::md5sum(path)), index$md5)) stop("Liftover index checksum mismatch")
  list(path = path, md5 = index$md5)
}

.gact_gwas_work_dir <- function(parent) {
  if (is.null(parent)) {
    parent <- Sys.getenv("SLURM_TMPDIR", unset = "")
    if (!nzchar(parent) || !dir.exists(parent)) parent <- tempdir()
  }
  if (!is.character(parent) || length(parent) != 1L || is.na(parent) || !nzchar(parent)) {
    stop("scratch_dir must name a directory")
  }
  parent <- path.expand(parent)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create scratch_dir: ", parent)
  }
  work <- tempfile("gact-gwas-", tmpdir = normalizePath(parent, winslash = "/", mustWork = TRUE))
  if (!dir.create(work, mode = "0700")) stop("Cannot create a GWAS work directory in ", parent)
  work
}

.gact_cpp_gwas_batch <- function(ga, descriptor, jobs, report, recipe_dir, work_dir, workers) {
  marker <- ga$dirs[["marker"]]
  markers <- file.path(marker, "markers.txt.gz")
  cpra <- file.path(marker, "cpra.txt.gz")
  stopifnot(file.exists(markers), file.exists(cpra))
  index <- file.path(work_dir, "marker-index.bin")
  files <- c(markers, cpra)
  fingerprint <- tools::md5sum(files)
  ids <- data.table::fread(markers, select = "rsids")$rsids
  if (!identical(as.character(ga$rsids), as.character(ids))) stop("Marker catalog and GAlist disagree")
  rm(ids)
  .Call(C_gact_gwas_index, markers, cpra, index)
  dir.create(recipe_dir, showWarnings = FALSE)
  if (!dir.exists(recipe_dir)) stop("Cannot create database recipe directory")
  active <- which(lengths(jobs) > 0L)
  for (i in active) {
    recipe <- jobs[[i]]$recipe
    path <- file.path(recipe_dir, paste0(recipe$accession, "-", recipe$input_md5, ".json"))
    if (!file.exists(path)) {
      temp <- tempfile(tmpdir = recipe_dir)
      tryCatch({
        jsonlite::write_json(recipe, temp, auto_unbox = TRUE, pretty = TRUE)
        if (!file.rename(temp, path)) stop("Cannot save verified recipe: ", path)
      }, finally = unlink(temp))
    }
    stored <- jsonlite::read_json(path)
    fields <- c("accession", "input_md5", "mapping", .gact_gwas_metadata_fields)
    if (!identical(stored[fields], recipe[fields]) ||
        (!is.null(stored$build_check) && !identical(stored$build_check, recipe$build_check)) ||
        !identical(stored$inspection$headers, recipe$inspection$headers)) {
      stop("Existing recipe disagrees with reviewed recipe: ", path)
    }
    report$recipe[i] <- path
    jobs[[i]]$recipe_path <- path
  }

  stage_one <- function(i) {
    job <- jobs[[i]]
    tryCatch({
      base <- file.path(work_dir, sprintf("study_%05d", i))
      file <- paste0(base, ".txt.gz")
      qc <- paste0(base, ".qc.json")
      recipe <- jsonlite::read_json(job$recipe_path)
      if (!identical(unname(tools::md5sum(job$file)), recipe$input_md5)) stop("Input checksum changed before staging")
      .Call(C_gact_gwas_stage, index, job$file, file, qc,
            unlist(recipe$inspection$headers, use.names = FALSE),
            .gact_cpp_gwas_mapping(recipe), identical(recipe$strand, "forward"),
            job$lift_index, identical(recipe$p_scale, "neglog10p"))
      metrics <- jsonlite::read_json(qc, simplifyVector = TRUE)
      checksum <- unname(tools::md5sum(file))
      stopifnot(!is.na(checksum), metrics$n_final > 0,
                metrics$n_input == metrics$n_final + metrics$n_dropped,
                metrics$n_dropped == sum(unlist(metrics$dropped_by_reason)),
                metrics$n_flipped <= metrics$n_final)
      c(job, list(staged = file, qc_file = qc, qc = metrics, checksum = checksum))
    }, error = function(e) list(error = conditionMessage(e)))
  }
  count <- min(workers, length(active))
  staged <- if (.Platform$OS.type == "unix" && count > 1L) {
    parallel::mclapply(active, stage_one, mc.cores = count, mc.preschedule = FALSE)
  } else lapply(active, stage_one)
  passed <- vapply(staged, function(value) is.list(value) && !is.null(value$staged), logical(1))
  for (j in which(!passed)) {
    value <- staged[[j]]
    report$issues[active[j]] <- if (is.list(value) && length(value$error)) value$error else
      "C++ worker exited before producing a validated result"
  }
  if (!any(passed)) return(list(GAlist = ga, report = report))
  if (!identical(fingerprint, tools::md5sum(files))) stop("Marker catalog changed during staging")
  result <- tryCatch(.gact_publish_cpp_gwas(ga, descriptor, staged[passed]),
                     error = function(e) e)
  if (inherits(result, "error")) {
    report$issues[active[passed]] <- conditionMessage(result)
    return(list(GAlist = ga, report = report))
  }
  for (j in seq_along(result$ids)) {
    i <- active[passed][j]
    item <- staged[passed][[j]]
    report$status[i] <- "imported"
    report$study_id[i] <- result$ids[j]
    for (field in c("n_input", "n_final", "n_dropped", "n_flipped")) {
      report[[field]][i] <- item$qc[[field]]
    }
    cleanup <- .gact_gwas_cleanup_download(item$recipe)
    if (!is.null(cleanup)) {
      report$cleanup[i] <- cleanup$status
      if (length(cleanup$warning)) {
        warnings <- c(report$warnings[i], cleanup$warning)
        report$warnings[i] <- paste(warnings[!is.na(warnings) & nzchar(warnings)], collapse = " | ")
      }
    }
  }
  list(GAlist = result$GAlist, report = report)
}

.gact_publish_cpp_gwas <- function(ga, descriptor, staged) {
  gstat <- ga$dirs[["gstat"]]
  inventory <- file.path(gstat, "GWAS_information.csv")
  stage <- file.path(gstat, ".gact-import.lock")
  if (!dir.create(stage, mode = "0700", showWarnings = FALSE)) stop("Another GWAS import holds the database lock")
  temp_rds <- tempfile(".gact-descriptor-", tmpdir = dirname(descriptor))
  backup <- file.path(stage, "previous-inventory.csv")
  published <- character()
  replaced <- committed <- FALSE
  on.exit({
    if (!committed) {
      if (replaced) {
        if (file.exists(backup)) {
          if (!file.rename(backup, inventory)) warning("Cannot restore the previous GWAS inventory")
        } else unlink(inventory)
      }
      unlink(published)
    }
    unlink(temp_rds)
    unlink(stage, recursive = TRUE)
  }, add = TRUE)
  current <- readRDS(descriptor)
  if (!identical(current$study, ga$study) || !identical(current$dirs, ga$dirs) ||
      !identical(current$rsids, ga$rsids)) stop("Database changed while studies were staging. Retry the batch")
  ids <- as.character(current$study$id)
  if (file.exists(inventory)) {
    if (!setequal(read.csv2(inventory, colClasses = "character")$id, ids) ||
        !file.copy(inventory, backup)) stop("Cannot verify and back up the current GWAS inventory")
  } else if (length(ids)) stop("GWAS inventory disappeared while studies were staging")
  last <- max(c(0L, suppressWarnings(as.integer(sub("^GWAS", "", ids)))), na.rm = TRUE)
  new_ids <- paste0("GWAS", last + seq_along(staged))
  source <- destination <- character()
  for (j in seq_along(staged)) {
    item <- staged[[j]]
    recipe <- item$recipe
    if (!is.null(item$lift_index)) recipe$build <- "GRCh38"
    id <- new_ids[j]
    at <- length(current$study$id) + 1L
    current$study <- lapply(current$study, function(column) { length(column) <- at; column })
    recipeID <- sub("^GWAS", "REC", id)
    name <- paste0(id, ".txt")
    qc <- paste0(id, ".qc.json")
    recipe_path <- file.path("..", "recipes", basename(item$recipe_path))
    counts <- c(n = "study_n_total", ncase = "study_n_case", ncontrol = "study_n_control")
    numbers <- lapply(names(counts), function(field) {
      value <- recipe[[field]]
      if (!length(value)) value <- recipe[[counts[[field]]]]
      if (!length(value)) NA_real_ else as.numeric(value)
    })
    names(numbers) <- c("n", "ncase", "ncontrol")
    total <- numbers$ncase + numbers$ncontrol
    neff <- if (identical(recipe$type, "quantitative")) numbers$n else
      if (is.finite(total) && total > 0) numbers$ncase * numbers$ncontrol / total else NA_real_
    row <- c(list(id = id, file = name, created_at = Sys.time(), recipeID = recipeID,
                  recipe_file = recipe_path, qc_file = qc, checksum = item$checksum,
                  accession = recipe$accession, input_md5 = recipe$input_md5,
                  input_url = recipe$input_url,
                  recipe_checksum = unname(tools::md5sum(item$recipe_path)), neff = neff),
             numbers, Filter(function(value) length(value) > 0L, recipe[intersect(names(recipe), c("source", "trait", "type", "category",
               "efo", "gender", "ancestry", "build", "reference", "comments",
               "study_n_total", "study_n_case", "study_n_control"))]))
    for (field in names(row)) {
      column <- current$study[[field]]
      if (is.null(column)) column <- rep(row[[field]][NA_integer_], at)
      column[at] <- row[[field]]
      current$study[[field]] <- column
    }
    for (pair in list(c(item$staged, paste0(name, ".gz")), c(item$qc_file, qc))) {
      from <- file.path(stage, pair[2L])
      if (!file.copy(pair[1L], from)) stop("Cannot stage ", pair[2L], " for publication")
      source <- c(source, from)
      destination <- c(destination, file.path(gstat, pair[2L]))
    }
  }
  current$studies <- as.data.frame(current$study)
  rownames(current$studies) <- current$studies$id
  current$studyfiles <- stats::setNames(file.path(gstat, paste0(current$study$file, ".gz")), current$study$id)
  temp_csv <- file.path(stage, basename(inventory))
  write.csv2(current$studies, temp_csv, row.names = FALSE)
  saveRDS(current, temp_rds, compress = FALSE)
  suspendInterrupts({
    for (i in seq_along(source)) {
      .publish_gwas_import_file(source[i], destination[i])
      published <- c(published, destination[i])
    }
    .publish_gwas_import_file(temp_csv, inventory, replace = TRUE)
    replaced <- TRUE
    .publish_gwas_import_file(temp_rds, descriptor, replace = TRUE)
    committed <- TRUE
  })
  list(GAlist = current, ids = new_ids)
}

# GWAS database quality control and study updates --------

#' @export
#'
#' @title Remove studies from a GACT database
#'
#' @description Deletes the selected summary-statistics files and their QC reports and updates the database study metadata.
#'
#' @param GAlist GACT database list containing directory paths, study metadata and mappings.
#' @param studyID Character vector of existing study IDs to remove.
#'
#' @return Updated GAlist. The summary files are removed and existing GWAS_information.csv metadata is rewritten.
#' @export
removeStatDB <- function(GAlist,studyID=NULL) {
 if (is.null(studyID)) stop("Please provide a valid studyID.")
 if (!all(studyID %in% GAlist$study$id))
  stop("One or more study IDs not found in the database.")
 # remove summary statistics files
 for (i in 1:length(studyID)) {
  file.remove(GAlist$studyfiles[studyID[i]])
 }
 qc_files <- GAlist$study$qc_file[match(studyID, GAlist$study$id)]
 if (length(qc_files)) {
  qc_files <- qc_files[!is.na(qc_files) & nzchar(qc_files)]
  qc_paths <- file.path(GAlist$dirs["gstat"], basename(qc_files))
  file.remove(qc_paths[file.exists(qc_paths)])
 }
 attr(GAlist, "qc_metrics") <- NULL
 # update GAlist
 #study <- as.data.frame(GAlist$study)
 #study <- study[!study$id%in%studyID,]
 #GAlist$study <- as.list(study)
 #studyfiles <- GAlist$studyfiles
 #studyfiles <- studyfiles[!names(studyfiles)%in%studyID]
 #GAlist$studyfiles <- studyfiles

 # Update GAlist metadata
 study <- as.data.frame(GAlist$study)
 study <- study[!study$id%in%studyID,]
 GAlist$study <- as.list(study)

 GAlist$studyfiles <- GAlist$studyfiles[!names(GAlist$studyfiles) %in% studyID]
 GAlist$studies <- study
 rownames(GAlist$studies) <- GAlist$studies$id

 # Update on-disk metadata
 file_stat_information <- file.path(GAlist$dirs["gstat"], "GWAS_information.csv")
 if (file.exists(file_stat_information))
  write.csv2(study, file_stat_information, row.names = FALSE)

 message("Removed studies: ", paste(studyID, collapse = ", "))
 return(GAlist)
}




.publish_gwas_import_file <- function(from, to, replace = FALSE) {
 if (!replace && file.exists(to)) stop("GWAS import destination already exists: ", to)
 if (!file.rename(from, to)) stop("Cannot publish GWAS import file: ", to)
}

#' Update and process genetic association study summary statistics
#'
#' This function updates and processes genetic association study summary statistics. The input to the function includes study information, summary statistics, and information about the sample.
#' Study files and metadata are staged before publication. If publication
#' fails, files published by that attempt are removed and the previous
#' inventory is preserved. Use the latest descriptor to add further studies.
#'
#' @param GAlist A list object providing information and infrastructure of the gact database.
#' @param stat A data frame with `chr`, `pos`, `ea`, `nea`, `b`, `seb`, and `p`
#'   when writing a study. Supply `b` and `seb` on the beta or log-odds scale.
#'   If `eaf` is omitted, QC uses the marker catalog's effect-allele frequency.
#' @param source Source of the summary statistics
#' @param trait Trait being studied
#' @param type Type of trait (binary or quantitative)
#' @param category Category of trait (complex, biomarker)
#' @param efo Experimental Factor Ontology (complex, biomarker)
#' @param gender Gender of study participants
#' @param ancestry Ancestry of study participants (EUR, AFR, ,,,,, Mixed)
#' @param build Genome build used in the study. For a marker catalog with a
#'   known build, a conflicting GRCh37 or GRCh38 study build is rejected.
#' @param n Study sample size. Used for variant `n` when the column is absent.
#' @param ncase Number of cases
#' @param ncontrol Number of controls
#' @param reference Reference for the study
#' @param comments Additional comments
#' @param recipe Optional ingestion recipe to save alongside the processed study.
#' @param writeStatDB A flag indicating whether to perform quality control and write the processed summary statistics to a file
#' @param strand Input allele strand. Use `"forward"` only when the source
#'   documents forward-strand alleles for the declared build. The default
#'   `"unknown"` drops strand-ambiguous A/T and C/G SNPs.
#' @param descriptor Optional path to the current GAlist RDS. When supplied, the
#'   updated descriptor is published with the study and inventory.
#' @return A list of updated study information, including the study id, file name, trait, type, gender, ancestry, build, sample size, case/control count, reference, and comments.
#'   The QC summary is attached as `qc_metrics` and saved as `GWAS<n>.qc.json`.
#'
#' @export
#'
updateStatDB <- function(GAlist=NULL,
                         stat=NULL,
                         source="unknown",
                         trait="unknown",
                         type="unknown",
                         category="unknown",
                         efo="unknown",
                         gender="unknown",
                         ancestry="unknown",
                         build="unknown",
                         n=NULL,
                         ncase=NULL,
                         ncontrol=NULL,
                         reference="unknown",
                         comments="none",
                         recipe=NULL,
                         writeStatDB=TRUE,
                         strand = c("unknown", "forward"),
                         descriptor = NULL) {

 strand <- match.arg(strand)
 attr(GAlist, "qc_metrics") <- NULL

 # Update study information only - useful if we need to add extra information
 message("Collecting information on external summary statistics")
 study_number <- length(GAlist$study$id)+1

 #studyID <- paste0("GWAS",study_number)

 # Extract numeric part of all current IDs - works if some studies have been removed
 existing_ids <- GAlist$study$id
 existing_nums <- as.integer(gsub("\\D", "", existing_ids))
 next_id_num <- ifelse(length(existing_nums) == 0, 1, max(existing_nums, na.rm = TRUE) + 1)
 studyID <- paste0("GWAS", next_id_num)
 recipeID <- paste0("REC", next_id_num)

 # Extend existing optional metadata as well as the standard fields. This
 # preserves a rectangular study table when successive studies omit fields.
 GAlist$study <- lapply(GAlist$study, function(column) {
  length(column) <- study_number
  column
 })

 GAlist$study$id[study_number] <- studyID
 GAlist$study$file[study_number] <- paste0(studyID,".txt")
 GAlist$study$trait[study_number] <- trait
 GAlist$study$type[study_number] <- type
 GAlist$study$category[study_number] <- category
 GAlist$study$efo[study_number] <- efo
 GAlist$study$gender[study_number] <- gender
 GAlist$study$ancestry[study_number] <- ancestry
 GAlist$study$build[study_number] <- build
 GAlist$study$n[study_number] <- if (is.null(n)) NA_real_ else n
 GAlist$study$ncase[study_number] <- if (is.null(ncase)) NA_real_ else ncase
 GAlist$study$ncontrol[study_number] <- if (is.null(ncontrol)) NA_real_ else ncontrol
 GAlist$study$comments[study_number] <- comments

 # After studyID creation
 GAlist$study$created_at[study_number] <- Sys.time()
 # Optional recipe hook
 GAlist$study$recipeID[study_number] <- recipeID

 #if(type=="binary") {
 # ntotal <- ncase + ncontrol
 # pcase <- ncase/(ncase+ncontrol)
 # GAlist$study$neff[study_number] <- ntotal*pcase*(1-pcase)
 #}
 #if(type=="quantitative") GAlist$study$neff[study_number] <- n
 if(type == "binary" && !is.null(ncase) && !is.null(ncontrol) &&
    !is.na(ncase) && !is.na(ncontrol)) {

  ntotal <- ncase + ncontrol
  pcase <- ncase / ntotal
  GAlist$study$neff[study_number] <- ntotal * pcase * (1 - pcase)

 } else if(type == "quantitative") {

  GAlist$study$neff[study_number] <- if (is.null(n)) NA_real_ else n

 } else {
  GAlist$study$neff[study_number] <- NA
 }
 GAlist$study$reference[study_number] <- reference
 GAlist$study$source[study_number] <- source
 if (!is.null(recipe$accession)) GAlist$study$accession[study_number] <- recipe$accession
 if (!is.null(recipe$input_md5)) GAlist$study$input_md5[study_number] <- recipe$input_md5

 gstat_dir <- unname(GAlist$dirs["gstat"])
 if (length(gstat_dir) != 1L || !dir.exists(gstat_dir)) {
  stop("GAlist must point to an existing gstat directory.")
 }
 file_stat_information <- file.path(gstat_dir, "GWAS_information.csv")
 GAlist$studyfiles <- file.path(gstat_dir, paste0(GAlist$study$file, ".gz"))
 names(GAlist$studyfiles) <- GAlist$study$id


 # processing of summary statistics
 if(writeStatDB) {
  message("Perform quality control of external summary statistics")

  catalog_build <- .gact_catalog_build(GAlist)
  builds <- c(grch37 = "GRCh37", hg19 = "GRCh37", grch38 = "GRCh38", hg38 = "GRCh38")
  study_build <- unname(builds[tolower(trimws(build))])
  if (strand == "forward" && (length(build) != 1L || is.na(build) ||
      tolower(trimws(build)) %in% c("", "unknown"))) {
   stop("Declare the GWAS build when using forward-strand alleles.")
  }
  if (length(study_build) == 1L && !is.na(study_build) &&
      length(catalog_build) == 1L && !is.na(catalog_build) && study_build != catalog_build) {
   stop("GWAS build ", study_build, " differs from marker catalog build ", catalog_build, ".")
  }
  if (!is.data.frame(stat)) stop("stat must be a data.frame")
  stat <- as.data.frame(stat)
  if ("p" %in% names(stat)) stat$p <- as.numeric(stat$p)
  validation <- validateStatSchema(stat)
  if (!validation$ok) stop(paste(validation$errors, collapse = " | "))
  if (!("n" %in% names(stat)) && length(n) == 1L && is.finite(n)) {
   stat$n <- n
  }

  stat <- qcStatDB(GAlist = GAlist, stat = stat, excludeMAFDIFF = 0.5, strand = strand)
  qc_metrics <- attr(stat, "qc_metrics")
  if (!nrow(stat)) {
   dropped <- unlist(qc_metrics$dropped_by_reason)
   dropped <- dropped[dropped > 0L]
   stop("No variants passed GWAS quality control. Dropped: ",
        paste(paste(names(dropped), dropped, sep = "="), collapse = ", "), ".")
  }
 }

 # A single staging directory also serializes inventory updates.
 stage <- file.path(gstat_dir, ".gact-import.lock")
 if (!dir.create(stage, showWarnings = FALSE, mode = "0700")) {
  stop("Cannot stage GWAS import at ", stage,
       ". Check for another import or a filesystem permission error.")
 }
 published <- character()
 committed <- FALSE
 staged_descriptor <- inventory_backup <- NULL
 inventory_published <- FALSE
 on.exit({
  if (!committed && inventory_published) {
   restored <- if (is.null(inventory_backup)) {
    unlink(file_stat_information) == 0L
   } else {
    file.rename(inventory_backup, file_stat_information)
   }
   if (!restored) warning("Could not restore GWAS inventory: ", file_stat_information)
  }
  if (!committed) unlink(published)
  if (!is.null(staged_descriptor)) unlink(staged_descriptor)
  unlink(stage, recursive = TRUE)
  if (!committed && any(file.exists(published))) {
   warning("Could not remove GWAS files from the failed import: ",
           paste(published[file.exists(published)], collapse = ", "))
  }
 }, add = TRUE)

 if (file.exists(file_stat_information)) {
  inventory <- read.csv2(file_stat_information, colClasses = "character")
  if (!"id" %in% names(inventory) || !setequal(inventory$id, as.character(existing_ids))) {
   stop("GAlist study IDs differ from GWAS_information.csv. Load the current descriptor before importing.")
  }
 }

 recipe_path <- qc_path <- NULL
 if (writeStatDB) {
  message(paste("Writing processed summary statistics to internal file:",
                GAlist$study$file[study_number]))
  file_stat <- unname(GAlist$studyfiles[study_number])
  if (file.exists(file_stat)) stop("GWAS summary statistics file already exists: ", file_stat)
  staged_stat <- file.path(stage, basename(file_stat))
  fwrite(stat, staged_stat, col.names = TRUE)
  if (!is.null(qc_metrics)) {
   qc_path <- file.path(gstat_dir, paste0(studyID, ".qc.json"))
   staged_qc <- file.path(stage, basename(qc_path))
   jsonlite::write_json(qc_metrics, staged_qc, auto_unbox = TRUE, pretty = TRUE)
   GAlist$study$qc_file[study_number] <- basename(qc_path)
   attr(GAlist, "qc_metrics") <- qc_metrics
  }
  if (!is.null(recipe)) {
   recipe$qc <- qc_metrics
   recipe_path <- file.path(gstat_dir, paste0(recipeID, ".json"))
   if (file.exists(recipe_path)) stop("GWAS recipe file already exists: ", recipe_path)
   staged_recipe <- file.path(stage, basename(recipe_path))
   jsonlite::write_json(recipe, staged_recipe, auto_unbox = TRUE, pretty = TRUE)
   GAlist$study$recipe_file[study_number] <- basename(recipe_path)
  }
  checksum <- unname(tools::md5sum(staged_stat))
  if (is.na(checksum)) stop("Cannot checksum staged GWAS file: ", staged_stat)
  GAlist$study$checksum[study_number] <- checksum
 }
 GAlist$studies <- as.data.frame(GAlist$study)
 rownames(GAlist$studies) <- GAlist$studies$id
 staged_inventory <- file.path(stage, basename(file_stat_information))
 write.csv2(GAlist$studies, file = staged_inventory, row.names = FALSE)
 if (!is.null(descriptor)) {
  descriptor <- normalizePath(descriptor, mustWork = TRUE)
  staged_descriptor <- tempfile(".gact-descriptor-", tmpdir = dirname(descriptor))
  saveRDS(GAlist, staged_descriptor, compress = FALSE)
  if (file.exists(file_stat_information)) {
   inventory_backup <- file.path(stage, "previous-inventory.csv")
   if (!file.copy(file_stat_information, inventory_backup)) {
    stop("Cannot back up the GWAS inventory: ", file_stat_information)
   }
  }
 }

 suspendInterrupts({
  if (writeStatDB) {
   .publish_gwas_import_file(staged_stat, file_stat)
   published <- c(published, file_stat)
  }
  if (!is.null(recipe_path)) {
   .publish_gwas_import_file(staged_recipe, recipe_path)
   published <- c(published, recipe_path)
  }
  if (!is.null(qc_path)) {
   .publish_gwas_import_file(staged_qc, qc_path)
   published <- c(published, qc_path)
  }
  .publish_gwas_import_file(staged_inventory, file_stat_information, replace = TRUE)
  inventory_published <- TRUE
  if (!is.null(staged_descriptor)) {
   .publish_gwas_import_file(staged_descriptor, descriptor, replace = TRUE)
  }
  committed <- TRUE
 })
 return(GAlist)
}


#' @export
#'
# qcStatDB <- function(GAlist=NULL, stat=NULL, excludeMAF=0.01, excludeMAFDIFF=0.05,
#                      excludeINFO=0.8, excludeCGAT=TRUE, excludeINDEL=TRUE,
#                      excludeDUPS=TRUE, excludeMHC=FALSE, excludeMISS=0.05,
#                      excludeHWE=1e-12) {
#
#  # we use cpra to link sumstats and Glist
#  cpra <- as.vector(fread(file = file.path(GAlist$dirs["marker"], "cpra.txt.gz"), header=FALSE, data.table=FALSE))[[1]]
#  if(cpra[1]=="cpra") cpra <- cpra[-1]    # this is a fix
#  rsids <- GAlist$rsids
#
#  # stat is a data.frame
#  if(!is.data.frame(stat)) stop("stat should be  a data frame")
#  if(!is.null(stat$rsids)) rownames(stat) <- stat$rsids
#
#  if (!"chr" %in% names(stat) || !"pos" %in% names(stat))
#   stop("Input `stat` must contain 'chr' and 'pos' columns for mapping.")
#
#
#  # Remove rows with missing values
#  stat <- na.omit(stat)
#
#  # Handle zero values in p column
#  stat$p <- as.numeric(stat$p)
#  stat$p[stat$p == 0] <- .Machine$double.xmin
#
#
#  # Internal summary statistic column format
#  # data.frame(rsids, chr, pos, a1, a2, af, b, seb, stat, p, n)     (single trait)
#  # list(marker=(rsids, chr, pos, a1, a2, af), b, seb, stat, p, n)  (multiple trait)
#
#  fm_internal <- c("marker","chr","pos","ea","nea","eaf","b","seb","p","n")
#
#  format <- "unknown"
#
#  if(all(fm_internal%in%colnames(stat))) format <- "internal"
#  #if(all(fm_internal[1:9]%in%colnames(stat))) format <- "internal"
#  if(all(fm_internal[1:5]%in%colnames(stat))) format <- "internal"
#
#  if(format=="unknown") {
#   message("Column headings for stat object not found")
#   message("Column headings for stat object should be:")
#   print(fm_internal)
#   stop("please revised your stat object accordingly")
#  }
#
#  cpra1 <- paste(stat[,"chr"],stat[,"pos"],toupper(stat[,"ea"]),toupper(stat[,"nea"]),sep="_")
#  cpra2 <- paste(stat[,"chr"],stat[,"pos"],toupper(stat[,"nea"]),toupper(stat[,"ea"]),sep="_")
#
#  mapped <- cpra1%in%cpra | cpra2%in%cpra
#  message("Map markers based on cpra")
#  message(paste("Number of markers in stat mapped to marker ids in GAlist:",sum(mapped)))
#  message(paste("Number of markers in stat not mapped to marker ids in GAlist:",sum(!mapped)))
#
#  stat <- stat[mapped,]
#  cpra1 <- cpra1[mapped]
#  cpra2 <- cpra2[mapped]
#  rws1 <- match(cpra1,cpra)
#  rws2 <- match(cpra2,cpra)
#
#  stat$marker[!is.na(rws1)] <- rsids[rws1[!is.na(rws1)]]
#  stat$marker[!is.na(rws2)] <- rsids[rws2[!is.na(rws2)]]
#
#  isdup <- duplicated(stat$marker)
#  if(any(isdup)) message("Removing markers with duplicated ids")
#  if(any(isdup)) message(paste("Number of markers duplicated in stat:",sum(isdup)))
#  stat <- stat[!isdup,]
#  rownames(stat) <- stat$marker
#
#  marker <- fread(file.path(GAlist$dirs["marker"],"markers.txt.gz"),
#                  data.table=FALSE)
#
#  rownames(marker) <- marker$rsids
#
#  #message("Filtering markers based on information in GAlist:")
#  #message("")
#
#  #message("Filtering markers based on information in stat:")
#  #message("")
#
#  if(!is.null(stat$marker)) marker_in_stat <- marker$rsids%in%stat$marker
#  message(paste("Number of markers in stat also found in bimfiles:", sum(marker_in_stat)))
#  message("")
#  if(sum(marker_in_stat)==0) stop("No marker ids found in bimfiles")
#
#  # align marker and stat object
#  marker <- marker[marker_in_stat,]
#  stat <- stat[marker$rsids,]
#
#  aligned <- stat$ea==marker$ea
#  message(paste("Number of effect alleles aligned with first allele in bimfiles:", sum(aligned)))
#  message(paste("Number of effect alleles not aligned with first allele in bimfiles:", sum(!aligned)))
#  message("")
#
#  #original
#  effect <- stat[,"b"]
#  effect_allele <- stat[,"ea"]
#  non_effect_allele <- stat[,"nea"]
#  if(!is.null(stat$eaf)) effect_allele_freq <- stat[,"eaf"]
#
#  # flip is not aligned
#  stat[!aligned,"b"] <- -effect[!aligned]
#  stat[!aligned,"ea"] <- non_effect_allele[!aligned]
#  stat[!aligned,"nea"] <- effect_allele[!aligned]
#  if(!is.null(stat$eaf)) stat[!aligned,"eaf"] <- 1-effect_allele_freq[!aligned]
#
#  if(is.null(stat$eaf)) {
#   message("No effect allele frequency (eaf) provided - using eaf in database")
#   stat$eaf <- marker$eaf
#  }
#
#
#  # exclude based on maf
#  excludeMAFDIFF <- abs(marker$eaf-stat$eaf) > excludeMAFDIFF
#  message(paste("Number of markers excluded by large difference between MAF difference:", sum(excludeMAFDIFF)))
#  message("")
#  stat <- stat[!excludeMAFDIFF,]
#  marker <- marker[!excludeMAFDIFF,]
#
#  if(is.null(stat$n)) stat$n <- neff(seb=stat$seb,af=stat$eaf)
#  colnames(stat)[1] <- "rsids"
#  return(stat)
# }
#' @title Map and align summary statistics to a GACT marker reference
#'
#' @description Maps chromosome-position-allele keys to the database, removes duplicate markers, aligns effect alleles and filters allele-frequency differences.
#'
#' @param GAlist GACT database list containing directory paths, study metadata and mappings.
#' @param stat Ordinary data frame with chr, pos, ea, nea and p columns; b is used when alleles are flipped. Optional eaf, seb and n columns supply frequencies, standard errors and sample sizes.
#' @param excludeMAFDIFF Maximum absolute difference between reported and reference effect-allele frequency.
#' @param excludeMAF,excludeINFO,excludeCGAT,excludeINDEL,excludeDUPS,excludeMHC,excludeMISS,excludeHWE Legacy compatibility arguments, currently unused.
#' @param return_metrics Attach the QC summary as the `qc_metrics` attribute.
#' @param strand Input allele strand. Use `"forward"` only for documented
#'   forward-strand input on the catalog build. The default `"unknown"`
#'   drops strand-ambiguous A/T and C/G SNPs.
#'
#' @details Matches either allele order against the marker catalog and keeps
#'   uniquely matched allele pairs. Unsplit allele lists and shared rsIDs are
#'   dropped. Duplicate estimates are compared after alignment, rounding EAF
#'   to 12 decimal places for that comparison. Conflicting groups are dropped.
#'   An absent EAF column is filled from the catalog. Output follows catalog
#'   order. The summary counts each dropped row once and counts flips only
#'   among retained rows.
#'
#' @return Aligned, filtered data frame with marker IDs and optional qc_metrics.
#' @export
qcStatDB <- function(GAlist = NULL,
                     stat = NULL,
                     excludeMAF = 0.01,
                     excludeMAFDIFF = 0.05,
                     excludeINFO = 0.8,
                     excludeCGAT = TRUE,
                     excludeINDEL = TRUE,
                     excludeDUPS = TRUE,
                     excludeMHC = FALSE,
                     excludeMISS = 0.05,
                     excludeHWE = 1e-12,
                     return_metrics = TRUE,
                     strand = c("unknown", "forward")) {
 strand <- match.arg(strand)
 if (!is.data.frame(stat)) stop("stat must be a data.frame")
 stat <- as.data.frame(stat)
 required <- c("chr", "pos", "ea", "nea", "p")
 missing <- setdiff(required, names(stat))
 if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

 cpra <- data.table::fread(
  file.path(GAlist$dirs["marker"], "cpra.txt.gz"),
  header = FALSE
 )[[1]]
 if (identical(cpra[1L], "cpra")) cpra <- cpra[-1L]
 catalog <- data.table::fread(
  file.path(GAlist$dirs["marker"], "markers.txt.gz"),
  select = c("rsids", "ea", "nea", "eaf"),
  data.table = FALSE
 )
 if (length(cpra) != nrow(catalog) ||
     !identical(as.character(GAlist$rsids), as.character(catalog$rsids))) {
  stop("The marker catalog and GAlist are out of sync.")
 }

 # Give each rejected row its first reason, in processing order.
 n_input <- nrow(stat)
 drop <- rep.int("", n_input)
 drop[!complete.cases(stat[required])] <- "missing_required"
 ea <- toupper(stat$ea)
 nea <- toupper(stat$nea)
 unsplit <- grepl(",", ea, fixed = TRUE) | grepl(",", nea, fixed = TRUE)
 drop[drop == "" & unsplit] <- "unsplit_multiallelic"
 if (strand == "unknown") {
  palindromic <- paste0(ea, nea) %in% c("AT", "TA", "CG", "GC") &
   nchar(ea) == 1L & nchar(nea) == 1L
  drop[drop == "" & palindromic] <- "unresolved_strand"
 }

 rows <- which(drop == "")
 key <- paste(stat$chr[rows], stat$pos[rows], ea[rows], nea[rows], sep = "_")
 swapped_key <- paste(stat$chr[rows], stat$pos[rows], nea[rows], ea[rows], sep = "_")
 direct_match <- match(key, cpra)
 swapped_match <- match(swapped_key, cpra)
 target <- pmin(direct_match, swapped_match, na.rm = TRUE)
 repeated_keys <- cpra[duplicated(cpra)]
 ambiguous <- (!is.na(direct_match) & !is.na(swapped_match)) |
  key %in% repeated_keys | swapped_key %in% repeated_keys
 drop[rows[is.na(target)]] <- "unmatched"
 drop[rows[ambiguous]] <- "ambiguous_catalog_match"

 repeated_ids <- catalog$rsids[duplicated(catalog$rsids)]
 shared_id <- catalog$rsids[target] %in% repeated_ids
 drop[rows[drop[rows] == "" & shared_id]] <- "shared_catalog_rsid"
 mapped <- drop[rows] == ""
 rows <- rows[mapped]
 target <- target[mapped]
 flipped <- !is.na(swapped_match[mapped])
 stat <- stat[rows, , drop = FALSE]
 marker <- catalog[target, , drop = FALSE]
 stat$p <- as.numeric(stat$p)
 stat$p[stat$p == 0] <- .Machine$double.xmin

 stat$marker <- stat$rsids <- marker$rsids
 stat$ea <- marker$ea
 stat$nea <- marker$nea
 for (field in intersect(c("b", "z", "z_score"), names(stat))) {
  stat[[field]][flipped] <- -as.numeric(stat[[field]][flipped])
 }
 if (!("eaf" %in% names(stat))) {
  stat$eaf <- marker$eaf
 } else {
  stat$eaf[flipped] <- 1 - stat$eaf[flipped]
 }

 duplicate_rows <- which(duplicated(target) | duplicated(target, fromLast = TRUE))
 if (length(duplicate_rows)) {
  fields <- intersect(c("rsids", "b", "z", "z_score", "seb", "p", "eaf", "n", "ncase",
                         "ncontrol", "info"), names(stat))
  estimates <- stat[duplicate_rows, fields, drop = FALSE]
  estimates$eaf <- round(estimates$eaf, 12L)
  distinct_ids <- unique(estimates)$rsids
  conflicting <- stat$rsids[duplicate_rows] %in% distinct_ids[duplicated(distinct_ids)]
  drop[rows[duplicate_rows[conflicting]]] <- "conflicting_duplicate"

  # Prefer an already aligned copy, then sort its fields to break ties.
  copies <- duplicate_rows[!conflicting]
  fields <- lapply(stat[copies, , drop = FALSE], function(x) {
   if (is.list(x)) as.character(x) else x
  })
  copies <- copies[do.call(order, c(list(flipped[copies]), unname(fields), list(method = "radix")))]
  drop[rows[copies[duplicated(target[copies])]]] <- "redundant_duplicate"
 }

 maf_diff <- abs(marker$eaf - stat$eaf)
 drop[rows[drop[rows] == "" & !is.na(maf_diff) & maf_diff > excludeMAFDIFF]] <- "eaf_mismatch"
 keep <- which(drop[rows] == "")
 keep <- keep[order(target[keep])]
 stat <- stat[keep, , drop = FALSE]
 rownames(stat) <- stat$rsids

 reasons <- c("missing_required", "unsplit_multiallelic", "unresolved_strand", "unmatched",
              "ambiguous_catalog_match", "shared_catalog_rsid", "conflicting_duplicate",
              "redundant_duplicate", "eaf_mismatch")
 counts <- setNames(tabulate(match(drop, reasons), length(reasons)), reasons)
 metrics <- list(
  strand = strand,
  n_input = n_input,
  n_mapped = length(rows),
  n_final = nrow(stat),
  n_dropped = sum(counts),
  n_flipped = sum(flipped[keep]),
  n_removed_maf = unname(counts["eaf_mismatch"]),
  dropped_by_reason = as.list(counts)
 )
 message("GWAS QC: ", metrics$n_final, " of ", n_input, " rows retained, ",
         metrics$n_dropped, " dropped, ", metrics$n_flipped, " retained rows flipped.")
 if (return_metrics) {
  attr(stat, "qc_metrics") <- metrics
 }
 stat
}



#' @export
#'
#' @title Suggest summary-statistic columns using legacy heuristics
#'
#' @description Inspects column classes, cardinalities and numeric ranges to suggest column names.
#'
#' @param stat Data frame of summary-statistic columns.
#'
#' @details The heuristics are intended for inspection and may fail on inputs without the expected column types. They do not validate or normalize a schema; use detectStatSchema(), normalizeStatSchema() and validateStatSchema() for the current ingestion workflow.
#'
#' @return Character vector of candidate column names; the vector is not a named schema mapping.
#' @export
columnStatDB <- function(stat=NULL) {
 col_names <- colnames(stat)
 types <- sapply(stat,class)
 isInteger <- types=="integer"
 isNumeric <- types=="numeric"
 isCharacter <- types=="character"
 isChr <- isPos <- isMarker <- NULL

 if(any(isInteger)) {
  nlevs <- sapply(stat[,isInteger],function(x){length(unique(x))})
  maxVal <- sapply(stat[,isInteger],function(x){max(unique(x))})
  isChr <- col_names[isInteger][22 <= nlevs & nlevs <= 24]
  isPos <- col_names[isInteger][maxVal>100000000]
 }
 if(any(isCharacter)) {
  nlevs <- sapply(stat[,isCharacter],function(x){length(unique(x))})
  if(is.null(isChr)) isChr <- col_names[isCharacter][22 <= nlevs & nlevs <= 24]
  #isAllele <- col_names[isCharacter][nlevs<5]
  isAllele <- sapply(stat[,isCharacter],function(x){sum(c("A","T","C","G")%in%toupper(unique(x)))==4})
  isAllele <- col_names[isCharacter][isAllele]
  if(length(isAllele)==1) isEA <- isNEA <- isAllele
  if(length(isAllele)==2) {
   isEA <- isAllele[1]
   isNEA <- isAllele[2]
  }
  isMarker <- col_names[isCharacter][nlevs>10000]
 }
 if(any(isNumeric)) {
  minVal <- sapply(stat[,isNumeric],min, na.rm=TRUE)
  maxVal <- sapply(stat[,isNumeric],max, na.rm=TRUE)
  meanVal <- sapply(stat[,isNumeric],mean, na.rm=TRUE)
  medianVal <- sapply(stat[,isNumeric],median, na.rm=TRUE)
  isB <- col_names[isNumeric][minVal<0]
  isSEB <- col_names[isNumeric][0<minVal & maxVal<0.2]
  isP <- col_names[isNumeric][meanVal>0 & minVal<0.00001 & meanVal<0.4]
  isN <- col_names[isNumeric][minVal>1000]
  isNca <- isN[1]
  isNco <- isN[1]
  isN <- isN[1]
  if(max(maxVal)>100000000) isPos <- col_names[isNumeric][maxVal>100000000 && minVal>0]
  isEAF <- col_names[isNumeric][0<=minVal & maxVal<1 & medianVal>0.1 & medianVal<0.7]
  if(length(isEAF)==0) isEAF <- "Unknown"
  isInfo <- col_names[isNumeric][0.3<minVal & maxVal<=1 & medianVal>0.7]
  if(length(isInfo)==0) isInfo <- "Unknown"
 }
 best_guess <- c(isMarker,isChr,isPos, isEA, isNEA, isEAF,
                 isB, isSEB, isP, isN, isNca, isNco, isInfo)
 #names(best_guess) <- c("marker","chr","pos","ea","nea","eaf",
 #                      "b","seb","p","n","ncase","ncontrols","info")
 best_guess <- as.character(na.omit(best_guess))
 best_guess
}
