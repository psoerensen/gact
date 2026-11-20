test_that("Pan-UKB scanning writes importable recipes from bounded previews", {
  skip_if_not_installed("yaml")
  root <- tempfile("scan-gwas-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  phenotype <- data.frame(filename = c("bio.tsv.bgz", "diabetes.tsv.bgz"),
    trait_type = c("biomarkers", "phecode"), phenocode = c("30750", "250"),
    pheno_sex = "both_sexes", pops = "EUR",
    description = c("Glycated haemoglobin (HbA1c)", "Diabetes mellitus"))
  manifest <- file.path(root, "phenotypes.tsv.gz")
  data.table::fwrite(phenotype, manifest, sep = "\t")
  catalog <- function(id, trait) list(accessionId = id, reportedTrait = trait,
    summaryStatistics = "https://ftp.ebi.ac.uk/example")
  records <- list(catalog("GCST90692780", "Glycated hemoglobin HbA1c levels"),
    catalog("GCST90691576", "Glycated hemoglobin HbA1c levels"),
    catalog("GCST90692999", "Diabetes mellitus (PheCode 250)"),
    catalog("GCST90692416", "Diabetes mellitus (PheCode 250)"),
    catalog("GCST90692998", "Diabetes mellitus (PheCode 249)"))
  profile <- gact:::.gact_panukb_catalog_profile("verified")
  fields <- c("1", "10", "A", "C", ".2", ".1", ".2", ".05", "rs1", "lo", "10", "1_10_C_A")
  # Include the provider quality flag and an unfinished trailing record.
  fields <- c(fields, "False")
  known <- charToRaw(paste0(paste(c(profile$headers, "low_confidence"), collapse = "\t"), "\n",
    paste(fields, collapse = "\t"), "\n", paste(fields, collapse = "\t"), "\npartial"))
  compressed <- file.path(root, "sample.tsv.gz")
  con <- gzfile(compressed, "wb")
  writeBin(known, con)
  close(con)
  gzip <- readBin(compressed, "raw", n = file.info(compressed)$size - 4L)
  writeBin(gzip, compressed)
  expect_equal(nrow(gact:::.gact_cpp_gwas_preview(compressed, 3L, 0L, prefix = TRUE)), 2L)
  calls <- character()
  local_mocked_bindings(
    .gact_panukb_catalog_studies = function() records,
    .gact_panukb_catalog_row = function(row, id, analysis, metadata_dir, ...) {
      if (id == "GCST90691576") gact:::.gact_panukb_mismatch("ancestry differs")
      data.frame(accession = id, filename = paste0(id, ".h.tsv.gz"),
        url = paste0("https://ftp.ebi.ac.uk/", id), md5 = paste(rep("a", 32L), collapse = ""),
        source = "verified", trait = row$description, type = "quantitative",
        build = "GRCh38", ancestry = "EUR", strand = "forward",
        reference = paste0("https://www.ebi.ac.uk/gwas/studies/", id))
    },
    .gact_fetch_gwas_prefix = function(url, path, max_bytes) {
      expect_equal(max_bytes, 2048)
      calls <<- c(calls, url)
      writeBin(gzip, path)
      length(gzip)
    }, .package = "gact")
  expect_error(scanGWAS(data.frame(filename = phenotype$filename, analysis = "hq"),
                       file.path(root, "combined"), manifest), "single-ancestry")
  result <- scanGWAS(data.frame(filename = phenotype$filename, analysis = "EUR"),
    file.path(root, "scan"), manifest, preview_rows = 2L, max_preview_bytes = 2048L)
  expect_identical(result$review$status, c("ready", "needs_review"))
  expect_match(result$review$issues[2L], "Multiple verified Catalog matches")
  expect_identical(calls, "https://ftp.ebi.ac.uk/GCST90692780")
  expect_identical(data.table::fread(result$manifest)$accession, "GCST90692780")
  book <- jsonlite::read_json(result$knowledge)
  expect_true(book$formats[[1L]]$verified)
  recipe <- jsonlite::read_json(result$review$recipe[1L])
  expect_equal(recipe$inspection$preview_rows, 2)
  expect_true(recipe$inspection$remote)
  expect_null(recipe$file_name)
  expect_identical(recipe$input_md5, paste(rep("a", 32L), collapse = ""))
  expect_true(recipe$verified)
  expect_true("low_confidence" %in% unlist(recipe$inspection$unmapped_columns))
  expect_length(list.files(file.path(root, "scan"), recursive = TRUE, pattern = "\\.gz$"), 0L)

  recipes <- file.path(root, "scan", "recipes")
  saved <- tools::md5sum(result$review$recipe[1L])
  resumed <- scanGWAS(phenotype_manifest = manifest, recipe_dir = recipes,
                      preview_rows = 2L, max_preview_bytes = 2048L)
  expect_identical(resumed$review$status, c("ready", "needs_review"))
  expect_true(resumed$review$reused[1L])
  expect_identical(resumed$review$preview_bytes[1L], 0L)
  expect_identical(tools::md5sum(resumed$review$recipe[1L]), saved)
  expect_length(calls, 1L)

  ga <- createDB("hsa.GRCh38.0.0.1", root)
  descriptor <- file.path(ga$dbdir, "GAlist.rds")
  saveRDS(ga, descriptor)
  expect_true(file.copy(result$review$recipe[1L], ga$dirs[["recipes"]]))
  report <- runGWASManifest(resumed$manifest, descriptor, import = FALSE, engine = "cpp")
  expect_identical(report$status, "ready")
  expect_length(list.files(ga$dirs[["download"]]), 0L)
  expect_length(list.files(ga$dirs[["gstat"]]), 0L)
  changed <- data.table::fread(resumed$manifest, data.table = FALSE)
  changed$ancestry <- "CSA"
  conflict <- runGWASManifest(changed, descriptor, import = FALSE, engine = "cpp")
  expect_identical(conflict$status, "needs_review")
  expect_match(conflict$issues, "metadata changed")
})

test_that("Pan-UKB matching resolves shared labels using code, sex and ancestry", {
  manifest <- data.frame(filename = c("diabetes-a", "diabetes-b", "diabetes-c"),
    trait_type = "phecode", phenocode = c("250", "250.2", "250"),
    description = "Diabetes mellitus", pheno_sex = c("both_sexes", "both_sexes", "females"),
    pops = "CSA,EUR", n_cases_CSA = 1695, n_controls_CSA = 7181,
    n_cases_EUR = 24348, n_controls_EUR = 396181, category = "Disease")
  ids <- c("GCST90692416", "GCST90691951", "GCST90691952", "GCST90692999")
  labels <- c("Diabetes (PheCode 250)", "Diabetes mellitus (PheCode 250)",
              "Diabetes mellitus (PheCode 250.2)", "Diabetes mellitus (PheCode 250)")
  records <- lapply(seq_along(ids), function(i) list(accessionId = ids[i],
    reportedTrait = labels[i], summaryStatistics = "https://ftp.ebi.ac.uk/example",
    discoverySampleAncestry = if (i %in% c(1L, 4L)) "8876 South Asian, Central Asian" else
      "420529 European"))
  catalog <- gact:::.gact_panukb_index_catalog(records)
  row <- manifest[1L, , drop = FALSE]
  csa <- gact:::.gact_panukb_candidates(row, catalog, "CSA")
  eur <- gact:::.gact_panukb_candidates(row, catalog, "EUR")
  expect_setequal(csa$accession, ids[c(1L, 4L)])
  expect_identical(eur$accession, ids[2L])
  expect_identical(gact:::.gact_panukb_catalog_code(
    c("Body mass index (UKB data field 21001)", "Diabetes (ICD-10 code E11.9)")),
    c("ukb:21001", "icd10:e119"))

  local_mocked_bindings(.gact_catalog_file = function(id, metadata_dir) {
    i <- match(id, ids)
    file <- paste0(id, ".h.tsv.gz")
    hash <- "5d140be767245495c51f255ec6effc4b"
    list(file = file, url = paste0("https://ftp.ebi.ac.uk/", file),
      base = "https://ftp.ebi.ac.uk/", fetch = function(...) paste(hash, file),
      meta = list(gwas_id = id, data_file_name = file, genome_assembly = "GRCh38",
        coordinate_system = "1-based", is_harmonised = TRUE, file_type = "GWAS-SSF v1.0",
        author_notes = "https://pan.ukbb.broadinstitute.org/", trait_description = labels[i],
        data_file_md5sum = hash, sex = if (i == 4L) "female" else "combined",
        samples = list(list(
          sample_size = if (i %in% c(1L, 4L)) 8876 else 420529, case_control_study = TRUE,
          sample_ancestry_category = if (i %in% c(1L, 4L)) c("South Asian", "Central Asian") else
            "European"))))
  }, .package = "gact")
  match_row <- gact:::.gact_panukb_scan_match
  a <- match_row(row, csa, "CSA", manifest[2L, , drop = FALSE])
  b <- match_row(row, eur, "EUR", manifest[2L, , drop = FALSE])
  expect_identical(c(a$accession, b$accession), ids[1:2])
  expect_identical(c(a$ancestry, b$ancestry), c("CSA", "EUR"))
  expect_identical(a$md5, b$md5)
  tied <- transform(row, filename = "other-diabetes-file")
  expect_error(match_row(row, csa, "CSA", tied), "also match other-diabetes-file")
})

test_that("Catalog scanning verifies metadata and saves a checksum-keyed recipe", {
  skip_if_not_installed("yaml")
  root <- tempfile("scan-catalog-")
  on.exit(unlink(root, recursive = TRUE))
  ids <- c("GCST90270927", "GCST90018877", "GCST90270928", "GCST90270929")
  md5 <- strrep("a", 32L)
  calls <- 0L
  local_mocked_bindings(
    .gact_catalog_file = function(id, catalog_path = NULL) {
      file <- if (is.null(catalog_path)) paste0(id, ".h.tsv.gz") else basename(catalog_path)
      hash <- paste0(substr(md5, 1L, 31L), match(id, ids))
      list(file = file, url = paste0("https://ftp.ebi.ac.uk/", file), base = "https://ftp.ebi.ac.uk/",
           meta = list(gwas_id = id, data_file_name = file, genome_assembly = "GRCh38",
             coordinate_system = "1-based", is_harmonised = TRUE,
             file_type = if (id == ids[2L]) "pre-GWAS-SSF" else "GWAS-SSF v1.0",
             author_notes = "Genome-wide association summary statistics",
             analysis_software = if (id == ids[3L]) "METAL" else "REGENIE",
             trait_description = list("Gastric ulcer"), data_file_md5sum = hash, sex = "combined",
             samples = list(list(case_control_study = TRUE, sample_size = 261647,
                                 sample_ancestry_category = list("East Asian")))),
           fetch = function(...) paste(hash, file))
    },
    .gact_fetch_gwas_prefix = function(url, path, max_bytes) {
      calls <<- calls + 1L
      headers <- gact:::.gact_panukb_catalog_profile("Catalog")$headers[9:1]
      row <- c("1", "10", "A", "C", ".2", ".1", ".2", ".05", "rs1")[9:1]
      con <- gzfile(path, "wt")
      writeLines(c(paste(headers, collapse = "\t"), paste(row, collapse = "\t")), con)
      close(con)
      file.info(path)$size
    }, .package = "gact")
  evidence <- stats::setNames(c("Reviewed single Japanese cohort", "Cannot override METAL"),
                             ids[c(1L, 3L)])
  ga <- createDB("hsa.GRCh38.0.0.1", root)
  scan <- function() scanGWAS(ids, GAlist = ga, repository = "gwas-catalog",
                              preview_rows = 1L, cohort_evidence = evidence)
  result <- scan()
  expect_identical(result$review$status, rep("ready", 4))
  expect_identical(dirname(result$review$recipe[1L]), ga$dirs[["recipes"]])
  expect_length(list.files(ga$dirs[["gstat"]]), 0L)
  expect_true(is.na(result$review$issues[2L]))
  expect_identical(jsonlite::read_json(result$review$recipe[3L])$study_design, "meta_analysis")
  expect_true(is.na(result$review$issues[4L]))
  recipe <- jsonlite::read_json(result$review$recipe[1L])
  expect_equal(recipe$study_n_total, 261647)
  expect_identical(recipe$ancestry, "EAS")
  expect_null(recipe$mapping[["n"]])
  expect_identical(calls, 4L)
  expect_true(scan()$review$reused[1L])
  expect_identical(calls, 4L)
  md5 <- strrep("b", 32L)
  updated <- scan()
  expect_false(updated$review$reused[1L])
  expect_false(identical(updated$review$recipe[1L], result$review$recipe[1L]))
  expect_identical(calls, 8L)

  index <- gact:::.gact_parse_catalog_index(c(
    "./GCST90270001-GCST90271000/GCST90270927/harmonised/study.h.tsv.gz",
    "./GCST90270001-GCST90271000/GCST90270928/harmonised/first.h.tsv.gz",
    "./GCST90270001-GCST90271000/GCST90270928/harmonised/second.h.tsv.gz",
    "./unexpected/path.txt"))
  expect_identical(index$ambiguous, c(FALSE, TRUE))
  local_mocked_bindings(.gact_catalog_index = function() index, .package = "gact")
  md5 <- strrep("c", 32L)
  z <- scanGWAS(GAlist = ga, repository = "gwas-catalog",
                cohort_evidence = evidence, preview_rows = 1L)
  expect_identical(z$review$status, c("ready", "needs_review"))
  expect_identical(data.table::fread(z$manifest)$filename, "study.h.tsv.gz")
  expect_match(z$review$issues[2L], "Multiple harmonised files")
})

test_that("parallel scanning reports failures in selection order", {
  skip_if_not_installed("yaml")
  root <- tempfile("scan-parallel-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  phenotype <- file.path(root, "phenotypes.tsv.gz")
  data.table::fwrite(data.frame(filename = "known.tsv.bgz", description = "Known",
    pheno_sex = "both_sexes"), phenotype, sep = "\t")
  ga <- createDB("hsa.GRCh38.0.0.1", root)
  selection <- data.frame(filename = c("b.tsv.bgz", "a.tsv.bgz"), analysis = "EUR",
                          accession = c("GCST90691576", "GCST90692416"))
  scan <- scanGWAS(selection, GAlist = ga, repository = "pan-ukb", workers = 2L,
                   phenotype_manifest = phenotype, preview_rows = 1L)
  expect_identical(scan$review$filename, selection$filename)
  expect_identical(scan$review$status, rep("needs_review", 2L))
  expect_match(scan$review$issues, "absent from the Pan-UKB phenotype manifest", all = TRUE)
  expect_length(list.files(ga$dirs[["recipes"]], pattern = "\\.json$"), 0L)
})

test_that("preview downloads require a bounded HTTP 206 response", {
  response <- 206L
  bytes <- 4L
  local_mocked_bindings(curl_fetch_stream = function(url, fun, handle) {
    fun(as.raw(c(31L, 139L, 8L, 0L, 0L)[seq_len(bytes)]))
    list(status_code = response,
      headers = charToRaw("HTTP/2 206\r\nContent-Range: bytes 0-3/400\r\n\r\n"))
  }, .package = "curl")
  path <- tempfile()
  on.exit(unlink(path))
  fetch <- function() gact:::.gact_fetch_gwas_prefix("https://ftp.ebi.ac.uk/example", path, 4L)
  expect_equal(fetch(), 4)
  response <- 200L
  expect_error(fetch(), "valid bounded byte range")
  bytes <- 5L
  expect_error(fetch(), "exceeded the preview byte limit")
})
