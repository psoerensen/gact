test_that("Pan-UKB scanning writes importable recipes from bounded previews", {
  skip_if_not_installed("yaml")
  root <- tempfile("scan-gwas-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  phenotype <- data.frame(filename = c("bio.tsv.bgz", "diabetes.tsv.bgz"),
    trait_type = c("biomarkers", "phecode"), phenocode = c("30750", "250"),
    pheno_sex = "both_sexes", description = c("Glycated haemoglobin (HbA1c)", "Diabetes mellitus"))
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
  # Include an unfinished record so the native prefix reader must discard it.
  known <- charToRaw(paste0(paste(profile$headers, collapse = "\t"), "\n",
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
    .gact_panukb_catalog_row = function(row, id, analysis, metadata_dir) {
      if (id == "GCST90691576") gact:::.gact_panukb_mismatch("ancestry differs")
      data.frame(accession = id, filename = paste0(id, ".h.tsv.gz"),
        url = paste0("https://ftp.ebi.ac.uk/", id), md5 = paste(rep("a", 32L), collapse = ""),
        source = "verified", trait = row$description, type = "quantitative",
        build = "GRCh38", ancestry = "Mixed", strand = "forward")
    },
    .gact_fetch_gwas_prefix = function(url, path, max_bytes) {
      expect_equal(max_bytes, 2048)
      calls <<- c(calls, url)
      writeBin(gzip, path)
      length(gzip)
    }, .package = "gact")
  result <- scanGWAS(data.frame(filename = phenotype$filename, analysis = "hq"),
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
  expect_length(list.files(file.path(root, "scan"), recursive = TRUE, pattern = "\\.gz$"), 0L)
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
