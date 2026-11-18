test_that("original Pan-UKB scans save reusable single-ancestry liftover recipes", {
  root <- tempfile("original-scan-")
  ga <- createDB("hsa.GRCh38.0.0.1", root)
  on.exit(unlink(root, recursive = TRUE))
  marker <- file.path(ga$dirs[["marker"]], "markers.txt.gz")
  catalog <- data.frame(rsids = paste0("rs", 1:12), chr = 1L, pos = 101:112,
                        ea = "C", nea = "A", eaf = .2)
  data.table::fwrite(catalog, marker)
  data.table::fwrite(data.frame(cpra = with(catalog, paste(chr, pos, ea, nea, sep = "_"))),
                     file.path(dirname(marker), "cpra.txt.gz"), col.names = FALSE)
  jsonlite::write_json(list(build = "GRCh38", effect_allele = "ALT",
    checksums = list("markers.txt.gz" = unname(tools::md5sum(marker)))),
    file.path(dirname(marker), "marker_manifest.json"), auto_unbox = TRUE)
  ga$rsids <- catalog$rsids
  descriptor <- file.path(root, "GAlist.rds")
  saveRDS(ga, descriptor)
  stat <- data.frame(chr = 1L, pos = 1:12, ref = "A", alt = "C",
    beta_EUR = .2, se_EUR = .1, neglog10_pval_EUR = 3, af_EUR = .2,
    beta_CSA = .3, se_CSA = .1, neglog10_pval_CSA = 4, af_CSA = .2,
    beta_meta = 99, low_confidence_EUR = FALSE)
  variants <- file.path(root, "variants.tsv")
  data.table::fwrite(stat, variants, sep = "\t")
  chain <- file.path(root, "example.chain")
  writeLines(c("chain 100 chr1 200 + 0 12 chr1 200 + 100 112 1", "12"), chain)
  map <- buildLiftoverMap(ga, variants, chain)
  index_before <- tools::md5sum(file.path(map, "index.bin"))
  filenames <- c("biomarkers-30750-both_sexes-irnt.tsv.bgz",
                 "phecode-250-both_sexes.tsv.bgz", "wrong-build.tsv.bgz",
                 "missing-se.tsv.bgz")
  upstream <- file.path(root, filenames)
  data.table::fwrite(stat, upstream[1L], sep = "\t")
  binary <- stat[setdiff(names(stat), c("af_EUR", "af_CSA"))]
  binary$af_cases_EUR <- .1
  binary$af_controls_EUR <- .2
  data.table::fwrite(binary, upstream[2L], sep = "\t")
  data.table::fwrite(transform(stat, pos = pos + 100L), upstream[3L], sep = "\t")
  data.table::fwrite(stat[setdiff(names(stat), "se_EUR")], upstream[4L], sep = "\t")
  phenotype <- data.frame(filename = filenames,
    trait_type = c("biomarkers", "phecode", "biomarkers", "biomarkers"),
    phenocode = c("30750", "250", "other", "other2"), pheno_sex = "both_sexes",
    pops = c("EUR,CSA", "EUR", "EUR", "EUR"),
    description = c("HbA1c", "Diabetes mellitus", "Wrong build", "Missing SE"),
    n_cases_EUR = c(1000, 20, 1000, 1000), n_controls_EUR = c(NA, 980, NA, NA),
    n_cases_CSA = c(100, NA, NA, NA), n_controls_CSA = NA,
    aws_path = paste0("s3://pan-ukb-us-east-1/sumstats_flat_files/", filenames),
    md5_hex = unname(tools::md5sum(upstream)))
  metadata <- file.path(root, "phenotypes.tsv.gz")
  data.table::fwrite(phenotype, metadata, sep = "\t", na = "NA", quote = FALSE)
  requests <- character()
  local_mocked_bindings(
    .gact_panukb_catalog_studies = function(...) stop("Unexpected Catalog discovery"),
    .gact_fetch_gwas_prefix = function(url, path, max_bytes) {
      requests <<- c(requests, url)
      stopifnot(file.copy(upstream[match(basename(url), filenames)], path))
      file.size(path)
    }, .package = "gact")
  scan <- function(...) scanGWAS(GAlist = ga, repository = "pan-ukb",
    files = "original", phenotype_manifest = metadata, preview_rows = 12L, ...)
  result <- scan()
  expect_identical(result$review$status, c("ready", "ready", "ready", "needs_review", "needs_review"))
  expect_match(result$review$issues[4L], "does not verify")
  expect_match(result$review$issues[5L], "Mapped|required|mapping", ignore.case = TRUE)
  manifest <- data.table::fread(result$manifest, data.table = FALSE)
  expect_identical(manifest$ancestry, c("EUR", "CSA", "EUR"))
  expect_identical(manifest$build, rep("GRCh37", 3))
  expect_equal(manifest$study_n_total, c(1000, 100, 1000))
  expect_equal(manifest$study_n_case[3L], 20)
  expect_equal(manifest$study_n_control[3L], 980)
  expect_identical(manifest$md5[1:2], rep(phenotype$md5_hex[1L], 2))
  expect_equal(length(unique(manifest$accession)), 3)
  recipes <- lapply(result$review$recipe[1:3], jsonlite::read_json)
  expect_true(all(vapply(recipes, function(x) isTRUE(x$verified), logical(1))))
  expect_identical(recipes[[2L]]$mapping$b, "beta_CSA")
  expect_identical(recipes[[2L]]$p_scale, "neglog10p")
  expect_identical(recipes[[1L]]$mapping$eaf, "af_EUR")
  expect_null(recipes[[3L]]$mapping$eaf)
  expect_identical(recipes[[1L]]$build_check$build, "GRCh37")
  expect_identical(recipes[[1L]]$liftover$index_md5, unname(index_before))
  expect_length(list.files(ga$dirs[["gstat"]]), 0L)
  expect_length(list.files(ga$dirs[["download"]]), 0L)
  count <- length(requests)
  again <- scan(selection = data.frame(filename = filenames[c(1, 1, 2)], analysis = c("EUR", "CSA", "EUR")))
  expect_true(all(again$review$reused))
  expect_length(requests, count)
  short <- scanGWAS(GAlist = ga, repository = "pan-ukb", files = "original",
    phenotype_manifest = metadata, selection = data.frame(filename = filenames[1L], analysis = "EUR"),
    recipe_dir = file.path(root, "short-recipes"), preview_rows = 3L)
  expect_identical(short$review$status, "needs_review")
  expect_match(short$review$issues, "does not verify")
  recipe_path <- result$review$recipe[1L]
  previous <- readBin(recipe_path, "raw", n = file.size(recipe_path))
  altered <- recipes[[1L]]
  altered$build_check$index_md5 <- strrep("0", 32)
  jsonlite::write_json(altered, recipe_path, auto_unbox = TRUE)
  stale <- scan(selection = data.frame(filename = filenames[1L], analysis = "EUR"))
  expect_identical(stale$review$status, "needs_review")
  expect_match(stale$review$issues, "Saved recipe build check differs")
  writeBin(previous, recipe_path)
  local_mocked_bindings(.gact_cpp_gwas_preview = function(...) stop("Unexpected new preview"),
                        .package = "gact")
  expect_identical(runGWASManifest(result$manifest, descriptor, engine = "cpp",
                                  import = FALSE)$status, rep("ready", 3))
  downloads <- character()
  local_mocked_bindings(curl_download = function(url, destfile, ...) {
    downloads <<- c(downloads, url)
    stopifnot(file.copy(upstream[match(basename(url), filenames)], destfile))
  }, .package = "curl")
  imported <- runGWASManifest(result$manifest, descriptor, engine = "cpp")
  expect_identical(imported$status, rep("imported", 3))
  expect_length(downloads, 2L)
  expect_equal(imported$n_final, rep(12, 3))
  saved <- readRDS(descriptor)
  for (i in seq_len(3)) {
    out <- getGWAS(saved, studyID = imported$study_id[i], rsids = catalog$rsids)
    expect_equal(out$b, rep(c(.2, .3, .2)[i], 12))
    expect_equal(out$pos, 101:112)
    expect_true(all(is.na(out$n)))
  }
  expect_identical(tools::md5sum(file.path(map, "index.bin")), index_before)
  expect_error(scan(selection = data.frame(filename = filenames[1L], analysis = "hq")),
               "single-ancestry")
  bad <- transform(manifest, md5 = c(strrep("0", 32), md5[2:3]))
  expect_match(runGWASManifest(bad, descriptor, engine = "cpp")$issues[1L],
               "duplicate filename|different input checksum")
})
