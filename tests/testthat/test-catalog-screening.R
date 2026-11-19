test_that("Catalog screening resolves missing coverage wording using provider file evidence", {
  meta <- list(genotyping_technology = "Genome-wide genotyping array")
  asset <- list(meta = meta, base = "https://example.org/harmonised/",
    fetch = function(...) "93.90% ( 8554961 of 9110781 ) sites successfully harmonised.")
  check <- gact:::.gact_catalog_coverage
  result <- check(asset, NULL, "GCST90432071")
  expect_identical(result$coverage, "genome_wide")
  expect_match(result$evidence, "8554961")
  expect_match(result$evidence, "not a full-file coverage audit")
  asset$fetch <- function(...) "100% ( 1200 of 1200 ) sites successfully harmonised."
  expect_error(check(asset, NULL, "GCST90432071"), class = "gact_coverage_review")
  asset$fetch <- function(...) stop("Unavailable log")
  expect_error(check(asset, NULL, "GCST90432071"), class = "gact_coverage_review")
  asset$meta$author_notes <- "cis-only genome-wide association results"
  asset$fetch <- function(...) stop("Restricted files must not fetch a log")
  expect_error(check(asset, "Reviewed file", "GCST90432071"), "restricted|Metadata indicates")
  asset$meta <- list()
  asset$fetch <- function(...) "100% ( 8000000 of 8000000 ) sites successfully harmonised."
  expect_error(check(asset, NULL, "GCST90432071"), class = "gact_coverage_review")
})

test_that("single ancestry is independent of literal cohort wording and cannot override meta-analysis", {
  check <- gact:::.gact_catalog_study_scope
  meta <- list(samples = list(list(sample_ancestry_category = "European")))
  expect_identical(check(meta, "GCST91007350")$ancestry, "EUR")
  for (software in c("METAL", "MTAG v0.9.0", "MR-MEGA", "METASOFT", "GWAMA")) {
    meta$analysis_software <- software
    expect_error(check(meta, "GCST91007350"), class = "gact_scan_excluded")
  }
  meta$analysis_software <- "REGENIE"
  meta$author_notes <- "European meta-analysis"
  expect_error(check(meta, "GCST91007350"), class = "gact_scan_excluded")
  meta$author_notes <- NULL
  meta$samples[[1L]]$sample_ancestry_category <- c("European", "East Asian")
  expect_error(check(meta, "GCST91007350"), class = "gact_scan_excluded")
  meta$samples[[1L]]$sample_ancestry_category <- "Unknown"
  expect_error(check(meta, "GCST91007350"), "single ancestry")
  for (ancestry in c("Mixed", "European, East Asian", "European and East Asian")) {
    meta$samples[[1L]]$sample_ancestry_category <- ancestry
    expect_error(check(meta, "GCST91007350"), class = "gact_scan_excluded")
  }
  meta$samples[[1L]]$sample_ancestry_category <- "African American or Afro-Caribbean"
  expect_identical(check(meta, "GCST91007350")$ancestry, "AFR")
  meta$samples[[1L]]$sample_ancestry_category <- "African American or Afro-Caribbean"
  expect_identical(check(meta, "GCST90274421")$ancestry, "AFR")
})

test_that("Catalog metadata reads UTF-8 author notes on Windows as well as Linux", {
  root <- tempfile("catalog-utf8-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  path <- file.path(root, "GCST90432071.h.tsv.gz-meta.yaml")
  writeLines(enc2utf8(c("gwas_id: GCST90432071", "author_notes: \"Each individual\u2019s ancestry\"")),
             path, useBytes = TRUE)
  meta <- gact:::.gact_catalog_file("GCST90432071", metadata_dir = root)$meta
  expect_identical(meta$gwas_id, "GCST90432071")
  expect_identical(meta$author_notes, "Each individual\u2019s ancestry")
})

test_that("standard and legacy Catalog previews preserve harmonised coordinates and effect scale", {
  sample <- data.frame(chromosome = 1, base_pair_location = 727717,
    effect_allele = "C", other_allele = "G", beta = -.02375591,
    standard_error = .05768172, effect_allele_frequency = .7116421,
    p_value = .6804532, rsid = "rs61769340", hazard_ratio = .976524)
  profile <- gact:::.gact_catalog_preview_profile(sample, "Non-GWAS-SSF")
  recipe <- gact:::.gact_gwas_scan_recipe(list(md5 = strrep("a", 32), url = "example"),
                                        sample, 100L, profile)
  expect_length(recipe$inspection_issues, 0L)
  expect_true(recipe$validation$ok)
  expect_identical(recipe$mapping$b, "beta")
  expect_identical(recipe$effect_scale, "beta")
  sample$hazard_ratio <- 2
  profile <- gact:::.gact_catalog_preview_profile(sample, "Non-GWAS-SSF")
  expect_match(profile$issues, "disagrees")
  sample$hazard_ratio <- NULL
  for (field in c("chr", "pos", "ea", "nea", "b", "eaf", "marker")) {
    old <- c(chr = "chromosome", pos = "base_pair_location", ea = "effect_allele",
             nea = "other_allele", b = "beta", eaf = "effect_allele_frequency", marker = "rsid")
    new <- c(chr = "hm_chrom", pos = "hm_pos", ea = "hm_effect_allele",
             nea = "hm_other_allele", b = "hm_beta", eaf = "hm_effect_allele_frequency", marker = "hm_rsid")
    sample[[new[field]]] <- sample[[old[field]]]
  }
  sample$base_pair_location <- 123 # Source coordinates must not override hm_*.
  profile <- gact:::.gact_catalog_preview_profile(sample, "legacy")
  expect_identical(profile$mapping$pos, "hm_pos")
  expect_identical(profile$mapping$b, "hm_beta")
  sample$hm_beta <- NA_real_
  profile <- gact:::.gact_catalog_preview_profile(sample, "legacy")
  expect_match(profile$issues, "No usable harmonised")
})

test_that("report-derived Catalog cases save recipes and recheck exclusions on resume", {
  root <- tempfile("catalog-screening-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  ids <- c("GCST90432071", "GCST91007350", "GCST000028", "GCST90270929", "GCST90270930")
  calls <- 0L
  software <- "REGENIE"
  scan_ids <- ids
  local_mocked_bindings(
    .gact_catalog_file = function(id, catalog_path = NULL) {
      i <- match(id, ids)
      file <- paste0(id, ".h.tsv.gz")
      hash <- strrep(as.character(i), 32)
      list(file = file, url = paste0("https://example.org/", file), base = "https://example.org/",
        meta = list(gwas_id = id, data_file_name = file, genome_assembly = "GRCh38",
          coordinate_system = "1-based", is_harmonised = TRUE,
          file_type = c("Non-GWAS-SSF", "GWAS-SSF v1.0", "", "GWAS-SSF v1.0", "GWAS-SSF v1.0")[i],
          analysis_software = if (i == 2L) "MTAG" else software,
          genotyping_technology = "Genome-wide genotyping array", data_file_md5sum = hash,
          trait_description = "Metabolite level", sex = "combined",
          samples = list(list(case_control_study = FALSE, sample_size = 6840,
            sample_ancestry_category = if (i == 5L) c("European", "East Asian") else "European"))),
        fetch = function(remote, local) if (endsWith(remote, ".running.log"))
          if (i == 4L) "Log unavailable" else "99% ( 8554961 of 9110781 ) sites successfully harmonised." else
          paste(hash, file))
    },
    .gact_fetch_gwas_prefix = function(url, path, max_bytes) {
      calls <<- calls + 1L
      con <- gzfile(path, "wt")
      on.exit(close(con))
      writeLines(c("chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\tp_value",
                   "1\t727717\tC\tG\t0.1\t0.2\t0.4"), con)
      100L
    }, .package = "gact")
  scan <- function() scanGWAS(scan_ids, out_dir = tempfile(tmpdir = root),
    recipe_dir = file.path(root, "recipes"), repository = "gwas-catalog", preview_rows = 1L)
  result <- scan()
  expect_true(all(is.na(result$review$issues[c(1L, 3L)])),
              info = paste(result$review$issues[c(1L, 3L)], collapse = " | "))
  expect_identical(result$review$status, c("ready", "excluded", "ready", "needs_review", "excluded"))
  expect_identical(calls, 2L)
  expect_identical(data.table::fread(result$manifest)$accession, ids[c(1L, 3L)])
  expect_length(list.files(file.path(root, "recipes"), pattern = "json$"), 2L)
  expect_identical(result$recipes_tsv,
    normalizePath(file.path(root, "recipes", "recipes.tsv"), winslash = "/"))
  expect_setequal(data.table::fread(result$recipes_tsv)$accession, ids[c(1L, 3L)])
  paths <- list.files(file.path(root, "recipes"), pattern = "json$", full.names = TRUE)
  hashes <- tools::md5sum(paths)
  resumed <- scan()
  expect_true(all(resumed$review$reused[c(1L, 3L)]))
  expect_identical(calls, 2L)
  scan_ids <- ids[1L]
  partial <- scan()
  expect_true(partial$review$reused)
  expect_identical(data.table::fread(partial$manifest)$accession, ids[1L])
  expect_setequal(data.table::fread(partial$recipes_tsv)$accession, ids[c(1L, 3L)])
  expect_identical(tools::md5sum(paths), hashes)
  scan_ids <- ids
  software <- "METAL"
  changed <- scan()
  expect_true(all(changed$review$status == "excluded"))
  expect_false(any(changed$review$reused))
  expect_identical(calls, 2L)
})
