test_that("coverage checks require explicit scope and keep molecular GWAS", {
  check <- gact:::.gact_gwas_coverage
  for (trait in c("Protein abundance", "Metabolite level", "Type 2 diabetes")) {
    meta <- list(trait_description = trait, author_notes = "Genome-wide association summary statistics",
                 genotyping_technology = "Genome-wide genotyping array")
    expect_identical(check(meta)$coverage, "genome_wide")
  }
  restrictions <- c("Genome-wide association study; cis-only results",
    "Targeted-region SNP results", "Candidate gene associations",
    "Exome-only association summary statistics", "Genome-wide gene-based tests",
    "Significant variants only", "Non-genome-wide results", "Not a genome-wide association file",
    "cis-pQTL summary statistics")
  for (notes in restrictions) {
    expect_identical(check(list(author_notes = notes), "Reviewed paper")$coverage, "restricted")
  }
  expect_identical(check(list(genotyping_technology = "Whole-exome sequencing"))$coverage, "restricted")
  expect_identical(check(list(genotyping_technology = "Exome genotyping array"))$coverage, "restricted")
  expect_identical(check(list(genotyping_technology = "Genome-wide genotyping array"))$coverage, "unclear")
  expect_identical(check(list(author_notes = "Full summary statistics"))$coverage, "unclear")
  expect_identical(check(list(author_notes = "cis-pQTLs discussed in the paper"))$coverage, "unclear")
  expect_identical(check(list(), "Reviewed complete genome-wide file in provider README")$coverage, "genome_wide")
})

test_that("Catalog coverage is checked before previews and before recipe reuse", {
  root <- tempfile("coverage-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  ga <- createDB("hsa.GRCh38.0.0.1", root)
  ids <- paste0("GCST", 91000001:91000005)
  scopes <- c("single-cohort genome-wide association summary statistics",
              "single-cohort genome-wide association summary statistics",
              "single-cohort cis-only results", "single-cohort study",
              "single-cohort whole-exome association analysis")
  traits <- c("Protein abundance", "Metabolite concentration", "Protein cis-pQTL",
              "Metabolite concentration", "Disease")
  previews <- character()
  profile <- gact:::.gact_panukb_catalog_profile("Catalog")
  local_mocked_bindings(
    .gact_catalog_file = function(id, catalog_path = NULL) {
      i <- match(id, ids)
      file <- paste0(id, ".h.tsv.gz")
      hash <- strrep(as.character(i), 32)
      list(file = file, url = paste0("https://example.org/", file),
        base = "https://example.org/",
        meta = list(gwas_id = id, data_file_name = file, genome_assembly = "GRCh38",
          data_file_md5sum = hash,
          coordinate_system = "1-based", is_harmonised = TRUE, file_type = "GWAS-SSF v1.0",
          author_notes = scopes[i], trait_description = traits[i], sex = "combined",
          genotyping_technology = "Genome-wide genotyping array",
          samples = list(list(case_control_study = FALSE, sample_size = 1000,
                              sample_ancestry_category = list("European")))),
        fetch = function(...) paste(hash, file))
    },
    .gact_fetch_gwas_prefix = function(url, path, max_bytes) {
      previews <<- c(previews, url)
      con <- gzfile(path, "wt")
      on.exit(close(con))
      row <- c("1", "10", "A", "C", ".2", ".1", ".2", ".05", "rs1", "lo", "10", "1_10_C_A")
      writeLines(c(paste(profile$headers, collapse = "\t"), paste(row, collapse = "\t")), con)
      100L
    }, .package = "gact")
  scan <- function(...) scanGWAS(ids, GAlist = ga, repository = "gwas-catalog",
                                 preview_rows = 1L, ...)
  result <- scan()
  expect_identical(result$review$issues[1:2], rep(NA_character_, 2L))
  expect_identical(result$review$status, c("ready", "ready", "excluded", "needs_review", "excluded"))
  expect_identical(result$review$coverage, c("genome_wide", "genome_wide", "restricted", "unclear", "restricted"))
  expect_length(previews, 2L)
  manifest <- data.table::fread(result$manifest)
  expect_identical(manifest$trait, traits[1:2])
  recipe <- jsonlite::read_json(result$review$recipe[1L])
  expect_identical(recipe$coverage, "genome_wide")
  expect_match(recipe$coverage_evidence, "Provider metadata")
  expect_length(list.files(ga$dirs[["gstat"]]), 0L)
  expect_length(list.files(ga$dirs[["download"]]), 0L)
  expect_true(all(scan()$review$reused[1:2]))
  expect_length(previews, 2L)
  scopes[1L] <- "single-cohort study"
  changed <- scan()
  expect_identical(changed$review$status[1L], "needs_review")
  expect_false(changed$review$reused[1L])
  expect_length(previews, 2L)
  evidence <- stats::setNames("Reviewed genome-wide results from the provider README", ids[4L])
  reviewed <- scan(coverage_evidence = evidence)
  expect_identical(reviewed$review$status[4L], "ready")
  expect_length(previews, 3L)
  expect_error(scan(coverage_evidence = "unnamed evidence"), "named by Catalog accession")
})
