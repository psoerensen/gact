test_that("recipe TSV includes all verified JSONs without modifying them", {
  root <- tempfile("recipe-inventory-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  save <- function(name, value) jsonlite::write_json(value, file.path(root, name),
    auto_unbox = TRUE, na = "null", null = "null", digits = NA)
  save("z.json", list(verified = TRUE, accession = "GCST90000002", filename = "b.tsv.gz",
    input_md5 = strrep("b", 32), input_url = "https://example.org/b.tsv.gz",
    trait = "Metabolite\tconcentration", study_n_total = 338919,
    study_n_case = NULL, type = "quantitative", ancestry = "EUR", build = "GRCh38"))
  save("a.json", list(verified = TRUE, accession = "PANUKB-EXAMPLE-AFR", filename = "a.tsv.bgz",
    input_md5 = strrep("a", 32), study_n_total = 7015, study_n_case = 100,
    study_n_control = 6915, trait = "Diabetes", type = "binary", ancestry = "AFR"))
  save("draft.json", list(verified = FALSE, accession = "DRAFT"))
  save("formats.json", list(formats = list(list(verified = TRUE))))
  dir.create(file.path(root, "reviews"))
  jsonlite::write_json(list(verified = TRUE, accession = "DRAFT"),
                       file.path(root, "reviews", "draft.json"))
  paths <- list.files(root, pattern = "json$", full.names = TRUE)
  hashes <- tools::md5sum(paths)
  output <- gact:::.gact_write_recipe_inventory(root)
  expect_identical(output, file.path(root, "recipes.tsv"))
  table <- data.table::fread(output, data.table = FALSE)
  expect_identical(table$accession, c("GCST90000002", "PANUKB-EXAMPLE-AFR"))
  expect_identical(table$recipe_file, c("z.json", "a.json"))
  expect_equal(table$study_n_total, c(338919, 7015))
  expect_equal(table$study_n_case, c(NA, 100))
  expect_equal(table$study_n_control, c(NA, 6915))
  expect_identical(table$trait[1L], "Metabolite\tconcentration")
  expect_identical(tools::md5sum(paths), hashes)
  before <- tools::md5sum(output)
  gact:::.gact_write_recipe_inventory(root)
  expect_identical(tools::md5sum(output), before)
  expect_length(list.files(root, all.files = TRUE, pattern = "^\\.recipes-"), 0L)
})

test_that("empty inventories have stable columns and malformed JSONs are reported", {
  root <- tempfile("empty-inventory-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  output <- gact:::.gact_write_recipe_inventory(root)
  expect_equal(nrow(data.table::fread(output)), 0L)
  fields <- names(data.table::fread(output))
  expect_true(all(c("accession", "study_n_total", "input_md5", "recipe_file") %in% fields))
  writeLines("{broken", file.path(root, "broken.json"))
  jsonlite::write_json(list(verified = TRUE, accession = "GCST90000003",
    input_md5 = strrep("c", 32), trait = c("a", "b")), file.path(root, "nonscalar.json"),
    auto_unbox = TRUE)
  expect_warning(gact:::.gact_write_recipe_inventory(root), "Skipped 2.*broken.json.*nonscalar.json")
  expect_identical(names(data.table::fread(output)), fields)
  expect_equal(nrow(data.table::fread(output)), 0L)
})

test_that("failed inventory publication preserves the previous TSV", {
  root <- tempfile("atomic-inventory-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  output <- gact:::.gact_write_recipe_inventory(root)
  before <- tools::md5sum(output)
  local_mocked_bindings(.publish_gwas_import_file = function(...) stop("Simulated publication failure"),
                        .package = "gact")
  expect_error(gact:::.gact_write_recipe_inventory(root), "Simulated publication failure")
  expect_identical(tools::md5sum(output), before)
  expect_length(list.files(root, all.files = TRUE, pattern = "^\\.recipes-"), 0L)
})

test_that("Catalog recipes retain explicit file-level case/control counts without inventing a split", {
  id <- "GCST90000001"
  hash <- strrep("a", 32)
  sample <- list(sample_ancestry_category = "European", sample_size = 1000,
                 case_control_study = TRUE, case_count = 100, control_count = 900)
  local_mocked_bindings(.gact_catalog_file = function(...) {
    file <- paste0(id, ".h.tsv.gz")
    list(file = file, url = paste0("https://example.org/", file), base = "https://example.org/",
      fetch = function(...) paste(hash, file),
      meta = list(gwas_id = id, data_file_name = file, data_file_md5sum = hash,
        genome_assembly = "GRCh38", coordinate_system = "1-based", is_harmonised = TRUE,
        trait_description = "Disease", sex = "combined", samples = list(sample),
        author_notes = "Genome-wide association summary statistics"))
  }, .package = "gact")
  row <- gact:::.gact_catalog_recipe_row(id, NULL)
  expect_equal(c(row$study_n_total, row$study_n_case, row$study_n_control), c(1000, 100, 900))
  sample$control_count <- 899
  expect_error(gact:::.gact_catalog_recipe_row(id, NULL), "counts disagree")
  sample$case_count <- -100
  expect_error(gact:::.gact_catalog_recipe_row(id, NULL), "Invalid Catalog case_count")
  sample$case_count <- sample$control_count <- NULL
  row <- gact:::.gact_catalog_recipe_row(id, NULL)
  expect_true(is.na(row$study_n_case))
  expect_true(is.na(row$study_n_control))
  expect_equal(row$study_n_total, 1000)
})
