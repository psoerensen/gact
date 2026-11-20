ingest_recipe_fixture <- function() {
  root <- tempfile("ingest-recipes-")
  dir.create(root)
  root <- normalizePath(root, winslash = "/")
  dirs <- setNames(file.path(root, c("marker", "gstat", "download", "recipes")),
                    c("marker", "gstat", "download", "recipes"))
  for (path in dirs) dir.create(path)
  markers <- data.frame(rsids = paste0("rs", 1:3), chr = 1, pos = 101:103,
                         ea = c("A", "C", "A"), nea = c("C", "T", "T"), eaf = c(.2, .3, .4))
  data.table::fwrite(markers, file.path(dirs[["marker"]], "markers.txt.gz"))
  data.table::fwrite(data.frame(cpra = with(markers, paste(chr, pos, ea, nea, sep = "_"))),
    file.path(dirs[["marker"]], "cpra.txt.gz"), col.names = FALSE)
  ga <- list(dirs = dirs, rsids = markers$rsids, study = list(id = character()),
             metadata = list(build = "GRCh38"))
  descriptor <- file.path(root, "GAlist.rds")
  saveRDS(ga, descriptor)
  paths <- files <- character()
  ids <- paste0("GCST", 90000001:90000002)
  for (i in 1:2) {
    source <- data.frame(chromosome = 1, base_pair_location = 101:103, effect_allele = c("A", "T", "A"),
      other_allele = c("C", "C", "T"), beta = c(.5, .2, .3) * i, standard_error = .1,
      effect_allele_frequency = c(.2, .7, .4), p_value = .05, n = c(NA, 800, 950))
    filename <- paste0(ids[i], ".h.tsv.gz")
    file <- file.path(dirs[["download"]], filename)
    data.table::fwrite(source, file, sep = "\t")
    row <- data.frame(accession = ids[i], filename = filename, url = "",
      md5 = unname(tools::md5sum(file)), source = "Synthetic Catalog", trait = paste("Trait", i),
      type = c("quantitative", "binary")[i], build = "GRCh38", ancestry = "EUR",
      gender = "both", category = "Test", efo = "EFO_TEST", strand = "forward",
      reference = gact:::.gact_gwas_catalog_url(ids[i]), study_n_total = 1000,
      study_n_case = if (i == 2) 100 else NA_real_, study_n_control = if (i == 2) 900 else NA_real_)
    recipe <- gact:::.gact_gwas_scan_recipe(row, source, 100L,
      gact:::.gact_catalog_preview_profile(source, row$source))
    recipe$verified <- TRUE
    path <- gact:::.gact_verified_recipe_path(dirs[["recipes"]], ids[i], row$md5)
    gact:::.gact_save_verified_recipe(recipe, path, row)
    paths <- c(paths, path)
    files <- c(files, file)
  }
  tsv <- gact:::.gact_write_recipe_inventory(dirs[["recipes"]])
  list(root = root, ga = ga, descriptor = descriptor, ids = ids, paths = paths, files = files, tsv = tsv)
}

test_that("ingestGWAS requires explicit selection and uses JSON metadata", {
  x <- ingest_recipe_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  before <- tools::md5sum(c(x$descriptor, x$paths, x$tsv))
  expect_error(ingestGWAS(x$ga), "explicitly")
  expect_error(ingestGWAS(x$ga, character()), "explicitly")
  expect_error(ingestGWAS(x$ga, "absent"), "exactly one")
  expect_error(ingestGWAS(x$ga, rep(x$ids[1L], 2)), "unique")
  selected <- data.table::fread(x$tsv, data.table = FALSE)[1L, ]
  selected$trait <- "TSV edits must not override JSON"
  local_mocked_bindings(.gact_gwas_manifest_file = function(...) stop("Unexpected download"), .package = "gact")
  report <- ingestGWAS(x$ga, selected, import = FALSE)
  expect_identical(report$status, "ready")
  from_table <- ingestGWAS(x$ga, data.table::as.data.table(selected), import = FALSE)
  expect_identical(from_table$status, "ready")
  manifest <- data.table::fread(file.path(attr(report, "out_dir"), "manifest.tsv"))
  expect_identical(manifest$trait, "Trait 1")
  expect_identical(tools::md5sum(c(x$descriptor, x$paths, x$tsv)), before)
  expect_length(list.files(x$ga$dirs[["gstat"]]), 0L)
})

test_that("selected imports assign local IDs, publish metadata and skip CSV history", {
  x <- ingest_recipe_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  before <- tools::md5sum(c(x$paths, x$tsv))
  first <- ingestGWAS(x$ga, data.table::fread(x$tsv)[2L])
  expect_identical(first$ingest_report$status, "imported")
  expect_identical(first$studies$id, "GWAS1")
  ga <- ingestGWAS(x$ga, basename(x$paths[1L])) # Stale caller list is reloaded.
  second <- ga$ingest_report
  expect_identical(second$status, "imported", info = second$issues)
  expect_identical(second$study_id, "GWAS2")
  saved <- readRDS(x$descriptor)
  expect_identical(ga[names(saved)], saved)
  expect_null(saved$ingest_report)
  csv_path <- file.path(x$ga$dirs[["gstat"]], "GWAS_information.csv")
  csv <- read.csv2(csv_path)
  expect_identical(csv$id, c("GWAS1", "GWAS2"))
  expect_identical(csv$accession, rev(x$ids))
  expect_identical(csv$trait, c("Trait 2", "Trait 1"))
  expect_equal(csv$n, c(1000, 1000))
  expect_equal(csv$ncase, c(100, NA))
  expect_equal(csv$ncontrol, c(900, NA))
  expect_equal(csv$neff, c(90, 1000))
  expect_identical(csv$category, rep("Test", 2))
  expect_identical(csv$efo, rep("EFO_TEST", 2))
  expect_equal(csv$recipe_checksum, unname(tools::md5sum(rev(x$paths))))
  for (id in csv$id) {
    stat <- getGWAS(ga, id)
    expect_equal(stat$n, c(NA, 800, 950))
    expect_equal(nrow(stat), 3)
  }
  matrix_dir <- file.path(x$ga$dirs[["gstat"]], "matrix-pilot")
  matrix <- writeMarkerStat(ga, out_dir = matrix_dir, progress = FALSE)
  expect_equal(matrix$dim, c(3, 2))
  z <- readBin(matrix$z, "double", n = 6L, size = 8L, endian = "little")
  expect_equal(z, c(10, -4, 6, 5, -2, 3))
  expect_identical(data.table::fread(matrix$studies)$study_id, c("GWAS1", "GWAS2"))
  unlink(x$files)
  after <- tools::md5sum(c(x$descriptor, csv_path))
  ga <- ingestGWAS(x$descriptor, x$ids)
  repeated <- ga$ingest_report
  expect_identical(ga$studies$id, csv$id)
  expect_identical(repeated$status, rep("already_imported", 2))
  expect_identical(repeated$study_id, c("GWAS2", "GWAS1"))
  expect_identical(tools::md5sum(c(x$descriptor, csv_path)), after)
  expect_identical(tools::md5sum(c(x$paths, x$tsv)), before)
})

test_that("changed sources, changed recipes and inconsistent CSVs require review", {
  x <- ingest_recipe_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  expect_identical(ingestGWAS(x$ga, x$ids[1L])$ingest_report$status, "imported")
  recipe <- jsonlite::read_json(x$paths[1L])
  recipe$input_md5 <- recipe$md5 <- strrep("f", 32)
  path <- gact:::.gact_verified_recipe_path(dirname(x$paths[1L]), recipe$accession, recipe$input_md5)
  jsonlite::write_json(recipe, path, auto_unbox = TRUE, na = "null", digits = NA)
  gact:::.gact_write_recipe_inventory(dirname(path))
  expect_error(ingestGWAS(x$ga, x$ids[1L]), "exactly one")
  before <- tools::md5sum(x$descriptor)
  changed <- ingestGWAS(x$ga, basename(path))$ingest_report
  expect_identical(changed$status, "needs_review")
  expect_match(changed$issues, "different input checksum")
  recipe <- jsonlite::read_json(x$paths[1L])
  recipe$comments <- "Changed after import"
  jsonlite::write_json(recipe, x$paths[1L], auto_unbox = TRUE, na = "null", digits = NA)
  altered <- ingestGWAS(x$ga, basename(x$paths[1L]))$ingest_report
  expect_identical(altered$status, "needs_review")
  expect_match(altered$issues, "different verified recipe")
  expect_identical(tools::md5sum(x$descriptor), before)
  csv_path <- file.path(x$ga$dirs[["gstat"]], "GWAS_information.csv")
  csv <- read.csv2(csv_path)
  csv$input_md5 <- strrep("0", 32)
  write.csv2(csv, csv_path, row.names = FALSE)
  expect_error(ingestGWAS(x$ga, basename(x$paths[1L])), "import identity differs")
})

test_that("a failed selected import publishes no history and can be retried", {
  x <- ingest_recipe_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  valid <- data.table::fread(x$files[1L], data.table = FALSE)
  invalid <- valid
  invalid$n[3L] <- Inf
  data.table::fwrite(invalid, x$files[1L], sep = "\t")
  recipe <- jsonlite::read_json(x$paths[1L])
  recipe$input_md5 <- recipe$md5 <- unname(tools::md5sum(x$files[1L]))
  path <- gact:::.gact_verified_recipe_path(dirname(x$paths[1L]), recipe$accession, recipe$input_md5)
  jsonlite::write_json(recipe, path, auto_unbox = TRUE, na = "null", digits = NA)
  gact:::.gact_write_recipe_inventory(dirname(path))
  before <- tools::md5sum(x$descriptor)
  ga <- ingestGWAS(x$ga, basename(path))
  expect_identical(ga$study$id, character())
  failed <- ga$ingest_report
  expect_identical(failed$status, "needs_review")
  expect_match(failed$issues, "Per-variant n")
  expect_identical(tools::md5sum(x$descriptor), before)
  expect_length(list.files(x$ga$dirs[["gstat"]], all.files = TRUE, no.. = TRUE), 0L)
  data.table::fwrite(valid, x$files[1L], sep = "\t")
  # Gzip timestamps can change; key the corrected recipe to its exact full input.
  recipe$input_md5 <- recipe$md5 <- unname(tools::md5sum(x$files[1L]))
  corrected <- gact:::.gact_verified_recipe_path(dirname(path), recipe$accession, recipe$input_md5)
  jsonlite::write_json(recipe, corrected, auto_unbox = TRUE, na = "null", digits = NA)
  gact:::.gact_write_recipe_inventory(dirname(path))
  ga <- ingestGWAS(ga, basename(corrected))
  retried <- ga$ingest_report
  expect_identical(ga$studies$id, "GWAS1")
  expect_identical(retried$status, "imported", info = retried$issues)
  expect_identical(retried$study_id, "GWAS1")
})

test_that("duplicate identities distinguish ancestries and flag incomplete legacy history", {
  study <- list(id = "GWAS1", accession = "GCST90000001", ancestry = "EUR", input_md5 = strrep("a", 32))
  row <- data.frame(accession = "GCST90000001", ancestry = "EAS")
  check <- gact:::.gact_gwas_inventory_match
  expect_length(check(study, row, strict = TRUE), 0L)
  row$ancestry <- "EUR"
  expect_identical(check(study, row, strict = TRUE), 1L)
  study$ancestry <- NA_character_
  expect_error(check(study, row, strict = TRUE), "Stored ancestry is missing")
  study$ancestry <- "EUR"
  study$input_md5 <- NA_character_
  expect_error(check(study, row, strict = TRUE), "Stored input checksum is missing")
})

test_that("full and Z-only imports append aligned columns and roll back together", {
  x <- ingest_recipe_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  gstat <- x$ga$dirs[["gstat"]]
  dir <- file.path(gstat, "matrix")
  z <- file.path(dir, "z.f64.bin")
  studies <- file.path(dir, "studies.tsv")
  first <- ingestGWAS(x$ga, x$ids[1])
  expect_identical(first$ingest_report$status, "imported", info = first$ingest_report$issues)
  expect_equal(readBin(z, "double", 3, endian = "little"), c(5, -2, 3))
  expect_true(file.exists(first$studyfiles["GWAS1"]))
  files <- c(z, studies, x$descriptor, file.path(gstat, "GWAS_information.csv"))
  before <- tools::md5sum(files)
  publish <- .publish_gwas_import_file
  fail <- TRUE
  local_mocked_bindings(.publish_gwas_import_file = function(from, to, ...) {
    if (fail && identical(to, x$descriptor)) stop("test descriptor failure")
    publish(from, to, ...)
  }, .package = "gact")
  failed <- ingestGWAS(first, x$ids[2], storage = "z")
  expect_match(failed$ingest_report$issues, "test descriptor failure")
  expect_identical(tools::md5sum(files), before)
  fail <- FALSE
  ga <- ingestGWAS(first, x$ids[2], storage = "z")
  expect_identical(ga$ingest_report$status, "imported", info = ga$ingest_report$issues)
  expect_identical(ga$studies$storage, c("full", "z"))
  expect_equal(readBin(z, "double", 6, endian = "little"), c(5, -2, 3, 10, -4, 6))
  expect_equal(data.table::fread(studies)$observed, c(3, 3))
  expect_false(file.exists(file.path(gstat, "GWAS2.txt.gz")))
  expect_true(file.exists(file.path(gstat, "GWAS2.qc.json")))
  expect_error(getGWAS(ga, "GWAS2"), "Only Z scores")
  expect_equal(as.numeric(getGWAS(ga, "GWAS2", what = "z")), c(10, -4, 6))
  subset <- ga
  subset$rsids <- c("rs3", "rs2")
  subset$studies <- subset$studies[2:1, ]
  named <- writeMarkerStat(subset, file.path(gstat, "pilot"), progress = FALSE)
  expect_equal(readBin(named$z, "double", 4, endian = "little"), c(6, -4, 3, -2))
  expect_identical(ingestGWAS(ga, x$ids)$ingest_report$status, rep("already_imported", 2))
  expect_equal(file.size(z), 48)
})

test_that("one worker can finish alignment while another is still downloading", {
  skip_on_os("windows")
  x <- ingest_recipe_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  scratch <- file.path(x$root, "scratch")
  dir.create(scratch)
  upstream <- file.path(x$root, basename(x$files))
  expect_true(all(file.rename(x$files, upstream)))
  for (i in seq_along(upstream)) {
    recipe <- jsonlite::read_json(x$paths[i])
    recipe$input_url <- paste0("https://example.org/", basename(upstream[i]))
    jsonlite::write_json(recipe, x$paths[i], auto_unbox = TRUE, digits = NA)
  }
  local_mocked_bindings(curl_download = function(url, destfile, mode) {
    i <- match(basename(url), basename(upstream))
    writeLines(as.character(Sys.getpid()), file.path(x$root, paste0("pid", i)))
    if (i == 2L) {
      until <- Sys.time() + 8
      while (!length(list.files(scratch, "^study_00001.qc.json$", recursive = TRUE))) {
        if (Sys.time() >= until) stop("First worker never staged its study during download")
        Sys.sleep(.02)
      }
    }
    if (!file.copy(upstream[i], destfile)) stop("Cannot copy fixture")
  }, .package = "curl")
  ga <- ingestGWAS(x$ga, x$ids, workers = 2, storage = "z", scratch_dir = scratch)
  expect_identical(ga$ingest_report$status, rep("imported", 2), info = paste(ga$ingest_report$issues, collapse = "; "))
  expect_length(unique(vapply(file.path(x$root, c("pid1", "pid2")), readLines, character(1))), 2)
  expect_identical(ga$ingest_report$cleanup, rep("removed", 2))
  expect_false(any(file.exists(x$files)))
})

test_that("absent and unusable effects occupy missing backbone rows in both modes", {
  for (mode in c("full", "z")) {
    x <- ingest_recipe_fixture()
    on.exit(unlink(x$root, recursive = TRUE), add = TRUE)
    input <- data.table::fread(x$files[1L])[2:1]
    input$beta[2L] <- NA_real_
    data.table::fwrite(input, x$files[1L], sep = "\t")
    recipe <- jsonlite::read_json(x$paths[1L])
    recipe$input_md5 <- recipe$md5 <- unname(tools::md5sum(x$files[1L]))
    path <- .gact_verified_recipe_path(dirname(x$paths[1L]), recipe$accession, recipe$input_md5)
    jsonlite::write_json(recipe, path, auto_unbox = TRUE, digits = NA)
    .gact_write_recipe_inventory(dirname(path))
    ga <- ingestGWAS(x$ga, basename(path), storage = mode)
    expect_identical(ga$ingest_report$status, "imported", info = ga$ingest_report$issues)
    dir <- file.path(ga$dirs[["gstat"]], "matrix")
    expect_equal(readBin(file.path(dir, "z.f64.bin"), "double", 3, endian = "little"),
                 c(NaN, -2, NaN))
    expect_equal(data.table::fread(file.path(dir, "studies.tsv"))$observed, 1)
    expect_equal(ga$ingest_report$n_final, 2)
    unlink(x$root, recursive = TRUE)
  }
})

test_that("retained rows without usable Z never publish in either storage mode", {
  x <- ingest_recipe_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  input <- data.table::fread(x$files[1L])
  input$beta <- NA_real_
  data.table::fwrite(input, x$files[1L], sep = "\t")
  recipe <- jsonlite::read_json(x$paths[1L])
  recipe$input_md5 <- unname(tools::md5sum(x$files[1L]))
  path <- .gact_verified_recipe_path(dirname(x$paths[1L]), recipe$accession, recipe$input_md5)
  jsonlite::write_json(recipe, path, auto_unbox = TRUE, digits = NA)
  .gact_write_recipe_inventory(dirname(path))
  for (storage in c("full", "z")) {
    result <- ingestGWAS(x$ga, basename(path), storage = storage)$ingest_report
    expect_identical(result$status, "needs_review")
    expect_match(result$issues, "No finite Z")
    expect_equal(c(result$n_input, result$n_final, result$n_observed_z,
                   result$n_unusable_effect, result$n_missing_backbone), c(3, 3, 0, 3, 0))
    qc <- jsonlite::read_json(result$qc_file)
    expect_equal(qc$n_observed_z, 0)
    expect_length(list.files(x$ga$dirs[["gstat"]], all.files = TRUE, no.. = TRUE), 0L)
  }
})

test_that("full and Z-only apply the same provider filters and retain study EAF", {
  observed <- list()
  for (storage in c("full", "z")) {
    x <- ingest_recipe_fixture()
    on.exit(unlink(x$root, recursive = TRUE), add = TRUE)
    input <- data.table::fread(x$files[1L])
    input$effect_allele_frequency[1L] <- NA_real_
    input$low_confidence_EUR <- FALSE
    input$low_confidence_EAS <- TRUE
    input$INFO <- .99
    input <- rbind(input, transform(input[3L], low_confidence_EUR = TRUE),
      transform(input[3L], INFO = .2))
    input$low_confidence_EUR[3L] <- TRUE
    data.table::fwrite(input, x$files[1L], sep = "\t")
    recipe <- jsonlite::read_json(x$paths[1L])
    recipe$input_md5 <- unname(tools::md5sum(x$files[1L]))
    recipe$inspection$headers <- names(input)
    recipe$quality_filters <- list(low_confidence = list(column = "low_confidence_EUR",
      evidence = "Fixture provider flag"), info = list(column = "INFO", min = .8,
      evidence = "Fixture provider quality specification"))
    path <- .gact_verified_recipe_path(dirname(x$paths[1L]), recipe$accession, recipe$input_md5)
    jsonlite::write_json(recipe, path, auto_unbox = TRUE, digits = NA)
    .gact_write_recipe_inventory(dirname(path))
    ga <- ingestGWAS(x$ga, basename(path), storage = storage, warn_coverage_below = .9)
    report <- ga$ingest_report
    expect_identical(report$status, "imported", info = report$issues)
    expect_equal(c(report$n_observed_z, report$n_provider_low_confidence, report$n_provider_info),
                 c(2, 2, 1))
    expect_match(report$warnings, "low_confidence_EAS")
    expect_match(report$warnings, "Usable Z coverage")
    observed[[storage]] <- readBin(file.path(ga$dirs[["gstat"]], "matrix", "z.f64.bin"),
                                   "double", 3L, endian = "little")
    if (storage == "full") {
      output <- data.table::fread(ga$studyfiles["GWAS1"])
      expect_equal(output$eaf, c(NA, .3))
      recipe$quality_filters$info$min <- .95
      jsonlite::write_json(recipe, path, auto_unbox = TRUE, digits = NA)
      retry <- ingestGWAS(ga, basename(path))$ingest_report
      expect_identical(retry$status, "needs_review")
      expect_match(retry$issues, "different verified recipe")
    }
  }
  expect_identical(observed$full, observed$z)
  expect_equal(observed$full, c(5, -2, NaN))
})

test_that("restart restores interrupted publication and honours the live lock", {
  x <- ingest_recipe_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  ga <- ingestGWAS(x$ga, x$ids[1L], storage = "z")
  gstat <- ga$dirs[["gstat"]]
  lock <- file.path(gstat, ".gact-import.lock")
  files <- c(x$descriptor, file.path(gstat, "GWAS_information.csv"),
    file.path(gstat, "matrix", c("studies.tsv", "z.f64.bin")))
  before <- tools::md5sum(files)
  publish <- .publish_gwas_import_file
  recover <- .gact_restore_gwas_transaction
  interrupted <- TRUE
  local_mocked_bindings(
    .publish_gwas_import_file = function(from, to, ...) {
      publish(from, to, ...)
      if (interrupted && identical(to, x$descriptor)) stop("Interrupted after descriptor rename")
    },
    .gact_restore_gwas_transaction = function(stage) {
      if (!interrupted) recover(stage)
    }, .package = "gact")
  result <- ingestGWAS(ga, x$ids[2L])
  expect_match(result$ingest_report$issues, "Interrupted")
  expect_true(file.exists(file.path(lock, "owner.rds"), file.path(lock, "journal.rds")) |> all())
  expect_true(file.exists(file.path(gstat, "GWAS2.txt.gz")))
  interrupted <- FALSE
  guard <- .gact_gwas_guard(gstat)
  tryCatch(expect_error(ingestGWAS(ga, x$ids[1L]), "holds the database lock"),
           finally = .Call(C_gact_gwas_unlock, guard))
  expect_true(dir.exists(lock))
  recovered <- ingestGWAS(x$descriptor, x$ids[1L])
  expect_identical(recovered$ingest_report$status, "already_imported")
  expect_identical(tools::md5sum(files), before)
  expect_false(file.exists(file.path(gstat, "GWAS2.txt.gz")))
  expect_false(file.exists(file.path(gstat, "GWAS2.qc.json")))
  expect_false(dir.exists(lock))
  expect_identical(ingestGWAS(recovered, x$ids[2L])$ingest_report$status, "imported")
})

test_that("Catalog meta-analysis recipes keep ancestry and design across a custom directory handoff", {
  x <- ingest_recipe_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  recipes <- file.path(x$root, "reviewed-recipes")
  ga <- x$ga
  ga$dirs[["recipes"]] <- file.path(x$root, "unused-default")
  saveRDS(ga, x$descriptor)
  hashes <- unname(tools::md5sum(x$files))
  local_mocked_bindings(
    .gact_catalog_file = function(id, catalog_path) {
      i <- match(id, x$ids)
      file <- basename(x$files[i])
      list(file = file, url = paste0("https://example.org/", file), base = "https://example.org/",
           meta = list(gwas_id = id, data_file_name = file, genome_assembly = "GRCh38",
             coordinate_system = "1-based", is_harmonised = TRUE,
             trait_description = paste("Trait", i), author_notes = "Genome-wide association meta-analysis",
             n_cohorts = 2L, meta_analysis_method = "Inverse variance",
             data_file_md5sum = hashes[i], samples = list(
               list(sample_ancestry_category = "European", case_control_study = FALSE),
               list(sample_ancestry_category = if (i == 1L) "European" else "East Asian",
                    case_control_study = FALSE))),
           fetch = function(...) "")
    },
    .gact_fetch_gwas_prefix = function(url, path, ...) {
      stopifnot(file.copy(x$files[match(basename(url), basename(x$files))], path))
      file.size(path)
    }, .package = "gact")
  scan <- scanGWAS(x$ids, GAlist = ga, recipe_dir = recipes,
                   repository = "gwas-catalog", preview_rows = 2L)
  expect_identical(scan$review$status, rep("ready", 2))
  inventory <- data.table::fread(scan$recipes_tsv)
  expect_identical(inventory$study_design, rep("meta_analysis", 2))
  expect_identical(inventory$acn, c("EUR", "MIXED"))
  expect_identical(inventory$n_cohorts, c(2L, 2L))
  expect_match(inventory$ancestry_composition[2L], "East Asian")
  expect_identical(jsonlite::read_json(scan$review$recipe[2L])$acn, "MIXED")
  result <- ingestGWAS(ga, x$ids, recipe_dir = scan$recipe_dir, storage = "z")
  expect_identical(result$ingest_report$status, rep("imported", 2),
                   info = paste(result$ingest_report$issues, collapse = "; "))
  gstat <- ga$dirs[["gstat"]]
  info <- read.csv2(file.path(gstat, "GWAS_information.csv"))
  expect_identical(info$acn, inventory$acn)
  expect_identical(info$ancestry, info$acn)
  expect_identical(info$study_design, inventory$study_design)
  expect_identical(info$meta_analysis_method, inventory$meta_analysis_method)
  expect_identical(info$n_cohorts, inventory$n_cohorts)
  expect_identical(dirname(info$recipe_file), rep(scan$recipe_dir, 2))
  expect_false(dir.exists(ga$dirs[["recipes"]]))
  columns <- read.delim(file.path(gstat, "matrix", "studies.tsv"), colClasses = "character")
  expect_identical(columns$ancestry, c("EUR", "MIXED"))
  expect_identical(columns$ld_map[2L], "")
})
