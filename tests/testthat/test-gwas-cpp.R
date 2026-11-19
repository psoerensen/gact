cpp_gwas_fixture <- function() {
  root <- tempfile("gact-cpp-")
  dir.create(root)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  dirs <- setNames(file.path(root, c("marker", "gstat", "download")),
                   c("marker", "gstat", "download"))
  for (dir in dirs) dir.create(dir)
  marker <- data.frame(rsids = c("rs1", "rs2", "rs3"), chr = 1L, pos = 101:103,
                       ea = c("A", "C", "A"), nea = c("C", "T", "T"), eaf = c(.2, .3, .4))
  data.table::fwrite(marker, file.path(dirs[["marker"]], "markers.txt.gz"))
  data.table::fwrite(data.frame(cpra = with(marker, paste(chr, pos, ea, nea, sep = "_"))),
    file.path(dirs[["marker"]], "cpra.txt.gz"), col.names = FALSE)
  ga <- list(dirs = dirs, rsids = marker$rsids, study = list(id = character()),
             metadata = list(build = "GRCh38"))
  descriptor <- file.path(root, "GAlist.rds")
  saveRDS(ga, descriptor)
  source <- data.frame(chromosome = 1L, base_pair_location = 101:103, effect_allele = c("A", "T", "A"),
    other_allele = c("C", "C", "T"), beta = c(.5, .2, .3), standard_error = .1,
    effect_allele_frequency = c(.2, .7, .4), p_value = .05, rsid = marker$rsids)
  filenames <- c("a.h.tsv.gz", "b.h.tsv.gz")
  files <- file.path(dirs[["download"]], filenames)
  data.table::fwrite(source, files[1L], sep = "\t")
  data.table::fwrite(transform(source, beta = beta * 2), files[2L], sep = "\t")
  mapping <- setNames(as.list(names(source)), c("chr", "pos", "ea", "nea", "b", "seb", "eaf", "p", "marker"))
  knowledge <- file.path(root, "formats.json")
  formats <- list(list(source = "test SSF", headers = names(source), mapping = mapping,
    effect_scale = "beta", p_scale = "p", verified = TRUE, evidence = "synthetic"))
  jsonlite::write_json(list(formats = formats), knowledge, auto_unbox = TRUE)
  manifest <- data.frame(accession = c("GCST123", "GCST456"), filename = filenames, url = "",
    md5 = unname(tools::md5sum(files)), source = "test SSF", trait = c("A", "B"),
    type = "quantitative", build = "GRCh38", ancestry = "Mixed", strand = "forward",
    study_n_total = c(1200L, 1000L))
  list(root = root, ga = ga, descriptor = descriptor, manifest = manifest, knowledge = knowledge,
       files = files, mapping = mapping)
}

test_that("native staging publishes verified recipes and reads both studies", {
  x <- cpp_gwas_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  scratch <- file.path(x$root, "scratch")
  dir.create(scratch)
  writeLines("owned by the caller", file.path(scratch, "keep.me"))
  inputs <- x$files
  inputs[1L] <- file.path(x$root, "upstream.tsv.gz")
  expect_true(file.rename(x$files[1L], inputs[1L]))
  x$manifest$url[1L] <- "https://example.org/a.h.tsv.gz"
  testthat::local_mocked_bindings(curl_download = function(url, destfile, mode) {
    expect_identical(url, x$manifest$url[1L])
    expect_true(file.copy(inputs[1L], destfile))
  }, .package = "curl")
  run <- function(name, import = TRUE) runGWASManifest(x$manifest, x$descriptor,
    out_dir = if (is.null(name)) NULL else file.path(x$root, name),
    knowledge = x$knowledge, engine = "cpp", workers = 2L, import = import,
    scratch_dir = scratch)
  original <- gact:::.publish_gwas_import_file
  fail <- TRUE
  writes <- 0L
  testthat::local_mocked_bindings(.publish_gwas_import_file = function(from, to, ...) {
    if (identical(to, x$descriptor)) {
      work <- setdiff(list.files(scratch, full.names = TRUE), file.path(scratch, "keep.me"))
      expect_length(work, 1L)
      expect_true(file.exists(file.path(work, "marker-index.bin")))
      expect_length(list.files(work, pattern = "\\.txt\\.gz$"), 2L)
      writes <<- writes + 1L
      if (fail) stop("Simulated batch descriptor failure")
    }
    original(from, to, ...)
  }, .package = "gact")
  expect_identical(run("preview", FALSE)$status, rep("ready", 2))
  expect_false(dir.exists(file.path(x$root, "recipes")))
  expect_identical(list.files(scratch), "keep.me")
  expect_true(all(file.exists(x$files)))
  failed <- run("failure")
  expect_identical(failed$status, rep("needs_review", 2))
  expect_match(failed$issues, "Simulated batch descriptor failure", all = TRUE)
  expect_length(list.files(x$ga$dirs[["gstat"]], all.files = TRUE, no.. = TRUE), 0L)
  expect_length(readRDS(x$descriptor)$study$id, 0L)
  expect_identical(list.files(scratch), "keep.me")
  expect_true(all(file.exists(x$files)))
  fail <- FALSE
  writes <- 0L
  result <- run(NULL)
  expect_identical(result$status, rep("imported", 2))
  expect_equal(result$n_flipped, c(1, 1))
  expect_equal(writes, 1L)
  expect_true(all(file.exists(result$recipe)))
  reviews <- normalizePath(file.path(x$root, "recipes", "reviews"), winslash = "/")
  reports <- attr(result, "out_dir")
  expect_identical(dirname(reports), reviews)
  expect_true(all(file.exists(file.path(reports, c("manifest.tsv", "manifest_review.tsv")))))
  expect_length(list.files(reports, recursive = TRUE, pattern = "\\.(gz|bin)$"), 0L)
  expect_identical(list.files(scratch), "keep.me")
  expect_identical(result$cleanup, c("removed", NA_character_))
  expect_identical(file.exists(x$files), c(FALSE, TRUE))
  ga <- readRDS(x$descriptor)
  expect_identical(ga$study$id, c("GWAS1", "GWAS2"))
  expect_equal(ga$study$study_n_total, c(1200, 1000))
  for (i in seq_along(inputs)) {
    source <- data.table::fread(inputs[i], data.table = FALSE)
    normalized <- normalizeStatSchema(source, x$mapping, effect_scale = "beta", p_scale = "p")
    expected <- gact:::qcStatDB(x$ga, normalized, excludeMAFDIFF = .5, strand = "forward")
    actual <- getGWAS(ga, ga$study$id[i], rsids = x$ga$rsids, format = "data.frame")
    expect_equal(lapply(names(expected), function(field) actual[[field]]),
                 lapply(names(expected), function(field) expected[[field]]))
    expect_true(all(is.na(actual$n)))
  }
  expect_identical(run("repeat")$status, rep("already_imported", 2))
})

test_that("native import consumes a saved scan recipe and verifies the full input", {
  x <- cpp_gwas_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  row <- x$manifest[1L, , drop = FALSE]
  row$reference <- paste0("https://www.ebi.ac.uk/gwas/studies/", row$accession)
  profile <- jsonlite::read_json(x$knowledge)$formats[[1L]]
  profile$headers <- unlist(profile$headers, use.names = FALSE)
  preview <- data.table::fread(x$files[1L], data.table = FALSE)
  recipe <- gact:::.gact_gwas_scan_recipe(row, preview, 1024L, profile)
  recipe$verified <- TRUE
  directory <- file.path(x$root, "recipes")
  dir.create(directory)
  path <- gact:::.gact_verified_recipe_path(directory, row$accession, row$md5)
  gact:::.gact_save_verified_recipe(recipe, path, row)
  local_mocked_bindings(inspectGWAS = function(...) stop("Unexpected new preview"), .package = "gact")

  backup <- file.path(x$root, "backup.tsv.gz")
  expect_true(file.copy(x$files[1L], backup))
  writeLines("different file", x$files[1L])
  bad <- runGWASManifest(row, x$descriptor, engine = "cpp")
  expect_identical(bad$status, "needs_review")
  expect_match(bad$issues, "Input MD5 checksum mismatch")
  expect_length(list.files(x$ga$dirs[["gstat"]]), 0L)
  expect_true(file.copy(backup, x$files[1L], overwrite = TRUE))
  good <- runGWASManifest(row, x$descriptor, engine = "cpp")
  expect_identical(good$status, "imported")
  expect_equal(good$n_final, 3L)
})

cpp_n_recipe <- function(x, source, column) {
  data.table::fwrite(source, x$files[1L], sep = "\t")
  row <- x$manifest[1L, , drop = FALSE]
  row$md5 <- unname(tools::md5sum(x$files[1L]))
  row$reference <- paste0("https://www.ebi.ac.uk/gwas/studies/", row$accession)
  profile <- gact:::.gact_catalog_preview_profile(source, row$source)
  recipe <- gact:::.gact_gwas_scan_recipe(row, source, 100L, profile)
  if (!length(recipe$inspection_issues)) {
    recipe$verified <- TRUE
    directory <- file.path(x$root, "recipes")
    dir.create(directory)
    gact:::.gact_save_verified_recipe(recipe,
      gact:::.gact_verified_recipe_path(directory, row$accession, row$md5), row)
  }
  list(row = row, recipe = recipe)
}

test_that("native n survives recipe, alignment, duplicates and default retrieval in any row order", {
  observed <- list()
  for (column in c("n", "N", "N_sample_size", "sample_size")) {
    for (reverse in c(FALSE, TRUE)) {
      x <- cpp_gwas_fixture()
      source <- data.table::fread(x$files[1L], data.table = FALSE)
      source[[column]] <- c(1200, NA_real_, 875)
      source <- rbind(source, source[1L, ])
      if (reverse) source <- source[nrow(source):1L, ]
      saved <- cpp_n_recipe(x, source, column)
      expect_length(saved$recipe$inspection_issues, 0L)
      expect_identical(saved$recipe$mapping$n, column)
      expect_true(is.null(saved$recipe$n)) # Study-level count is not a variant default.
      result <- runGWASManifest(saved$row, x$descriptor, engine = "cpp")
      expect_identical(result$status, "imported")
      ga <- readRDS(x$descriptor)
      actual <- getGWAS(ga, "GWAS1", format = "data.frame")
      expect_identical(actual$rsids, x$ga$rsids)
      expect_equal(actual$n, c(1200, NA, 875))
      expect_equal(actual$b, c(.5, -.2, .3))
      expect_equal(ga$study$study_n_total, 1200)
      expect_true(is.na(ga$study$n))
      qc <- jsonlite::read_json(file.path(x$ga$dirs[["gstat"]], ga$study$qc_file))
      expect_equal(c(qc$n_input, qc$n_final, qc$n_dropped, qc$n_flipped), c(4, 3, 1, 1))
      expect_equal(qc$dropped_by_reason$redundant_duplicate, 1)
      expect_equal(c(qc$n_with_sample_size, qc$n_missing_sample_size), c(2, 1))
      observed[[paste(column, reverse)]] <- actual[c("rsids", "b", "seb", "p", "eaf", "n")]
      unlink(x$root, recursive = TRUE)
    }
  }
  for (actual in observed) expect_equal(actual, observed[[1L]])
})

test_that("native duplicate conflicts include n and missing n", {
  for (conflict in c(1201, NA_real_)) {
    x <- cpp_gwas_fixture()
    source <- data.table::fread(x$files[1L], data.table = FALSE)
    source$n <- c(1200, 1100, 900)
    source <- rbind(source, transform(source[1L, ], n = conflict))
    saved <- cpp_n_recipe(x, source, "n")
    result <- runGWASManifest(saved$row, x$descriptor, engine = "cpp")
    expect_identical(result$status, "imported")
    ga <- readRDS(x$descriptor)
    actual <- getGWAS(ga, "GWAS1", format = "data.frame")
    expect_identical(actual$rsids, c("rs2", "rs3"))
    expect_equal(actual$n, c(1100, 900))
    qc <- jsonlite::read_json(file.path(x$ga$dirs[["gstat"]], ga$study$qc_file))
    expect_equal(c(qc$n_input, qc$n_final, qc$n_dropped), c(4, 2, 2))
    expect_equal(qc$dropped_by_reason$conflicting_duplicate, 2)
    unlink(x$root, recursive = TRUE)
  }
})

test_that("invalid n fails preview and full streaming without publishing a study", {
  x <- cpp_gwas_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  valid <- data.table::fread(x$files[1L], data.table = FALSE)
  valid$n <- c(1200, 1100, 900)
  for (bad in c("0", "-1", "Inf", "NaN", "not-a-count")) {
    preview <- valid
    preview$n <- c(bad, "1100", "900")
    expect_error({
      normalized <- normalizeStatSchema(preview, c(unlist(x$mapping), n = "n"))
      stopifnot(validateStatSchema(normalized)$ok)
    })
    # A clean leading preview does not bypass validation of later full-file rows.
    source <- valid
    source$n <- c("1200", "1100", bad)
    data.table::fwrite(source, x$files[1L], sep = "\t")
    row <- x$manifest[1L, , drop = FALSE]
    row$md5 <- unname(tools::md5sum(x$files[1L]))
    row$reference <- paste0("https://www.ebi.ac.uk/gwas/studies/", row$accession)
    profile <- gact:::.gact_catalog_preview_profile(valid, row$source)
    recipe <- gact:::.gact_gwas_scan_recipe(row, valid[1:2, ], 100L, profile)
    recipe$verified <- TRUE
    directory <- file.path(x$root, "recipes")
    dir.create(directory, showWarnings = FALSE)
    gact:::.gact_save_verified_recipe(recipe,
      gact:::.gact_verified_recipe_path(directory, row$accession, row$md5), row)
    before <- tools::md5sum(x$descriptor)
    result <- runGWASManifest(row, x$descriptor, engine = "cpp")
    expect_identical(result$status, "needs_review")
    expect_match(result$issues, "Per-variant n")
    expect_identical(tools::md5sum(x$descriptor), before)
    expect_length(list.files(x$ga$dirs[["gstat"]], all.files = TRUE, no.. = TRUE), 0L)
  }
})

test_that("wholly missing n stays missing and ambiguous n needs review", {
  x <- cpp_gwas_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  source <- data.table::fread(x$files[1L], data.table = FALSE)
  source$n <- NA_real_
  saved <- cpp_n_recipe(x, source, "n")
  expect_length(saved$recipe$inspection_issues, 0L)
  result <- runGWASManifest(saved$row, x$descriptor, engine = "cpp")
  expect_identical(result$status, "imported")
  ga <- readRDS(x$descriptor)
  actual <- getGWAS(ga, "GWAS1", format = "data.frame")
  expect_equal(nrow(actual), 3L)
  expect_true(all(is.na(actual$n)))
  qc <- jsonlite::read_json(file.path(x$ga$dirs[["gstat"]], ga$study$qc_file))
  expect_equal(c(qc$n_with_sample_size, qc$n_missing_sample_size), c(0, 3))
  source$N <- 1200
  profile <- gact:::.gact_catalog_preview_profile(source, "Catalog")
  expect_match(profile$issues, "Ambiguous per-variant n")
  expect_error(gact:::.gact_cpp_gwas_mapping(list(mapping = c(x$mapping, n = "n"),
    inspection = list(headers = names(source)), effect_scale = "beta", p_scale = "p")),
    "unambiguous")
})
