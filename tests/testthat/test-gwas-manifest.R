manifest_fixture <- function() {
  root <- tempfile("gwas-manifest-")
  dir.create(root)
  dirs <- stats::setNames(file.path(root, c("download", "marker", "gstat")),
                         c("download", "marker", "gstat"))
  for (dir in dirs) dir.create(dir)
  stat <- data.frame(chr = 1L, pos = 101:104, ea = "A", nea = "C",
                     b = 0.2, seb = 0.1, p = 0.05, eaf = 0.2)
  markers <- transform(stat[c("chr", "pos", "ea", "nea", "eaf")],
                       rsids = paste0("rs", stat$pos))
  data.table::fwrite(markers, file.path(dirs[["marker"]], "markers.txt.gz"))
  cpra <- with(markers, paste(chr, pos, ea, nea, sep = "_"))
  data.table::fwrite(data.frame(cpra), file.path(dirs[["marker"]], "cpra.txt.gz"),
                     col.names = FALSE)
  ga <- list(metadata = list(build = "GRCh37"), dirs = dirs, rsids = markers$rsids)
  ga <- updateStatDB(ga, stat, build = "GRCh37",
                     reference = "https://www.ebi.ac.uk/gwas/studies/GCST100/")
  descriptor <- file.path(root, "GAlist.rds")
  saveRDS(ga, descriptor)

  filename <- c("first.tsv", "second.tsv")
  paths <- file.path(dirs[["download"]], filename)
  for (path in paths) data.table::fwrite(stat, path, sep = "\t")
  manifest <- data.frame(accession = c("GCST100", "GCST200"), filename = filename,
    url = "", md5 = unname(tools::md5sum(paths)), trait = c("existing", "new"))
  defaults <- list(source = "test release", type = "quantitative",
                   ancestry = "EUR", build = "GRCh37")
  book <- file.path(root, "formats.json")
  jsonlite::write_json(list(formats = list(list(source = defaults$source,
    headers = names(stat), mapping = as.list(stats::setNames(names(stat), names(stat))),
    effect_scale = "beta", p_scale = "p", verified = TRUE, evidence = "Synthetic fixture"))),
    book, auto_unbox = TRUE)
  list(root = root, dirs = dirs, descriptor = descriptor, manifest = manifest,
       defaults = defaults, book = book)
}

test_that("manifest reuses verified formats, persists accession, and skips reruns", {
  x <- manifest_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  data <- x$dirs[["download"]]
  before <- tools::md5sum(x$descriptor)
  args <- list(manifest = x$manifest, descriptor = x$descriptor, data_dir = data,
               knowledge = x$book, defaults = x$defaults)
  preview <- do.call(runGWASManifest, c(args, list(
    out_dir = file.path(x$root, "preview"), import = FALSE)))
  expect_identical(preview$status, c("already_imported", "ready"))
  expect_true(jsonlite::read_json(preview$recipe[2])$inspection$profile_reused)
  expect_identical(tools::md5sum(x$descriptor), before)
  batch <- file.path(x$root, "batch")
  run <- do.call(runGWASManifest, c(args, list(out_dir = batch)))
  expect_identical(run$status, c("already_imported", "imported"))
  expect_equal(run$n_final[2], 4)
  expect_true(file.exists(file.path(batch, "manifest_review.tsv")))
  updated <- readRDS(x$descriptor)
  expect_identical(updated$study$accession, c(NA_character_, "GCST200"))
  expect_identical(updated$study$reference[2], "https://www.ebi.ac.uk/gwas/studies/GCST200")
  expect_equal(nrow(getGWAS(updated, "GWAS2", what = "p")), 4)
  unlink(file.path(data, x$manifest$filename))
  before <- tools::md5sum(x$descriptor)
  again <- do.call(runGWASManifest, c(args, list(out_dir = file.path(x$root, "again"))))
  expect_identical(again$status, rep("already_imported", 2))
  expect_identical(again$study_id, c("GWAS1", "GWAS2"))
  expect_identical(tools::md5sum(x$descriptor), before)
})

test_that("manifest reports review reasons and preserves the database", {
  x <- manifest_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  x$manifest$accession[1] <- "GCST300"
  x$manifest$md5[1] <- paste(rep("0", 32), collapse = "")
  before <- tools::md5sum(c(x$descriptor, list.files(x$dirs[["gstat"]], full.names = TRUE)))
  run <- runGWASManifest(x$manifest, x$descriptor, x$dirs[["download"]],
    file.path(x$root, "review"), defaults = x$defaults)
  expect_identical(run$status, rep("needs_review", 2))
  expect_match(run$issues[1], "checksum mismatch")
  expect_match(run$issues[2], "Review inferred mapping")
  wrong <- transform(x$manifest[2, ], build = "GRCh38")
  run <- runGWASManifest(wrong, x$descriptor, x$dirs[["download"]],
    file.path(x$root, "build"), knowledge = x$book, defaults = x$defaults)
  expect_match(run$issues, "same genome build")
  expect_identical(tools::md5sum(names(before)), before)
})

test_that("a descriptor publication failure leaves a manifest import safe to retry", {
  x <- manifest_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  original <- gact:::.publish_gwas_import_file
  fail <- TRUE
  testthat::local_mocked_bindings(.publish_gwas_import_file = function(from, to, ...) {
    if (fail && identical(basename(to), "GAlist.rds")) stop("Simulated descriptor write failure")
    original(from, to, ...)
  }, .package = "gact")
  paths <- c(x$descriptor, list.files(x$dirs[["gstat"]], full.names = TRUE))
  before <- tools::md5sum(paths)
  run <- runGWASManifest(x$manifest, x$descriptor, x$dirs[["download"]],
    file.path(x$root, "failure"), knowledge = x$book, defaults = x$defaults)
  expect_identical(run$status, c("already_imported", "needs_review"))
  expect_match(run$issues[2], "Simulated descriptor write failure")
  expect_identical(tools::md5sum(paths), before)
  expect_setequal(list.files(x$dirs[["gstat"]], full.names = TRUE, all.files = TRUE, no.. = TRUE),
                  paths[-1])
  fail <- FALSE
  retry <- runGWASManifest(x$manifest, x$descriptor, x$dirs[["download"]],
    file.path(x$root, "retry"), knowledge = x$book, defaults = x$defaults)
  expect_identical(retry$status, c("already_imported", "imported"))
  expect_identical(retry$study_id[2], "GWAS2")
})

test_that("downloads are published to the cache only after checksum validation", {
  x <- manifest_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  data <- x$dirs[["download"]]
  row <- x$manifest[2, ]
  row$url <- "https://example.org/study.tsv"
  cached <- file.path(data, row$filename)
  source <- file.path(x$root, "source.tsv")
  expect_true(file.rename(cached, source))
  testthat::local_mocked_bindings(curl_download = function(url, destfile, mode) {
    expect_identical(url, row$url)
    expect_identical(mode, "wb")
    expect_true(file.copy(source, destfile))
  }, .package = "curl")
  run <- runGWASManifest(transform(row, md5 = paste(rep("0", 32), collapse = "")),
    x$descriptor, data, file.path(x$root, "bad-download"), defaults = x$defaults)
  expect_match(run$issues, "Input MD5 checksum mismatch")
  expect_false(file.exists(cached))
  expect_setequal(list.files(data, all.files = TRUE, no.. = TRUE), x$manifest$filename[1])
  run <- runGWASManifest(row, x$descriptor, data, file.path(x$root, "downloaded"),
    knowledge = x$book, defaults = x$defaults, import = FALSE)
  expect_identical(run$status, "ready")
  expect_identical(unname(tools::md5sum(cached)), row$md5)
})
