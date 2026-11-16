marker_build_fixture <- function() {
  root <- tempfile("markers38-")
  dir.create(root)
  input <- file.path(root, "input")
  dir.create(input)
  for (chr in seq_len(22L)) {
    markers <- data.frame(rsids = paste0("rs", chr), chr = chr, pos = 100L,
                          ref = "A", alt = "C", af = 0.25)
    if (chr == 1L) {
      extras <- transform(markers[rep(1L, 3L), ], rsids = c("rs91", "rs92", "rs2"),
                          pos = c(101L, 101L, 102L))
      markers <- rbind(markers, extras)
    }
    if (chr == 2L) markers <- rbind(markers, transform(markers, rsids = "rs93", pos = 101L))
    path <- file.path(input, paste0("chr", chr, ".tsv.gz"))
    con <- gzfile(path, "wt")
    write.table(markers, con, sep = "\t", row.names = FALSE, quote = FALSE)
    close(con)
    receipt <- list(build = "GRCh38", chromosome = chr,
      counts = list(retained = nrow(markers)),
      output = list(file = basename(path), md5 = unname(tools::md5sum(path))))
    jsonlite::write_json(receipt, file.path(input, paste0("chr", chr, ".json")),
                         auto_unbox = TRUE)
  }
  list(root = root, input = input)
}

test_that("GRCh38 marker build publishes an aligned catalog and descriptor", {
  x <- marker_build_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  legacy <- file.path(x$root, "hsa.0.0.1")
  dir.create(legacy)
  legacy_file <- file.path(legacy, "untouched")
  writeLines("GRCh37", legacy_file)
  expect_error(gact(version = "hsa.GRCh38.0.0.1", dbdir = x$root), "buildMarkerDB")
  expect_false(dir.exists(file.path(x$root, "hsa.GRCh38.0.0.1")))
  ga <- buildMarkerDB(x$input, x$root)
  expect_identical(ga$version, "hsa.GRCh38.0.0.1")
  expect_identical(ga$metadata$build, "GRCh38")
  expect_identical(ga$assembly, "GRCh38")
  expect_identical(readRDS(file.path(ga$dbdir, "GAlist_hsa.GRCh38.0.0.1.rds")), ga)
  catalog <- data.table::fread(file.path(ga$dirs[["marker"]], "markers.txt.gz"))
  cpra <- data.table::fread(file.path(ga$dirs[["marker"]], "cpra.txt.gz"), header = FALSE)[[1L]]
  expect_identical(catalog$chr, seq_len(22L))
  expect_identical(catalog$rsids, ga$rsids)
  expect_identical(cpra, with(catalog, paste(chr, pos, ea, nea, sep = "_")))
  expect_true(all(catalog$ea == "C" & catalog$nea == "A" & catalog$eaf == 0.25))
  expect_equal(ga$metadata$catalog$n_repeated_position, 2)
  expect_equal(ga$metadata$catalog$n_repeated_rsid, 2)
  expect_equal(ga$metadata$catalog$n_dropped, 4)
  expect_error(buildMarkerDB(x$input, x$root), "already exists")
  expect_error(addLDRegions(ga), "differs")
  expect_identical(readLines(legacy_file), "GRCh37")
  stat <- transform(as.data.frame(catalog), b = 0.2, seb = 0.1, p = 0.05)
  expect_equal(nrow(qcStatDB(ga, stat)), 22)
  expect_error(updateStatDB(ga, stat, build = "GRCh37"), "differs from marker catalog")
})

test_that("catalog build rejects receipts with a wrong build or checksum", {
  x <- marker_build_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  file <- file.path(x$input, "chr1.json")
  receipt <- jsonlite::read_json(file)
  receipt$build <- "GRCh37"
  jsonlite::write_json(receipt, file, auto_unbox = TRUE)
  expect_error(buildMarkerDB(x$input, x$root), "receipt or checksum for chromosome 1")
  receipt$build <- "GRCh38"
  receipt$output$md5 <- "bad"
  jsonlite::write_json(receipt, file, auto_unbox = TRUE)
  expect_error(buildMarkerDB(x$input, x$root), "receipt or checksum for chromosome 1")
  expect_false(dir.exists(file.path(x$root, "hsa.GRCh38.0.0.1")))
})

test_that("database setup lists and accepts only build-specific release names", {
  root <- tempfile("versions-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  invisible(capture.output(versions <- list_gact_versions()))
  expect_identical(versions$version, c("hsa.GRCh37.0.0.1", "hsa.GRCh38.0.0.1"))
  expect_identical(versions$build, c("GRCh37", "GRCh38"))
  for (version in c("hsa.0.0.1", "hsa.0.0.2")) {
    expect_error(gact(version = version, dbdir = root), "not defined")
  }
  expect_length(list.files(root), 0L)
  seen <- character()
  testthat::local_mocked_bindings(downloadDB = function(GAlist, what) {
    expect_identical(GAlist$metadata$build, "GRCh37")
    seen <<- c(seen, what)
    GAlist
  }, .package = "gact")
  ga <- gact(version = "hsa.GRCh37.0.0.1", dbdir = root)
  expect_identical(basename(ga$dbdir), "hsa.GRCh37.0.0.1")
  expect_identical(seen, c("ldetect", "annotation", "marker", "gsea", "gstat",
                           "gbayes", "atc", "dgi"))
})
