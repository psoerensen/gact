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
  expect_identical(ga$dirs[["reference"]], file.path(ga$dbdir, "reference"))
  expect_true(dir.exists(ga$dirs[["reference"]]))
  prepared <- file.path(ga$dirs[["reference"]], "prepared_grch38")
  files <- list.files(x$input)
  expect_setequal(list.files(prepared), files)
  expect_identical(unname(tools::md5sum(file.path(prepared, files))),
                   unname(tools::md5sum(file.path(x$input, files))))
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

test_that("reference packaging uses release inputs and supports the old build cache", {
  x <- marker_build_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  cache <- file.path(x$root, "reference")
  lookup_dir <- file.path(cache, "gnomad-rsids")
  dir.create(lookup_dir, recursive = TRUE)
  paths <- file.path(lookup_dir, gact:::.gact_gnomad_lookup_files(1L))
  con <- gzfile(paths[1L], "wt")
  writeLines(c("#CHROM\tPOS\tID\tREF\tALT", "chr1\t100\trs1\tA\tC"), con)
  close(con)
  hash <- unname(tools::md5sum(paths[1L]))
  upstream <- list(version = "4.1", chromosome = 1L, url = "https://example.org/chr1.vcf.bgz",
                   source_md5 = strrep("a", 32L), output_md5 = hash)
  jsonlite::write_json(upstream, paths[2L], auto_unbox = TRUE)
  receipt_file <- file.path(x$input, "chr1.json")
  receipt <- jsonlite::read_json(receipt_file)
  receipt$sources <- list(list(), list(), list(file = basename(paths[1L]), md5 = hash,
                                              upstream = upstream))
  jsonlite::write_json(receipt, receipt_file, auto_unbox = TRUE)
  ga <- buildMarkerDB(x$input, x$root)
  reference <- ga$dirs[["reference"]]
  copied <- file.path(reference, "gnomad-rsids", basename(paths))
  expect_identical(unname(tools::md5sum(copied)), unname(tools::md5sum(paths)))
  preserved <- list.files(reference, recursive = TRUE, full.names = TRUE)
  prior <- tools::md5sum(preserved)
  bundles <- packageReference(ga)
  expect_identical(names(bundles), c("tables", "gnomad"))
  expect_identical(dirname(bundles$tables$archive), reference)
  expect_identical(tools::md5sum(preserved), prior)
  expect_setequal(utils::untar(bundles$tables$archive, list = TRUE, tar = "internal"),
                  list.files(x$input))
  expect_match(readLines(bundles$tables$checksum_file), bundles$tables$sha256)
  installed <- installGnomadLookups(bundles$gnomad$archive, file.path(x$root, "installed"))
  expect_identical(unname(tools::md5sum(file.path(installed, basename(paths)))),
                   unname(tools::md5sum(paths)))
  expect_error(packageReference(ga), "already exists")

  # Simulate a release built before reference inputs were copied into it.
  legacy_tables <- file.path(cache, "grch38-prep.auto", "tables")
  dir.create(dirname(legacy_tables))
  expect_true(file.rename(x$input, legacy_tables))
  unlink(file.path(reference, c("prepared_grch38", "gnomad-rsids")), recursive = TRUE)
  legacy_ga <- ga
  legacy_ga$dirs <- ga$dirs[names(ga$dirs) != "reference"]
  from_cache <- packageReference(legacy_ga, out_dir = file.path(x$root, "legacy-artifacts"))
  expect_identical(names(from_cache), c("tables", "gnomad"))
  writeLines("corrupt", file.path(legacy_tables, "chr1.tsv.gz"))
  expect_error(packageReference(ga, out_dir = file.path(x$root, "bad")),
               "differ from the built catalog")
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

test_that("marker release packaging preserves source files and verifies archive contents", {
  x <- marker_build_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  ga <- buildMarkerDB(x$input, x$root)
  paths <- list.files(ga$dbdir, full.names = TRUE, recursive = TRUE)
  before <- tools::md5sum(paths)
  out <- file.path(x$root, "release")
  release <- packageMarkerDB(ga$dbdir, out)
  expect_identical(tools::md5sum(paths), before)
  expect_setequal(list.files(out), basename(c(release$archive, release$checksum_file)))
  expect_identical(readLines(release$checksum_file),
    paste(gact:::.gact_sha256_file(release$archive), basename(release$archive), sep = "  "))
  members <- c("marker/markers.txt.gz", "marker/cpra.txt.gz", "marker/marker_manifest.json")
  expect_setequal(utils::untar(release$archive, list = TRUE, tar = "internal"), c(members, "SHA256SUMS"))
  unpack <- file.path(x$root, "unpack")
  expect_equal(utils::untar(release$archive, exdir = unpack, tar = "internal"), 0)
  hashes <- vapply(file.path(ga$dbdir, members), gact:::.gact_sha256_file, "", USE.NAMES = FALSE)
  expect_identical(readLines(file.path(unpack, "SHA256SUMS")), paste(hashes, members, sep = "  "))
  expect_error(packageMarkerDB(ga$dbdir, out), "already exists")
  expect_error(packageMarkerDB(ga$dbdir, file.path(ga$dbdir, "release")), "outside")
  writeLines("altered", file.path(ga$dirs[["marker"]], "cpra.txt.gz"))
  bad <- file.path(x$root, "bad")
  expect_error(packageMarkerDB(ga$dbdir, bad), "checksums do not match")
  expect_length(list.files(bad, all.files = TRUE, no.. = TRUE), 0L)
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
