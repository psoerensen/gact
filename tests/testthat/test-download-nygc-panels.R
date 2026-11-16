test_that("GRCh38 downloads install verified NYGC triplets and preserve existing files", {
  root <- tempfile("nygc-download-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  ga <- createDB("hsa.GRCh38.0.0.1", root)
  ga$metadata$sources <- list()
  marker <- ga$dirs[["marker"]]
  writeLines("catalog", file.path(marker, "markers.txt.gz"))
  scratch <- file.path(root, "archives")
  dir.create(scratch)
  ancestries <- c("afr", "amr", "eas", "eur", "sas")
  files <- as.vector(outer(paste0("g1000_", ancestries), c(".bed", ".bim", ".fam"), paste0))
  archives <- paste0("gact-nygc30x-grch38-", ancestries, ".tar.gz")
  checksums <- (function() {
    old <- setwd(scratch)
    on.exit(setwd(old))
    for (name in files) writeLines(name, name)
    for (i in seq_along(archives)) {
      members <- paste0("g1000_", ancestries[i], c(".bed", ".bim", ".fam"))
      utils::tar(archives[i], members, compression = "gzip", tar = "internal")
    }
    stats::setNames(vapply(c(files, archives), gact:::.gact_sha256_file, character(1)),
                    c(files, archives))
  })()
  seen <- character()
  testthat::local_mocked_bindings(
    .gact_nygc_release_hashes = function(record) checksums,
    .gact_download_nygc_archive = function(record, archive, dest) {
      expect_identical(record, "24292610")
      seen <<- c(seen, archive)
      stopifnot(file.copy(file.path(scratch, archive), dest))
    }, .package = "gact"
  )
  checksums[archives[5L]] <- strrep("0", 64L)
  expect_error(downloadDB(ga, "1000G"), "SHA-256 mismatch.*sas")
  expect_identical(list.files(marker, all.files = TRUE, no.. = TRUE), "markers.txt.gz")
  checksums[archives[5L]] <- gact:::.gact_sha256_file(file.path(scratch, archives[5L]))
  checksums[files[1L]] <- strrep("0", 64L)
  expect_error(downloadDB(ga, "1000G"), "SHA-256 mismatch for extracted file")
  expect_identical(list.files(marker, all.files = TRUE, no.. = TRUE), "markers.txt.gz")
  checksums[files[1L]] <- gact:::.gact_sha256_file(file.path(scratch, files[1L]))

  expect_true(file.copy(file.path(scratch, files[1L]), marker))
  before <- file.info(file.path(marker, files[1L]))$mtime
  seen <- character()
  installed <- downloadDB(ga, "1000G")
  expect_setequal(list.files(marker), c(files, "markers.txt.gz"))
  expect_identical(seen, archives)
  expect_identical(file.info(file.path(marker, files[1L]))$mtime, before)
  expect_identical(installed$downloads$status, "success")
  expect_identical(installed$downloads$dest_file, marker)
  expect_identical(readLines(file.path(marker, "markers.txt.gz")), "catalog")
  seen <- character()
  downloadDB(ga, "1000G")
  expect_length(seen, 0L)
  writeLines("conflict", file.path(marker, files[1L]))
  expect_error(downloadDB(ga, "1000G"), "Existing NYGC destination differs")
  expect_identical(readLines(file.path(marker, files[1L])), "conflict")
})

test_that("GRCh37 continues to download its legacy Zenodo records", {
  root <- tempfile("nygc-grch37-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  ga <- createDB("hsa.GRCh37.0.0.1", root)
  seen <- character()
  testthat::local_mocked_bindings(download_zenodo = function(doi, path, parallel) {
    expect_identical(unname(path), ga$dirs[["marker"]])
    seen <<- c(seen, doi)
  }, .package = "gact")
  downloadDB(ga, "1000G")
  expect_identical(seen, paste0("10.5281/zenodo.",
    c(10462403, 14141438, 14141523, 14141597, 14141774, 14141885)))
})
