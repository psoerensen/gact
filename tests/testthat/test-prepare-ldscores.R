test_that("legacy GRCh37 tables remain readable and GRCh38 requires provenance", {
  x <- ldscore_fixture("GRCh37")
  on.exit(unlink(x$root, recursive = TRUE))
  path <- file.path(x$root, "markers_1000G_eur_filtered.txt.gz")
  data.table::fwrite(data.frame(rsids = "rs1", chr = 1, ldscores = 1), path)
  expect_equal(getLDscores(x$ga, ancestry = "EUR", version = "1000G"), c(rs1 = 1))
  x$ga$metadata$build <- "GRCh38"
  expect_error(getLDscores(x$ga, ancestry = "EUR", version = "1000G"), "require provenance")
  # Earlier-format tables still have their original names.
  unlink(path)
  expect_error(getLDscores(x$ga, ancestry = "EUR", version = "1000G"), "prepareLD")
})


test_that("selection uses only the requested BED rows and is order independent", {
  x <- ldscore_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  unlink(file.path(x$root, "markers.txt.gz"))
  calls <- 0L
  prep <- gbase::gprep
  testthat::local_mocked_bindings(gprep = function(...) {
    calls <<- calls + 1L
    prep(...)
  }, .package = "gbase")
  full <- prepareLD(x$ga, progress = FALSE)
  expect_identical(calls, 1L)
  sel <- c("rs4", "rs1", "rs2")
  sub <- prepareLD(full, rsids = sel, r2 = 0.00001, progress = FALSE)
  expect_identical(calls, 2L)
  gl <- readRDS(sub$ldPanels$EUR$descriptor)
  expect_identical(gl$rsidsLD[[1]], c("rs1", "rs2", "rs4"))
  scores <- gact::getLDscores(sub)
  expect_equal(scores, gbase::getLDscores(gl))
  expect_equal(unname(scores["rs1"]), 2)
  expect_false(identical(sub$ldPanels$EUR, full$ldPanels$EUR))
  expect_identical(prepareLD(full, rsids = rev(sel), r2 = 1e-5, progress = FALSE), sub)
  expect_identical(calls, 2L)
  expect_length(getLDscores(x$ga, version = "1000G"), 6L)
  expect_error(prepareLD(full, rsids = c("rs1", "unknown"), progress = FALSE), "not found")
  expect_error(prepareLD(full, r2 = -1), "r2")
  expect_error(prepareLD(full, window_cm = 0), "window_cm")
})
