test_that("single-study retrieval keeps otherwise complete rows with missing n", {
  file <- tempfile(fileext = ".tsv")
  on.exit(unlink(file))
  stat <- data.frame(rsids = c("rs1", "rs2", "rs3"), b = c(0.2, -0.1, NA),
                     seb = 0.1, p = 0.05, n = c(NA_real_, 120, NA_real_))
  ga <- list(rsids = c(stat$rsids, "rs4"),
             study = list(id = "GWAS1", neff = 100),
             studyfiles = c(GWAS1 = file))
  data.table::fwrite(stat, file, sep = "\t")

  out <- getGWAS(ga, "GWAS1")
  expect_identical(out$rsids, c("rs1", "rs2"))
  expect_equal(out$n, c(NA_real_, 120))
  expect_identical(getGWAS(ga, "GWAS1", rsids = c("rs4", "rs2", "rs1"))$rsids,
                   c("rs2", "rs1"))
  expect_identical(getGWAS(ga, "GWAS1", rm.na = FALSE)$rsids, ga$rsids)

  stat$n <- NULL
  data.table::fwrite(stat, file, sep = "\t")
  expect_equal(getGWAS(ga, "GWAS1")$n, c(100, 100))
  ga$study$neff <- NULL
  out <- getGWAS(ga, "GWAS1")
  expect_identical(out$rsids, c("rs1", "rs2"))
  expect_identical(out$n, rep(NA_real_, 2L))
})
