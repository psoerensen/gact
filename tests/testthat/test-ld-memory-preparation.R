test_that("a small positive r2 removes a weak nonzero pair from storage and scores", {
  x <- ldscore_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  n <- 1000L
  a <- rep(c(0L, 2L), n / 2)
  b <- rep(c(0L, 0L, 2L, 2L), n / 4)
  b[2] <- 1L
  pack <- function(v) {
    code <- c(3L, 2L, 0L)[v + 1L]
    as.raw(colSums(matrix(code, nrow = 4L) * c(1L, 4L, 16L, 64L)))
  }
  writeBin(c(as.raw(c(0x6c, 0x1b, 1)), pack(a), pack(b)), paste0(x$prefix, ".bed"))
  writeLines(paste("f", seq_len(n), 0, 0, 0, -9), paste0(x$prefix, ".fam"))
  writeLines(c("1 rs1 0 10 A G", "1 rs2 1 2000000 A G"), paste0(x$prefix, ".bim"))
  full <- prepareLD(x$ga, progress = FALSE)
  cut <- prepareLD(x$ga, r2 = 0.00001, progress = FALSE)
  expected <- cor(a, b)^2
  expect_gt(expected, 0)
  expect_lt(expected, 0.00001)
  expect_equal(unname(getLDscores(full, version = "1000G")), rep(1 + expected, 2), tolerance = 1e-10)
  expect_equal(unname(getLDscores(cut, version = "1000G")), c(1, 1))
  expect_equal(readRDS(full$ldPanels$EUR$descriptor)$sparseLD$nnz, 1)
  expect_equal(readRDS(cut$ldPanels$EUR$descriptor)$sparseLD$nnz, 0)
})

test_that("unusable genetic maps fail before producing an LD resource", {
  x <- ldscore_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  bim <- paste0(x$prefix, ".bim")
  rows <- read.table(bim)
  rows$V3 <- 0
  write.table(rows, bim, row.names = FALSE, col.names = FALSE, quote = FALSE)
  expect_error(prepareLD(x$ga, progress = FALSE), "cM map is all zero")
  rows$V3 <- c(0, 1, 0.5, 2, 0, 1)
  write.table(rows, bim, row.names = FALSE, col.names = FALSE, quote = FALSE)
  expect_error(prepareLD(x$ga, progress = FALSE), "nondecreasing")
  expect_length(list.files(x$root, pattern = "row_ptr|[.]lock$"), 0L)
})
