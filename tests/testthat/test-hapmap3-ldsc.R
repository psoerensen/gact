test_that("local LDSC uses corrected reference and regression sums and a panel M", {
  x <- ldscore_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  hm <- c("rs1", "rs4", "rs5")
  ga <- prepareLDSC(x$ga, hm3_rsids = hm)
  file <- ga$ldscPanels$EUR$file
  result <- readRDS(file)
  expect_identical(result$rows$rsids, hm)
  # rs1 and rs2 have identical genotypes, rs4 has one missing observation.
  # rs2 tags rs4 inside 1 cM, while rs1 is out of range of rs4.
  a <- c(0, 1, 2, 0)
  b <- c(1.0, 1, 2, 0)
  r2 <- cor(a, b)^2
  expect_equal(result$rows$ld_score, c(2, 1 + r2 - (1-r2)/2, 1))
  expect_equal(result$rows$weight_ld, c(1, 1, 1))
  expect_equal(result$marker_count, 5)
  expect_equal(getLDscores(ga, version = "HapMap3", chr = 2), c(rs5 = 1))
  stamp <- file.mtime(file)
  expect_equal(prepareLDSC(x$ga, hm3_rsids = rev(hm))$ldscPanels, ga$ldscPanels)
  expect_identical(file.mtime(file), stamp)
  expect_error(prepareLDSC(ga, hm3_rsids = c(hm, "rs2")), "inputs changed")
})

test_that("MAF boundaries, nonvariable SNPs and the extended MHC have distinct roles", {
  x <- ldscore_fixture("GRCh37")
  on.exit(unlink(x$root, recursive = TRUE))
  # 50 diploid individuals: one, two, five and six minor alleles in 100 copies.
  het <- function(k) c(rep(2L, k), rep(0L, 50L-k))
  # Reverse the allele coding at the 5% boundary (PLINK homozygous code 3).
  code <- cbind(het(1), het(2), c(rep(2L,5), rep(3L,45)), het(6), het(25), rep(2L,50))
  ids <- paste0("rs", seq_len(ncol(code)))
  x$ga$rsids <- ids
  writeLines(paste("f", paste0("i", 1:50), 0, 0, 0, -9), paste0(x$prefix,".fam"))
  writeLines(paste(6,ids,0:5,c(1:4,26000000L,28000000L),"A","G"), paste0(x$prefix,".bim"))
  bytes <- apply(code, 2, function(z) as.raw(crossprod(c(1,4,16,64),
    matrix(c(z,0L,0L),nrow=4L))))
  writeBin(c(as.raw(c(108,27,1)),as.raw(bytes)),paste0(x$prefix,".bed"))
  ref <- readRDS(prepareLDSC(x$ga, hm3_rsids=ids)$ldscPanels$EUR$file)
  expect_identical(ref$rows$rsids, ids[2:4])
  expect_equal(ref$marker_count, 2)
  expect_identical(ref$settings$build, "GRCh37")
})
