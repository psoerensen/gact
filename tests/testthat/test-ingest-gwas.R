test_that("detectStatSchema identifies common GWAS columns", {
  stat <- data.frame(
    SNP = paste0("rs", 1:100),
    CHR = rep(1:2, each = 50),
    BP = 1:100,
    A1 = rep(c("A", "C"), 50),
    A2 = rep(c("G", "T"), 50),
    BETA = rnorm(100, 0, 0.1),
    SE = runif(100, 0.01, 0.2),
    P = runif(100, 0, 1),
    N = rep(10000, 100),
    stringsAsFactors = FALSE
  )

  out <- detectStatSchema(stat)

  expect_true(is.list(out))
  expect_equal(out$mapping[["chr"]], "CHR")
  expect_equal(out$mapping[["pos"]], "BP")
  expect_equal(out$mapping[["ea"]], "A1")
  expect_equal(out$mapping[["nea"]], "A2")
  expect_equal(out$mapping[["p"]], "P")
})

test_that("normalizeStatSchema normalizes key fields", {
  stat <- data.frame(
    snp = c("rs1", "rs2"),
    chromosome = c("chrX", "chr2"),
    position = c("100", "200"),
    effect_allele = c("a", "c"),
    other_allele = c("g", "t"),
    pvalue = c(0, 0.5),
    stringsAsFactors = FALSE
  )

  schema <- c(
    marker = "snp",
    chr = "chromosome",
    pos = "position",
    ea = "effect_allele",
    nea = "other_allele",
    p = "pvalue"
  )

  out <- normalizeStatSchema(stat, schema = schema)

  expect_equal(out$chr, c(23L, 2L))
  expect_equal(out$ea, c("A", "C"))
  expect_equal(out$nea, c("G", "T"))
  expect_true(out$p[1] > 0)
})

test_that("validateStatSchema rejects invalid p-values", {
  stat <- data.frame(
    chr = c(1, 1),
    pos = c(100, 200),
    ea = c("A", "C"),
    nea = c("G", "T"),
    p = c(0.1, 1.2)
  )

  out <- validateStatSchema(stat)
  expect_false(out$ok)
  expect_true(any(grepl("Invalid p-values", out$errors)))
})

test_that("ingestStatDB stops on low-confidence required mappings in strict mode", {
  stat <- data.frame(
    chr = c(1, 1),
    pos = c(100, 200),
    ea = c("A", "C"),
    nea = c("G", "T"),
    p = c(0.1, 0.2)
  )

  schema <- list(
    mapping = c(chr = "chr", pos = "pos", ea = "ea", nea = "nea", p = "p"),
    confidence = c(chr = 0.05, pos = 1, ea = 1, nea = 1, p = 1)
  )

  expect_error(
    ingestStatDB(GAlist = list(), stat = stat, schema = schema, strict = TRUE, min_confidence = 0.1),
    "Low-confidence required mappings"
  )
})
