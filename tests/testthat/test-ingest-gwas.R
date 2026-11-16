test_that("detectStatSchema identifies common GWAS columns", {
  set.seed(20261005)
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
  expect_equal(out$mapping[["marker"]], "SNP")
  expect_equal(out$mapping[["seb"]], "SE")
  expect_equal(out$mapping[["n"]], "N")
})

test_that("detectStatSchema recognizes effect_se in harmonized GWAS files", {
  stat <- data.frame(
    rsID = c("rs1", "rs2"), chr = c(1L, 1L), position = c(100L, 200L),
    N_sample_size = c(10000, 10000),
    effect_allele = c("A", "C"), non_effect_allele = c("G", "T"),
    effect_size = c(-0.2, 0.3), effect_se = c(0.04, 0.08), pvalue = c(0.01, 0.02)
  )
  schema <- detectStatSchema(stat)
  expect_identical(schema$mapping[["seb"]], "effect_se")
  expect_identical(normalizeStatSchema(stat, schema)$seb, stat$effect_se)

  stat$effect_se <- NA_real_
  schema <- detectStatSchema(stat)
  expect_identical(schema$mapping[["seb"]], "effect_se")
  expect_error(normalizeStatSchema(stat, schema), "Missing required canonical fields.*seb")
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

  out <- normalizeStatSchema(stat, schema = schema, require_fields = names(schema))

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

test_that("ingestion requires beta and standard error columns", {
  stat <- data.frame(chr = 1L, pos = 100L, ea = "A", nea = "G", p = 0.1)
  expect_match(validateStatSchema(stat)$errors, "b, seb")
  expect_error(ingestStatDB(list(), stat, schema = c(p = "p")), "b, seb")
  stat$b <- 0.2
  stat$seb <- NA_real_
  expect_match(validateStatSchema(stat)$errors, "Empty required columns: seb")
  stat$p <- NULL
  stat$eaf <- 0.1
  expect_true(is.na(detectStatSchema(stat)$mapping[["p"]]))
})

test_that("odds-ratio conversion follows the mapped name or explicit scale", {
  stat <- data.frame(chr = 1L, pos = 100:101, ea = "A", nea = "G",
    p = 0.1, beta = c(0.75, 1.25), odds_ratio = c(0.2, 0.4), seb = 0.1)

  beta <- normalizeStatSchema(stat, c(b = "beta"))
  expect_identical(beta$b, stat$beta)
  expect_null(attr(beta, "transformations"))

  expect_warning(or <- normalizeStatSchema(stat, c(b = "odds_ratio")), "SE may not")
  expect_equal(or$b, log(stat$odds_ratio))
  expect_warning(explicit <- normalizeStatSchema(stat, c(b = "beta"), effect_scale = "or"),
    "SE may not")
  expect_equal(explicit$b, log(stat$beta))
  expect_equal(normalizeStatSchema(stat, c(b = "odds_ratio"), effect_scale = "beta")$b,
    stat$odds_ratio)
  stat$odds_ratio[1L] <- 0
  expect_error(normalizeStatSchema(stat, c(b = "odds_ratio")), "strictly positive")
})

test_that("validation rejects unusable beta, SE, and EAF values", {
  stat <- data.frame(chr = 1L, pos = 100:105, ea = "A", nea = "G", p = 0.1,
    b = c(0, -0.1, Inf, -Inf, NaN, NA_real_),
    seb = c(0.1, 0, -0.1, Inf, NaN, NA_real_),
    eaf = c(0, 1, -0.1, 1.1, Inf, NA_real_))

  validation <- validateStatSchema(stat)
  expect_false(validation$ok)
  expect_match(validation$errors[1], "Non-finite beta values \\(b\\): 3")
  expect_match(validation$errors[2], "standard errors|Standard errors")
  expect_match(validation$errors[2], "strictly positive: 4")
  expect_match(validation$errors[3], "outside \\[0,1\\]: 3")
  expect_true(validateStatSchema(stat[c(1, 6), ])$ok)
  expect_true(validateStatSchema(transform(stat[1:2, ], seb = 0.1))$ok)
  expect_error(ingestStatDB(list(), transform(stat, eaf = c("0", "bad", rep(NA, 4))),
    schema = c(p = "p")), "Non-numeric values in eaf")
})

test_that("p values are converted by the mapped column name or explicit scale", {
  stat <- data.frame(chr = 1L, pos = 100:104, ea = "A", nea = "G", b = 0.1, seb = 0.1,
    neglog10p = c(0, 0.3, 0.7, 400, NA_real_))
  expected <- c(1, 10^-0.3, 10^-0.7, .Machine$double.xmin, NA_real_)
  norm <- normalizeStatSchema(stat, c(p = "neglog10p"))
  expect_equal(norm$p, expected)
  expect_identical(attr(norm, "transformations"), "-log10(p) -> p")

  names(stat)[7] <- "Score"
  expect_equal(normalizeStatSchema(stat, c(p = "Score"), p_scale = "neglog10p")$p, expected)
  expect_false(validateStatSchema(normalizeStatSchema(stat, c(p = "Score")))$ok)

  ordinary <- transform(stat, neglog10p = c(0.8, 0.3, 0.7, 0.2, NA_real_))
  expect_identical(normalizeStatSchema(ordinary, c(p = "neglog10p"), p_scale = "p")$p,
    ordinary$neglog10p)
  expect_error(normalizeStatSchema(transform(stat, Score = -0.1), c(p = "Score"),
    p_scale = "neglog10p"), "finite and nonnegative")
  expect_error(normalizeStatSchema(transform(stat, Score = Inf), c(p = "Score"),
    p_scale = "neglog10p"), "finite and nonnegative")
  expect_error(normalizeStatSchema(stat, c(b = "b")), "Missing required canonical fields.*p")
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

test_that("alias dictionary boosts schema mapping", {
  stat <- data.frame(
    chromosome_name = c(1, 1, 2, 2),
    base_pair_location = c(101, 102, 201, 202),
    effect_allele = c("A", "C", "A", "G"),
    other_allele = c("G", "T", "C", "A"),
    p_value = c(0.1, 0.2, 0.01, 0.9),
    stringsAsFactors = FALSE
  )

  alias_dict <- list(
    version = 1,
    fields = list(
      chr = list(aliases = c("chromosomename"), regex = character(0)),
      pos = list(aliases = c("basepairlocation"), regex = character(0)),
      ea = list(aliases = c("effectallele"), regex = character(0)),
      nea = list(aliases = c("otherallele"), regex = character(0)),
      p = list(aliases = c("pvalue"), regex = character(0))
    )
  )

  out <- detectStatSchema(stat, alias_dict = alias_dict)
  expect_equal(out$mapping[["chr"]], "chromosome_name")
  expect_equal(out$mapping[["pos"]], "base_pair_location")
  expect_equal(out$mapping[["p"]], "p_value")
  expect_true(length(out$alias_hits[["chr"]]) > 0)
})
