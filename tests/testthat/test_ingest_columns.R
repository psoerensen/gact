library(testthat)

# -------------------------------
# Helper: run full pipeline
# -------------------------------
run_pipeline <- function(df) {
 schema <- detectStatSchema(df)
 norm <- normalizeStatSchema(df, schema)
 val <- validateStatSchema(norm)
 list(schema = schema, norm = norm, val = val)
}

# -------------------------------
# UK Biobank (BOLT-LMM style)
# -------------------------------
test_that("UKB-style data works", {

 df <- data.frame(
  SNP = paste0("rs", 1:1000),
  CHR = sample(1:22, 1000, TRUE),
  BP = sample(1e6:1e7, 1000),
  A1 = sample(c("A","C","G","T"), 1000, TRUE),
  A2 = sample(c("A","C","G","T"), 1000, TRUE),
  BETA = rnorm(1000, 0, 0.1),
  SE = runif(1000, 0.01, 0.2),
  P_BOLT_LMM = runif(1000),
  N = sample(50000:100000, 1000, TRUE)
 )

 res <- run_pipeline(df)

 expect_true(res$val$ok)
 expect_true(all(c("chr","pos","ea","nea","p") %in% names(res$norm)))
})

# -------------------------------
# FinnGen-style (OR + -log10(p))
# -------------------------------
test_that("FinnGen-style OR and mlogp conversion works", {

 df <- data.frame(
  rsid = paste0("rs", 1:1000),
  chrom = sample(1:22, 1000, TRUE),
  pos = sample(1e6:1e7, 1000),
  alt = sample(c("A","C","G","T"), 1000, TRUE),
  ref = sample(c("A","C","G","T"), 1000, TRUE),
  odds_ratio = runif(1000, 0.5, 1.5),
  neglog10p = runif(1000, 10, 50),
  af = runif(1000),
  n = sample(10000:50000, 1000, TRUE)
 )

 #res <- run_pipeline(df)
 #expect_warning(
 # run_pipeline(df),
 # "SE may not correspond"
 #)

 expect_warning(
  res <- run_pipeline(df),
  "SE may not correspond"
 )

 # Check conversion happened
 expect_true(all(res$norm$p <= 1))
 expect_true(any(res$norm$b < 0))  # log(OR) introduces negatives

 expect_true(res$val$ok)

 # Check metadata
 tr <- attr(res$norm, "transformations")
 expect_true(any(grepl("log", tr)))
 expect_true(any(grepl("p", tr)))
})

# -------------------------------
# SAIGE-style
# -------------------------------
test_that("SAIGE-style data works", {

 df <- data.frame(
  MarkerID = paste0("rs", 1:1000),
  CHR = sample(1:22, 1000, TRUE),
  POS = sample(1e6:1e7, 1000),
  Allele1 = sample(c("A","C","G","T"), 1000, TRUE),
  Allele2 = sample(c("A","C","G","T"), 1000, TRUE),
  AF_Allele2 = runif(1000),
  BETA = rnorm(1000, 0, 0.2),
  SE = runif(1000, 0.01, 0.3),
  p.value = runif(1000),
  N = sample(10000:50000, 1000, TRUE)
 )

 res <- run_pipeline(df)

 expect_true(res$val$ok)
 expect_false(any(is.na(res$norm$ea)))
 expect_false(any(is.na(res$norm$nea)))
})

# -------------------------------
# Edge case: missing columns
# -------------------------------
test_that("Missing required fields triggers error", {

 df <- data.frame(
  SNP = paste0("rs", 1:100),
  BETA = rnorm(100)
 )

 expect_error(
  normalizeStatSchema(df, detectStatSchema(df)),
  "Missing required"
 )
})

# -------------------------------
# Edge case: ambiguous columns
# -------------------------------
test_that("Ambiguous mapping produces warning", {

 df <- data.frame(
  pval = runif(100),
  pvalue = runif(100),
  chr = sample(1:22, 100, TRUE),
  pos = sample(1e6:1e7, 100),
  ea = sample(c("A","C","G","T"), 100, TRUE),
  nea = sample(c("A","C","G","T"), 100, TRUE)
 )

 schema <- detectStatSchema(df)

 expect_true(any(grepl("Ambiguous", schema$warnings)))
})

# -------------------------------
# Edge case: allele validation
# -------------------------------
test_that("Non-ACGT alleles produce warnings", {

 df <- data.frame(
  chr = 1,
  pos = 12345,
  ea = "DEL",
  nea = "A",
  p = 0.05
 )

 val <- validateStatSchema(df)

 expect_true(length(val$warnings) > 0)
})

# -------------------------------
# Robustness: large dataset sampling
# -------------------------------
test_that("Sampling does not break detection", {

 df <- data.frame(
  SNP = paste0("rs", 1:200000),
  CHR = sample(1:22, 200000, TRUE),
  BP = sample(1e6:1e7, 200000, TRUE),
  A1 = sample(c("A","C","G","T"), 200000, TRUE),
  A2 = sample(c("A","C","G","T"), 200000, TRUE),
  BETA = rnorm(200000),
  SE = runif(200000, 0.01, 0.2),
  P = runif(200000)
 )

 schema <- detectStatSchema(df, sample_n = 5000)

 expect_true(!all(is.na(schema$mapping)))
})
