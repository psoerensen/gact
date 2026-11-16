test_that("successive studies persist complete checksum and recipe metadata", {
  root <- tempfile("study-metadata-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  local_mocked_bindings(qcStatDB = function(GAlist, stat, ...) stat,
                        .package = "gact")
  ga <- list(dirs = c(gstat = root), study = list(id = character()))
  stat <- data.frame(chr = 1L, pos = 10L, ea = "A", nea = "G",
    eaf = 0.2, b = 0.1, seb = 0.05, p = 0.04)
  for (i in seq_len(3)) {
    ga <- updateStatDB(ga, stat, trait = paste0("trait", i),
                       recipe = if (i == 1) list(source = "fixture") else NULL)
    saved <- read.csv2(file.path(root, "GWAS_information.csv"))
    expect_equal(nrow(saved), i)
    expect_identical(saved$id, paste0("GWAS", seq_len(i)))
    expect_identical(saved$checksum, unname(tools::md5sum(ga$studyfiles)))
    expect_equal(unname(lengths(ga$study)), rep(i, length(ga$study)))
    expect_true(all(is.na(saved$n)))
    if (i > 1) expect_true(is.na(saved$recipe_file[i]))
  }
  ga <- updateStatDB(ga, writeStatDB = FALSE)
  saved <- read.csv2(file.path(root, "GWAS_information.csv"))
  expect_equal(nrow(saved), 4)
  expect_true(is.na(saved$checksum[4]))
  before <- readLines(file.path(root, "GWAS_information.csv"))
  expect_error(updateStatDB(ga, data.frame(b = 1)), "Missing required columns")
  expect_error(updateStatDB(ga, transform(stat, b = Inf, seb = 0, eaf = 1.1)),
    "Non-finite beta.*strictly positive.*effect-allele frequencies")
  expect_identical(readLines(file.path(root, "GWAS_information.csv")), before)
  expect_false(file.exists(file.path(root, "GWAS5.txt.gz")))
})

test_that("HbA1c-style ingestion is readable with catalog EAF and study N", {
  root <- tempfile("gact-hba1c-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  markers <- data.frame(rsids = "rs1", chr = 1L, pos = 10L,
    ea = "A", nea = "G", eaf = 0.2)
  data.table::fwrite(markers, file.path(root, "markers.txt.gz"))
  data.table::fwrite(data.frame(cpra = "1_10_A_G"), file.path(root, "cpra.txt.gz"))
  ga <- list(version = "hsa.0.0.1", dirs = c(marker = root, gstat = root), rsids = "rs1")
  raw <- data.frame(`#CHROM` = "chr1", POS = 10L, ALT = "G", REF = "A",
    Effect = 0.8, StdErr = 0.1, `P-value` = 0.04, check.names = FALSE)
  schema <- c(chr = "#CHROM", pos = "POS", ea = "ALT", nea = "REF",
    b = "Effect", seb = "StdErr", p = "P-value")
  expect_error(ingestStatDB(ga, raw, schema = schema, build = "hg38"), "differs from")
  expect_false(file.exists(file.path(root, "GWAS_information.csv")))

  result <- ingestStatDB(ga, raw, schema = schema, effect_scale = "beta",
    build = "GRCh37", n = 1000, type = "quantitative")
  study <- getGWAS(result$GAlist, "GWAS1")
  expect_identical(study$rsids, "rs1")
  expect_equal(study$b, -0.8)
  expect_equal(study$eaf, 0.2)
  expect_equal(study$n, 1000)

  log_stat <- transform(raw, `P-value` = NULL, Score = 0.3)
  schema[["p"]] <- "Score"
  result <- ingestStatDB(result$GAlist, log_stat, schema = schema, p_scale = "neglog10p",
    build = "GRCh37", n = 1000, type = "quantitative")
  expect_equal(getGWAS(result$GAlist, "GWAS2")$p, 10^-0.3)
})
