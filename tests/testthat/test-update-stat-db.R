test_that("successive studies persist complete checksum and recipe metadata", {
  root <- tempfile("study-metadata-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  local_mocked_bindings(qcStatDB = function(GAlist, stat, ...) stat,
                        .package = "gact")
  ga <- list(dirs = c(gstat = root), study = list(id = character()))
  stat <- data.frame(eaf = 0.2, b = 0.1, seb = 0.05, p = 0.04)
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
  expect_identical(readLines(file.path(root, "GWAS_information.csv")), before)
})
