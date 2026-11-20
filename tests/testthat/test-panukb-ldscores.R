test_that("Pan-UKB score matching uses unique rsIDs and release positions", {
  root <- tempfile()
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  path <- file.path(root, "input.tsv")
  catalog <- file.path(root, "markers.tsv")
  data.table::fwrite(data.frame(SNP = c("rs1", "rs2", "rs2", "rs3", "rs4", "rs5", "bad"),
                                L2 = c(5, 8, 9, 2, 6, 7, 3)), path, sep = "\t")
  data.table::fwrite(data.frame(rsids = c("rs1", "rs3", "rs3", "rs5"), chr = c(2, 1, 1, 1),
                                pos = c(30, 20, 22, 55)), catalog, sep = "\t")
  result <- .gact_panukb_table(path, catalog)
  expect_identical(result$rows$rsids, c("rs5", "rs1"))
  expect_equal(result$rows$ldscores, c(7, 5))
  expect_equal(result$rows$pos, c(55, 30))
  expect_equal(unlist(result$counts), c(input = 7, invalid = 1, ambiguous_source = 2,
               outside_backbone = 1, ambiguous_backbone = 1, matched = 2))
  expect_error(.gact_panukb_panel("SAS"), "Pan-UKB ancestry")
  expect_equal(.gact_panukb_panel("EUR")$marker_count, 6805960)
  expect_identical(.gact_panukb_panel("PanUKBB-EUR"), .gact_panukb_panel("EUR"))
  expect_error(getLDscores(list(dirs = c(marker = root)), ancestry = "PanUKBB-EUR",
                           version = "1000G"), "conflicts")
})
