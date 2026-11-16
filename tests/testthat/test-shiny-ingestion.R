test_that("Shiny downloads the descriptor after successful ingestion", {
  root <- tempfile("gact-shiny-import-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  markers <- data.frame(rsids = "rs1", chr = 1L, pos = 10L, ea = "A", nea = "G", eaf = 0.2)
  data.table::fwrite(markers, file.path(root, "markers.txt.gz"))
  data.table::fwrite(data.frame(cpra = "1_10_A_G"), file.path(root, "cpra.txt.gz"))
  galist_path <- file.path(root, "GAlist.rds")
  saveRDS(list(version = "hsa.0.0.1", dirs = c(marker = root, gstat = root),
    rsids = markers$rsids), galist_path)
  original_md5 <- tools::md5sum(galist_path)
  upload <- file.path(root, "study.tsv")
  invalid_upload <- file.path(root, "invalid.tsv")
  raw <- transform(markers, b = 0.2, seb = 0.1, p = 0.04)
  data.table::fwrite(raw, upload, sep = "\t")
  data.table::fwrite(transform(raw, seb = 0), invalid_upload, sep = "\t")
  app_path <- system.file("shiny", "ingest_app", "app.R", package = "gact")
  app <- suppressPackageStartupMessages(source(app_path, local = TRUE)$value)

  shiny::testServer(app, {
    session$setInputs(db_file = list(name = "GAlist.rds", datapath = galist_path),
      files = list(name = "invalid.tsv", datapath = invalid_upload), file_select = "invalid.tsv")
    session$setInputs(map_chr = "chr", map_pos = "pos", map_ea = "ea", map_nea = "nea",
      map_b = "b", map_seb = "seb", map_p = "p", map_eaf = "eaf",
      effect_scale = "beta", p_scale = "p", source = "fixture", trait = "test",
      type = "quantitative", category = "test", efo = "unknown", gender = "both",
      ancestry = "EUR", build = "GRCh37", n = 1000, ncase = NA, ncontrol = NA,
      reference = "unknown", comments = "")
    session$setInputs(ingest = 1)
    expect_false(has_imported())
    expect_false(file.exists(file.path(root, "GWAS_information.csv")))

    session$setInputs(files = list(name = "study.tsv", datapath = upload), file_select = "study.tsv")
    session$setInputs(ingest = 2)
    expect_true(has_imported())
    expect_match(output$download_galist_ui$html, "Download updated GAlist")
    downloaded <- output$download_galist
    expect_identical(basename(downloaded), "GAlist-updated.rds")
    reopened <- readRDS(downloaded)
    expect_identical(reopened$studies$id, "GWAS1")
    expect_equal(getGWAS(reopened, "GWAS1")$b, 0.2)
    expect_identical(tools::md5sum(galist_path), original_md5)
  })
})
