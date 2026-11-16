test_that("ingestion reconciles QC counts and is independent of input order", {
  root <- tempfile("gact-import-qc-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  catalog <- data.frame(chr = 1L,
    pos = c(10L, 20L, 30L, 31L, 40L, 50L, 50L, 60L, 60L, 70L, 70L, 80L, 80L, 90L),
    rsids = c("rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs7", "shared", "shared",
      "rs8", "rs9", "rs10", "rs11", "rs12"),
    ea = c("A", "C", "A", "C", "A", "A", "A", "A", "A", "G", "G", "G", "T", "A"),
    nea = c("G", "T", "T", "G", "C", "G", "C", "G", "C", "T", "T", "T", "G", "G"),
    eaf = c(rep(0.2, 13), 0.1))
  data.table::fwrite(catalog, file.path(root, "markers.txt.gz"))
  keys <- with(catalog, paste(chr, pos, ea, nea, sep = "_"))
  data.table::fwrite(data.frame(cpra = keys), file.path(root, "cpra.txt.gz"))
  ga <- list(version = "hsa.0.0.1", dirs = c(marker = root, gstat = root),
    rsids = catalog$rsids)
  row <- function(pos, ea, nea, b = 0.2, eaf = 0.2, p = 0.04, note = "") {
    data.frame(chr = 1L, pos = pos, ea = ea, nea = nea, b = b, eaf = eaf,
      seb = 0.1, p = p, z = b / 0.1, note = note)
  }
  raw <- rbind(
    row(10, "A", "G", note = "second"),
    row(10, "A", "G", note = "first"),
    row(10, "G", "A", b = -0.2, eaf = 0.8),
    row(20, "T", "C", eaf = 0.8),
    row(20, "G", "A"),
    row(30, "T", "A", eaf = 0.8),
    row(31, "G", "C", eaf = 0.8),
    row(40, "A", "C"),
    row(40, "C", "A", b = -0.3, eaf = 0.8),
    row(50, "A", "G"),
    row(50, "A", "C"),
    row(60, "A", "G"),
    row(60, "A", "C"),
    row(70, "G", "T"),
    row(80, "G", "T"),
    row(90, "G", "A", eaf = 0.1),
    row(100, "A", "G"),
    row(50, "A", "G,C"),
    row(10, "A", "G", p = NA_real_)
  )
  schema <- c(p = "p")
  first <- ingestStatDB(ga, raw, schema = schema, build = "GRCh37", n = 1000)
  qc <- first$qc
  expect_identical(qc$n_input, 19L)
  expect_identical(qc$n_final, 4L)
  expect_identical(qc$n_flipped, 1L)
  expect_identical(qc$n_dropped, 15L)
  expect_equal(qc$n_final + qc$n_dropped, qc$n_input)
  expect_equal(sum(unlist(qc$dropped_by_reason)), qc$n_dropped)
  expect_equal(unlist(qc$dropped_by_reason), c(missing_required = 1,
    unsplit_multiallelic = 1, unresolved_strand = 2, unmatched = 2,
    ambiguous_catalog_match = 2, shared_catalog_rsid = 2, conflicting_duplicate = 2,
    redundant_duplicate = 2, eaf_mismatch = 1))
  qc_path <- file.path(root, first$GAlist$study$qc_file[1L])
  expect_identical(basename(qc_path), "GWAS1.qc.json")
  expect_equal(jsonlite::read_json(qc_path), qc)
  expect_identical(read.csv2(file.path(root, "GWAS_information.csv"))$qc_file, "GWAS1.qc.json")

  publish <- gact:::.publish_gwas_import_file
  local({
    local_mocked_bindings(.publish_gwas_import_file = function(from, to, replace = FALSE) {
      if (basename(to) == "GWAS_information.csv") stop("Injected inventory failure")
      publish(from, to, replace)
    }, .package = "gact")
    expect_error(ingestStatDB(first$GAlist, raw, schema = schema), "Injected inventory failure")
  })
  expect_false(file.exists(file.path(root, "GWAS2.qc.json")))

  reversed <- raw[rev(seq_len(nrow(raw))), ]
  second <- ingestStatDB(first$GAlist, reversed, schema = schema, build = "GRCh37", n = 1000)
  expect_identical(second$qc, qc)
  stored <- function(id) data.table::fread(second$GAlist$studyfiles[id], data.table = FALSE)
  expect_identical(stored("GWAS1"), stored("GWAS2"))
  expect_identical(stored("GWAS1")$rsids, c("rs1", "rs2", "rs6", "rs7"))
  expect_identical(stored("GWAS1")$note[1L], "first")
  expect_equal(stored("GWAS1")$b, c(0.2, -0.2, 0.2, 0.2))

  expect_error(ingestStatDB(second$GAlist, raw, schema = schema, strand = "forward"),
    "Declare the GWAS build")
  third <- ingestStatDB(second$GAlist, raw, schema = schema, build = "GRCh37", strand = "forward",
    n = 1000)
  expect_equal(third$qc$n_final, 6)
  expect_equal(third$qc$n_flipped, 3)
  expect_equal(third$qc$dropped_by_reason$unresolved_strand, 0)
  expect_identical(getGWAS(third$GAlist, "GWAS3", what = "b")[c("rs3", "rs4"), 1],
    c(rs3 = -0.2, rs4 = -0.2))
  ga <- removeStatDB(third$GAlist, "GWAS3")
  expect_false(file.exists(file.path(root, "GWAS3.qc.json")))
  expect_error(ingestStatDB(ga, raw[6:7, ], schema = schema), "unresolved_strand=2")
  expect_false(file.exists(file.path(root, "GWAS3.txt.gz")))
})
