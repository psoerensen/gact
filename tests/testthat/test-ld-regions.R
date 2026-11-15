test_that("assocMatrix exports ancestry-specific half-open GRCh37 regions", {
  root <- tempfile("assoc-regions-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  marker_dir <- file.path(root, "marker")
  dir.create(marker_dir)
  ga <- downloadDB(list(dirs = c(marker = marker_dir)), what = "ldetect")
  published <- ga$ldRegions$GRCh37
  files <- vapply(published, function(x) x$file, character(1L))
  expect_identical(unname(dirname(files)), rep(normalizePath(marker_dir, winslash = "/"), 2L))
  expect_identical(unname(tools::md5sum(files)), unname(vapply(published, function(x) x$md5, "")))
  expect_identical(downloadDB(ga, what = "ldetect"), ga)
  intervals <- lapply(files, utils::read.table, header = TRUE)
  eur_stop <- intervals$EUR$stop[1L]
  eas_stop <- intervals$EAS$stop[1L]
  expect_identical(names(published), c("EUR", "EAS"))
  expect_identical(vapply(intervals, nrow, integer(1)), c(EUR = 1703L, EAS = 1445L))
  expect_identical(c(eur_stop, eas_stop), c(1892607L, 1976558L))
  expect_error(downloadDB(list(assembly = "GRCh38"), what = "ldetect"), "differs")
  created <- suppressMessages(gact:::createDB("hsa.0.0.1", root))
  expect_identical(unname(tools::md5sum(vapply(created$ldRegions$GRCh37, function(x) x$file, ""))),
                   unname(tools::md5sum(files)))

  markers <- data.frame(rsids = paste0("rs", 1:8), chr = c(rep("1", 7), "X"),
    pos = c(10583L, eur_stop - 1L, eur_stop, eur_stop + 1L, eas_stop, eas_stop + 1L, eas_stop + 2L, 20L),
    ea = "A", nea = "C")
  data.table::fwrite(markers, file.path(marker_dir, "markers.txt.gz"))
  eur <- markers[1:7, ]
  eur$z <- -(1:7)
  eas <- markers[1:8, ]
  eas$z <- (1:8) * 10
  sas <- eas
  for (panel in c("eur", "eas", "sas"))
    data.table::fwrite(get(panel), file.path(root, paste0(panel, ".tsv")), sep = "\t")
  ga$rsids <- markers$rsids[c(5, 3, 2, 1, 8, 4, 7)]
  ga$studies <- data.frame(id = c("eas", "sas", "eur"), ancestry = c("EAS", "SAS", "EUR"))
  ga$studyfiles <- setNames(file.path(root, paste0(ga$studies$id, ".tsv")), ga$studies$id)
  make_panel <- function(name, idx) {
    prefix <- file.path(marker_dir, paste0("g1000_", name, ".subset.aligned"))
    bed <- paste0(prefix, ".bed")
    bim <- paste0(prefix, ".bim")
    fam <- paste0(prefix, ".fam")
    writeBin(as.raw(c(0x6c, 0x1b, 0x01, ifelse(idx == 2L, 0x38, 0xcb))), bed)
    writeLines(paste("f", paste0("i", 1:4), 0, 0, 0, -9), fam)
    writeLines(with(markers[idx, ], paste(chr, rsids, 0, pos, ea, nea)), bim)
    gl <- gbase::gprep(study = name, bedfiles = bed, bimfiles = bim, famfiles = fam)
    gbase::gprep(gl, task = "sparseld", max_distance_bp = 0L,
      max_distance_variants = 5L, progress = FALSE)
  }
  panels <- list(EUR = make_panel("eur", c(1:3, 5:6, 8)),
    EAS = make_panel("eas", c(2:3, 5:6, 8)), SAS = make_panel("sas", c(1:3, 5)))
  bad <- panels
  bad$EUR$assembly <- "GRCh38"
  expect_error(assocMatrix(ga, file.path(root, "bad"), bad), "assembly differs")
  expect_false(dir.exists(file.path(root, "bad")))
  out <- gact::writeMarkerStat(ga)
  expect_identical(out$dim, c(7L, 3L))
  expect_identical(dirname(out$z), file.path(normalizePath(root, winslash = "/"), "gwas_matrix"))
  ld <- read.delim(out$ld)
  expect_identical(ld$ancestry, c("EUR", "EAS", "SAS"))
  expect_identical(ld$ld_markers, vapply(unname(panels), function(panel) {
    normalizePath(panel$ldfiles[["markers"]], winslash = "/")
  }, character(1L)))
  index <- read.delim(out$ld_regions)
  expect_identical(index$ancestry, c("EUR", "EAS"))
  expect_identical(index$assembly, rep("GRCh37", 2))
  expect_null(assocMatrix(ga, file.path(root, "sas-matrix"), panels["SAS"])$ld_regions)
  pairs <- function(i) {
    path <- file.path(dirname(out$z), index$row_pairs[i])
    readBin(path, "integer", n = file.size(path) / 4, size = 4, endian = "little")
  }
  expect_identical(pairs(1), c(3L, 0L, 2L, 1L, 1L, 2L, 0L, 3L))
  expect_identical(pairs(2), c(2L, 0L, 1L, 1L, 0L, 2L))
  legacy <- ga
  legacy$ldRegions$GRCh37 <- lapply(intervals, function(x) list(intervals = x))
  old <- writeMarkerStat(legacy, out_dir = file.path(root, "legacy-regions"))
  expect_identical(unname(tools::md5sum(file.path(dirname(old$z), index$row_pairs))),
                   unname(tools::md5sum(file.path(dirname(out$z), index$row_pairs))))
  for (i in 1:2) {
    regions <- read.delim(file.path(dirname(out$z), index$regions[i]))
    expect_identical(regions$region, seq_len(nrow(regions)) - 1L)
    expect_identical(regions$count[1:3], if (i == 1) c(2L, 2L, 0L) else c(2L, 1L, 0L))
    expect_identical(regions$offset[1:3], if (i == 1) c(0L, 2L, 4L) else c(0L, 2L, 3L))
  }

  drone <- Sys.which("gminer_impute")
  if (!nzchar(drone)) {
    if (identical(Sys.getenv("GACT_TEST_GMINER"), "1")) stop("gminer_impute not found")
    return(invisible(NULL))
  }
  run <- function(assembly, ancestry, id, output, dense = FALSE) {
    args <- c("--assoc-region", shQuote(dirname(out$z)), assembly, ancestry, id, shQuote(output))
    if (dense) args <- c(args, "--dense-ld")
    suppressWarnings(system2(drone, args, stdout = TRUE, stderr = TRUE))
  }
  for (ancestry in c("EUR", "EAS")) {
    output <- file.path(root, paste0("slice-", ancestry))
    result <- run("GRCh37", ancestry, 0, output, dense = TRUE)
    expect_null(attr(result, "status"), info = paste(result, collapse = "\n"))
    is_eur <- ancestry == "EUR"
    rows <- read.delim(file.path(output, "rows.tsv"))
    expect_identical(rows$z_row, if (is_eur) c(3L, 2L) else c(2L, 1L))
    expect_identical(rows$ld_row, c(0L, 1L))
    studies <- read.delim(file.path(output, "studies.tsv"))
    expect_identical(studies$study_id, tolower(ancestry))
    expect_identical(studies$z_column, if (is_eur) 2L else 0L)
    expect_identical(readBin(file.path(output, "z.bin"), "double", n = 2, endian = "little"),
      if (is_eur) c(-1, -2) else c(20, 30))
    expect_equal(readBin(file.path(output, "ld.bin"), "double", n = 4, endian = "little"),
      c(1, -1, -1, 1))
  }
  second <- file.path(root, "second")
  expect_null(attr(run("GRCh37", "EUR", 1, second), "status"))
  expect_identical(read.delim(file.path(second, "rows.tsv"))$ld_row, c(2L, 3L))
  expect_identical(readBin(file.path(second, "z.bin"), "double", n = 2, endian = "little"), c(-3, -5))
  expect_false(file.exists(file.path(second, "ld.bin")))
  empty <- file.path(root, "empty")
  expect_null(attr(run("GRCh37", "EUR", 2, empty), "status"))
  expect_equal(file.size(file.path(empty, "z.bin")), 0)
  unmapped <- run("GRCh37", "SAS", 0, file.path(root, "no-regions"))
  expect_identical(attr(unmapped, "status"), 1L)
  expect_match(paste(unmapped, collapse = "\n"), "No assocMatrix regions for GRCh37 SAS")
  wrong <- run("GRCh38", "EUR", 0, file.path(root, "wrong-build"))
  expect_identical(attr(wrong, "status"), 1L)
  expect_match(paste(wrong, collapse = "\n"), "No assocMatrix regions for GRCh38 EUR")
})
