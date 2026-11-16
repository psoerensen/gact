test_that("writeMarkerStat defaults to aligned observed Z columns", {
  root <- tempfile("z-matrix-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  markers <- data.frame(
    rsids = c("003", "rs2", "rs3", "rs4"),
    chr = c("1", "1", "2", "2"), pos = c(10, 20, 30, 40),
    ea = c("A", "C", "T", "G"), nea = c("G", "T", "C", "A")
  )
  marker_dir <- file.path(root, "marker")
  dir.create(marker_dir)
  gstat_dir <- file.path(root, "gstat")
  dir.create(gstat_dir)
  data.table::fwrite(markers, file.path(marker_dir, "markers.txt.gz"), sep = "\t")
  first <- file.path(gstat_dir, "one.txt.gz")
  second <- file.path(gstat_dir, "two.csv")
  data.table::fwrite(data.frame(rsids = c("rs3", "other", "003"), chr = c(2, 2, 1),
    pos = c(30, 55, 10), ea = c("C", "G", "A"), nea = c("T", "A", "G"),
    b = c(2, 7, 2), seb = c(1, 1, 0)), first, sep = "\t")
  writeLines(c("marker,chr,pos,ea,nea,z,b,seb,note", 'rs2,chr1,20,c,t,-8,4,0.5,"note, here"',
    "rs3,2,30,T,C,Inf,4,2,NA", "003,1,10,A,G,NA,1,-1,NA",
    "rs4,2,40,A,G,3,NA,NA,NA"), second)
  ga <- list(
    dirs = c(marker = marker_dir, gstat = gstat_dir),
    rsids = c("rs3", "003", "rs2", "rs4"),
    studies = data.frame(id = c("study2", "study1"), ancestry = c("EAS", "EUR")),
    studyfiles = c(study1 = first, study2 = second)
  )

  expect_message(
    expect_output(out <- gact::writeMarkerStat(ga), "2/2 study1.*1 observed.*3 missing"),
    "No matching LD for study2 \\(EAS\\), study1 \\(EUR\\)"
  )
  expect_identical(dirname(out$z), file.path(normalizePath(gstat_dir, winslash = "/"), "matrix"))
  expect_identical(out$dim, c(4L, 2L))
  expect_equal(file.size(out$z), 64)
  expect_identical(readBin(out$z, "double", n = 8, size = 8, endian = "little"),
    c(2, NaN, 8, -3, -2, NaN, NaN, NaN))
  catalog <- read.delim(out$markers, colClasses = c(rsids = "character"))
  expect_identical(catalog$rsids, ga$rsids)
  expect_equal(catalog$pos, c(30, 10, 20, 40))
  metadata <- read.delim(out$studies)
  expect_identical(metadata$study_id, c("study2", "study1"))
  expect_identical(metadata$ancestry, c("EAS", "EUR"))
  expect_identical(metadata$observed, c(3L, 1L))
  expect_identical(metadata$missing, c(1L, 3L))
  expect_null(out$ld)
  expect_error(writeMarkerStat(ga, dirname(out$z)), "already exists")
  expect_error(writeMarkerStat(ga, file.path(root, "beta"), stat = "b"), 'stat must be "z"')
  expect_false(dir.exists(file.path(root, "beta")))

  custom_descriptor <- ga
  custom_descriptor$dirs <- c(marker = marker_dir)
  expect_error(writeMarkerStat(custom_descriptor), "gstat directory or out_dir")
  expect_silent(explicit <- writeMarkerStat(custom_descriptor,
    out_dir = file.path(root, "explicit"), progress = FALSE))
  expect_identical(unname(tools::md5sum(c(explicit$z, explicit$markers, explicit$studies))),
    unname(tools::md5sum(c(out$z, out$markers, out$studies))))
  expect_error(writeMarkerStat(ga, studyID = "study1"), "unused argument")
  insufficient <- ga
  rejected <- file.path(root, "insufficient")
  for (count in 0:1) {
    insufficient$studies <- ga$studies[seq_len(count), , drop = FALSE]
    expect_error(writeMarkerStat(insufficient, out_dir = rejected), "at least two studies")
    expect_false(dir.exists(rejected))
  }

  ga$rsids <- markers$rsids
  out_ordered <- writeMarkerStat(ga, file.path(root, "ordered"))
  expect_identical(read.delim(out_ordered$markers,
    colClasses = c(rsids = "character"))$rsids, markers$rsids)

  cat("rs2,1,20,C,T,1,1,1,NA\n", file = second, append = TRUE)
  bad <- file.path(root, "failed")
  expect_error(writeMarkerStat(ga, bad), "Repeated marker rs2")
  expect_false(dir.exists(bad))
  expect_length(list.files(root, pattern = "^\\.gact-z-", all.files = TRUE), 0L)
})

test_that("one Z allele matches all ancestry panels with different marker coverage", {
  root <- tempfile("z-alignment-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  markers <- data.frame(rsids = paste0("rs", 1:5), chr = "1",
    pos = (1:5) * 100L, ea = c("A", "C", "T", "G", "A"), nea = c("G", "T", "C", "A", "C"))
  data.table::fwrite(markers, file.path(root, "markers.txt.gz"))
  stats <- markers[1:4, ]
  stats$ea[3] <- "C"
  stats$nea[3] <- "T"
  stats$b <- c(1, 2, -3, 4)
  stats$seb <- 0.5
  stats$p <- 0.1
  studies <- list(id = c("GWAS6", "GWAS7", "GWAS8"),
    ancestry = c("European", "EAS", "South Asian"), neff = rep(1000, 3))
  study_file <- file.path(root, "stats.csv")
  data.table::fwrite(stats, study_file)
  ga <- list(dirs = c(marker = root), rsids = markers$rsids[c(3, 1, 4, 2)],
    study = studies, studyfiles = setNames(rep(study_file, 3), studies$id))

  make_panel <- function(name, rows, swapped = character()) {
    panel <- markers[rows, ]
    prefix <- file.path(root, name)
    bed <- paste0(prefix, ".bed")
    bim <- paste0(prefix, ".bim")
    fam <- paste0(prefix, ".fam")
    writeBin(as.raw(c(0x6c, 0x1b, 0x01,
      ifelse(panel$rsids %in% swapped, 0x38, 0xcb))), bed)
    writeLines(paste("f", paste0("i", 1:4), 0, 0, 0, -9), fam)
    a1 <- ifelse(panel$rsids %in% swapped, panel$nea, panel$ea)
    a2 <- ifelse(panel$rsids %in% swapped, panel$ea, panel$nea)
    writeLines(paste(panel$chr, panel$rsids, 0, panel$pos, a1, a2), bim)
    gl <- gbase::gprep(study = name, bedfiles = bed, bimfiles = bim, famfiles = fam)
    gbase::gprep(gl, task = "sparseld", max_distance_variants = 3L, progress = FALSE)
  }
  panels <- list(EUR = make_panel("g1000_eur.subset.aligned", c(1, 2, 3)),
    EAS = make_panel("g1000_eas.subset.aligned", c(2, 3, 4, 5)),
    SAS = make_panel("g1000_sas.subset.aligned", c(1, 3, 4)))
  registered <- gact::addLDPanels(ga, panels)
  expect_identical(names(registered$ldPanels), names(panels))
  expect_identical(names(registered$ldPanels$EUR), c("markers", "assembly"))
  path <- file.path(root, "GAlist.rds")
  saveRDS(registered, path)
  out <- writeMarkerStat(readRDS(path), out_dir = file.path(root, "aligned"))
  z <- matrix(readBin(out$z, "double", n = 12, size = 8, endian = "little"), nrow = 4)
  catalog <- read.delim(out$ld)
  expect_identical(catalog$ancestry, names(panels))
  expect_identical(catalog$ld_rows, c(3L, 4L, 3L))
  expect_identical(catalog$matched, rep(3L, 3))
  expect_equal(unname(file.size(out$ld_maps)), 4 * catalog$ld_rows)
  expect_identical(read.delim(out$studies)$ld_map, basename(unname(out$ld_maps)))
  expect_identical(z, matrix(rep(c(6, 2, 8, 4), 3), nrow = 4))
  for (j in seq_along(panels)) {
    gl <- panels[[j]]
    stat <- getMarkerStat(GAlist = ga, studyID = studies$id[j])
    stat$chr <- as.character(stat$chr)
    # Compare allele alignment through the R reference workflow.
    aligned <- gbase::checkStat(Glist = gl, stat = stat, excludeMAF = NULL,
      excludeMAFDIFF = NULL, excludeMISS = NULL, excludeHWE = NULL)
    map <- readBin(out$ld_maps[[j]], "integer", n = catalog$ld_rows[j],
      size = 4L, endian = "little")
    expect_identical(ga$rsids[map], aligned$rsids)
    expect_identical(unlist(gl$rsidsLD, use.names = FALSE)[map != 0L], aligned$rsids)
    expect_identical(z[map, j], aligned$b / aligned$seb)
  }
  expect_identical(readBin(out$ld_maps[["EAS"]], "raw", n = 16L),
    as.raw(c(4, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0)))

  custom <- make_panel("custom-ld", c(2, 3))
  overrides <- addLDPanels(ga, list(EAS = custom))
  alternate <- writeMarkerStat(overrides, out_dir = file.path(root, "custom"))
  alternate_ld <- read.delim(alternate$ld)
  expect_setequal(alternate_ld$ancestry, names(panels))
  expect_identical(alternate_ld$ld_markers[alternate_ld$ancestry == "EAS"],
    normalizePath(custom$ldfiles[["markers"]], winslash = "/"))
  expect_identical(alternate_ld$ld_rows[alternate_ld$ancestry == "EAS"], 2L)
  expect_identical(alternate_ld$ld_markers[alternate_ld$ancestry == "EUR"],
    catalog$ld_markers[catalog$ancestry == "EUR"])

  bad <- file.path(root, "bad")
  reversed <- panels
  reversed$EAS <- make_panel("reversed-eas", c(2, 3, 4, 5), "rs3")
  # The two files represent the same genotypes. Recoding rs3 swaps the
  # homozygote codes in BED and changes the sign of its incident LD entries.
  read_ld <- function(gl) {
    list(columns = readBin(gl$ldfiles[["col_idx"]], "integer", n = 6L,
      size = 4L, endian = "little"),
    values = readBin(gl$ldfiles[["values"]], "numeric", n = 6L,
      size = 4L, endian = "little"))
  }
  before <- read_ld(reversed$EAS)
  after <- read_ld(panels$EAS)
  expect_identical(before$columns, c(1L, 2L, 3L, 2L, 3L, 3L))
  expect_identical(after$columns, before$columns)
  expect_equal(before$values, after$values * c(-1, 1, 1, -1, -1, 1))
  expect_true(all(after$values > 0))
  expect_error(writeMarkerStat(addLDPanels(ga, reversed), bad), "Reversed LD dosage allele for rs3")
  expect_false(dir.exists(bad))

  path <- panels$EUR$ldfiles[["markers"]]
  ld <- read.delim(path)
  ld$allele1[ld$marker_id == "rs2"] <- "A"
  data.table::fwrite(ld, path, sep = "\t")
  expect_error(writeMarkerStat(ga, bad), "Allele mismatch for rs2")
  expect_false(dir.exists(bad))
  ld$position_bp[ld$marker_id == "rs2"] <- 205L
  data.table::fwrite(ld, path, sep = "\t")
  expect_error(writeMarkerStat(ga, bad), "Coordinate mismatch for rs2")

  expect_true(file.remove(panels$EAS$ldfiles[["values"]]))
  expect_error(writeMarkerStat(ga, out_dir = bad),
    "Incomplete sparse LD for EAS.*values.f32.bin")
  expect_false(dir.exists(bad))
})
