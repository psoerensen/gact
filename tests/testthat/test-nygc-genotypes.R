test_that("NYGC merging preserves chromosome and sample order and checks input integrity", {
  root <- tempfile("nygc-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  parts <- file.path(root, "parts")
  dir.create(parts)
  groups <- c(afr = 661L, amr = 347L, eas = 504L, eur = 503L, sas = 489L)
  header <- as.raw(c(0x6c, 0x1b, 0x01))
  regions <- data.frame(chromosome = c(22L, 1L), start = 10L, end = 20L)

  write_receipt <- function(chr, files) {
    paths <- file.path(parts, files$files)
    jsonlite::write_json(list(build = "GRCh38", release = "20220422_3202_phased_SNV_INDEL_SV",
      chromosome = chr, region = c(10L, 20L),
      panels = data.frame(ancestry = toupper(names(groups)), markers = 1L),
      inputs = as.list(stats::setNames(rep(strrep("a", 32), 4L),
                                      c("vcf", "index", "sample_panel", "catalog"))),
      outputs = as.list(stats::setNames(unname(tools::md5sum(paths)), files$files))
    ), file.path(parts, files$receipt), auto_unbox = TRUE)
  }
  for (chr in regions$chromosome) {
    files <- gact:::.gact_nygc_chromosome_files(chr, c(10L, 20L))
    for (i in seq_along(groups)) {
      prefix <- file.path(parts, files$prefixes[i])
      n <- groups[[i]]
      writeBin(c(header, rep(as.raw(chr), ceiling(n / 4))), paste0(prefix, ".bed"))
      writeLines(paste0("chr", chr, "\trs", chr, "\t0\t10\tA\tC"), paste0(prefix, ".bim"))
      writeLines(paste(0, paste0("sample", seq_len(n)), 0, 0, 0, -9), paste0(prefix, ".fam"))
    }
    write_receipt(chr, files)
  }

  out <- file.path(root, "marker")
  summary <- mergeNYGCGenotypes(parts, out, regions)
  expect_equal(summary$samples, unname(groups))
  expect_equal(summary$markers, rep(2L, 5L))
  expect_length(list.files(out, pattern = "\\.(bed|bim|fam)$"), 15L)
  for (i in seq_along(groups)) {
    prefix <- file.path(out, summary$prefix[i])
    width <- ceiling(groups[[i]] / 4)
    expect_identical(readBin(paste0(prefix, ".bed"), "raw", n = 3L + 2L * width),
                     c(header, rep(as.raw(1), width), rep(as.raw(22), width)))
    expect_identical(sub("\t.*", "", readLines(paste0(prefix, ".bim"))), c("chr1", "chr22"))
  }
  expect_error(mergeNYGCGenotypes(parts, out, regions), "new directory")

  bad_out <- file.path(root, "failed")
  files <- gact:::.gact_nygc_chromosome_files(22L, c(10L, 20L))
  fam <- file.path(parts, paste0(files$prefixes[1L], ".fam"))
  writeLines(rev(readLines(fam)), fam)
  expect_error(mergeNYGCGenotypes(parts, bad_out, regions), "differ from")
  write_receipt(22L, files)
  expect_error(mergeNYGCGenotypes(parts, bad_out, regions), "Inconsistent NYGC FAM")
  writeLines(rev(readLines(fam)), fam)
  bed <- file.path(parts, paste0(files$prefixes[5L], ".bed"))
  width <- ceiling(groups[[5L]] / 4)
  writeBin(c(header[1:2], as.raw(0), rep(as.raw(22), width)), bed)
  write_receipt(22L, files)
  expect_error(mergeNYGCGenotypes(parts, bad_out, regions), "variant-major BED header")
  expect_false(file.exists(bad_out))
  expect_false(file.exists(paste0(bad_out, ".lock")))
})
