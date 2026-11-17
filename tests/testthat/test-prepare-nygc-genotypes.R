test_that("NYGC preparation preserves catalog identity and ancestry-specific GT filtering", {
  skip_if_not_installed("gsim")
  skip_if_not("gsim_scan_vcf" %in% getNamespaceExports("gsim"))
  root <- tempfile("nygc-native-")
  dir.create(root)
  out <- file.path(root, "out")
  dir.create(out)
  vcf <- file.path(root, "1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
  index <- paste0(vcf, ".tbi")
  samples <- sprintf("id%04d", seq_len(2504L))
  ancestry <- rep(c("AFR", "AMR", "EAS", "EUR", "SAS"), length.out = length(samples))
  panel <- data.frame(sample = samples, super_pop = ancestry, gender = "female")
  data.table::fwrite(panel, file.path(root, "integrated_call_samples_v3.20130502.ALL.panel"),
                     sep = "\t")
  calls <- rep("0|1", length(samples))
  calls_missing <- calls
  calls_missing[which(ancestry == "EUR")[1L]] <- ".|."
  record <- function(pos, id, alt, filter, cm, gt = calls) {
    paste(c("chr1", pos, id, "A", alt, ".", filter, paste0("CM=", cm), "GT", gt),
          collapse = "\t")
  }
  connection <- gzfile(vcf, "wt")
  writeLines(c(
    "##fileformat=VCFv4.2",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
            samples), collapse = "\t"),
    record(10, "v10", "C", "PASS", "0.1", rep("1|1", length(samples))),
    record(20, "v20", "G", ".", "0.2", rep("0|0", length(samples))),
    record(30, "v30", "T", "PASS", "0.3", calls_missing),
    record(40, "v40", "T", "FAIL", "0.4"),
    record(50, "v50", "C", "PASS", "0.5"),
    record(55, "v55", "C,G", "PASS", "0.55"),
    record(60, "v60", "G", "PASS", "0.6")
  ), connection)
  close(connection)
  writeLines("Synthetic release checksum fixture", index)
  writeLines(paste(basename(c(vcf, index)), tools::md5sum(c(vcf, index))),
             file.path(root, "20220804_manifest.txt"))
  catalog <- file.path(root, "markers.txt.gz")
  positions <- c(10, 20, 30, 40, 50, 55, 60)
  data.table::fwrite(data.frame(rsids = paste0("rs", positions), chr = 1L,
    pos = positions, ea = c("C", "G", "T", "T", "G", "C,G", "G"), nea = "A"),
    catalog, compress = "gzip")

  downloads <- character()
  testthat::local_mocked_bindings(curl_download = function(url, destfile, quiet, mode) {
    expect_false(quiet)
    downloads <<- c(downloads, basename(url))
    stopifnot(file.copy(file.path(root, basename(url)), destfile))
  }, .package = "curl")
  cache <- file.path(root, "downloads")
  input <- downloadReferenceInputs("NYGC", cache, chromosomes = 1L)
  expect_length(downloads, 4L)
  expect_identical(downloadReferenceInputs("NYGC", cache, chromosomes = 1L), input)
  expect_length(downloads, 4L)
  result <- prepareNYGCGenotypes(input$input_dir, catalog, 1L, out, region = c(10, 55))
  expect_equal(result$markers, c(3, 3, 3, 2, 3))
  expect_equal(result$samples, unname(as.integer(table(ancestry))))
  read_bim <- function(pop) data.table::fread(
    file.path(out, paste0("g1000_", pop, ".chr1.10-55.bim")), header = FALSE)
  eur <- read_bim("eur")
  afr <- read_bim("afr")
  expect_identical(eur$V2, c("rs10", "rs20"))
  expect_identical(afr$V2, c("rs10", "rs20", "rs30"))
  expect_equal(afr$V3, c(.1, .2, .3))
  expect_identical(afr$V5, c("C", "G", "T"))
  expect_identical(afr$V6, rep("A", 3L))
  bed_path <- file.path(out, "g1000_afr.chr1.10-55.bed")
  bytes <- readBin(bed_path, "raw", n = file.size(bed_path))
  stride <- ceiling(result$samples[1L] / 4)
  expect_identical(bytes[4 + (0:2) * stride], as.raw(c(0x00, 0xff, 0xaa)))
  expect_equal(jsonlite::read_json(file.path(out, "nygc-chr1-10-55.json"))$matched, 3L)
  expect_error(prepareNYGCGenotypes(root, catalog, 2L, out),
               "Missing NYGC preparation inputs for chr2")
  writeLines("corrupt", file.path(input$input_dir, basename(index)))
  expect_error(downloadReferenceInputs("NYGC", cache, chromosomes = 1L),
               "Release checksum mismatch")
  expect_length(downloads, 4L)
})
