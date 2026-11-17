test_that("IGSR preparation preserves allele and rsID rules and writes verified receipts", {
  skip_if_not_installed("gsim")
  skip_if_not("gsim_scan_vcf" %in% getNamespaceExports("gsim"))
  root <- tempfile("igsr-")
  dirs <- file.path(root, c("panel", "fa", "rsid", "prepared"))
  for (dir in dirs) dir.create(dir, recursive = TRUE)
  panel_name <- "ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
  panel <- file.path(dirs[1L], panel_name)
  fasta <- file.path(dirs[2L], "chr1.fa.gz")
  rsid <- file.path(dirs[3L], "gnomad_rsid_chr1.tsv.bgz")
  write_gz <- function(lines, path) {
    connection <- gzfile(path, "wt")
    on.exit(close(connection))
    writeLines(lines, connection)
  }
  write_gz(c(
    "##fileformat=VCFv4.2", "##reference=GRCh38",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "1\t1\trs1\tA\tG\t.\tPASS\tAC=2;AN=100",
    "1\t2\t.\tA\tC\t.\tPASS\tAN=100;AC=12",
    "1\t3\trs3\tA\tG\t.\tPASS\tAC=15;AN=100",
    "1\t4\trs4\tA\tC\t.\tPASS\tAC=0;AN=100",
    "1\t5\trs5\tC\tT\t.\tPASS\tAC=2;AN=100",
    "1\t6\trs6\tA\tC,G\t.\tPASS\tAC=2,3;AN=100",
    "1\t7\trs7\tA\tC\t.\t.\tAC=2;AN=100",
    "1\t8\trs8\tA\tG\t.\tPASS\tAC=21;AN=100",
    "1\t8\t.\tA\tC\t.\tPASS\tAC=22;AN=100",
    "1\t9\t.\tA\tC\t.\tPASS\tAC=1;AN=100",
    "1\t10\trs10\tA\tC\t.\tPASS\tAC=2;AN=100"
  ), panel)
  write_gz(c(">chr1", "AAAAAAAAAA"), fasta)
  write_gz(c("#CHROM\tPOS\tID\tREF\tALT",
             "chr1\t1\trs111\tA\tT",
             "chr1\t2\trs2;rs2;other\tA\tC",
             "chr1\t3\trs33\tA\tG",
             "chr1\t8\trs8\tA\tG",
             "chr1\t8\trs88\tA\tC",
             "chr1\t8\trs88\tA\tC",
             "chr1\t9\trs9;rs99\tA\tC"), rsid)
  writeLines(paste0("./", panel_name, " 1 ", tools::md5sum(panel)),
             file.path(dirs[1L], "20190312_biallelic_SNV_and_INDEL_MANIFEST.txt"))
  writeLines(paste(tools::md5sum(fasta), basename(fasta)),
             file.path(dirs[2L], "md5sum.txt"))

  gnomad <- file.path(dirs[3L], "gnomad.genomes.v4.1.sites.chr1.vcf.bgz")
  lookup <- data.table::fread(rsid, header = TRUE, data.table = FALSE)
  write_gz(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
             do.call(paste, c(unname(as.list(lookup)),
                             list(".", "PASS", ".", sep = "\t")))), gnomad)
  md5 <- unname(tools::md5sum(gnomad))
  hash <- openssl::base64_encode(as.raw(strtoi(substring(md5, seq(1, 31, 2),
                                                        seq(2, 32, 2)), base = 16L)))
  downloads <- character()
  testthat::local_mocked_bindings(curl_download = function(url, destfile, mode) {
    name <- basename(url)
    downloads <<- c(downloads, name)
    source <- if (name == basename(gnomad)) dirs[3L] else if (
      grepl("hgdownload.soe.ucsc.edu", url, fixed = TRUE)) dirs[2L] else dirs[1L]
    stopifnot(file.copy(file.path(source, name), destfile))
  }, curl_fetch_memory = function(url) {
    stopifnot(grepl("release%2F4.1%2Fvcf%2Fgenomes", url, fixed = TRUE))
    list(status_code = 200L,
         content = charToRaw(jsonlite::toJSON(list(md5Hash = hash), auto_unbox = TRUE)))
  }, .package = "curl")
  cache <- file.path(root, "downloads")
  input <- downloadReferenceInputs("IGSR", cache, chromosomes = 1L)
  expect_length(downloads, 5L)
  expect_identical(downloadReferenceInputs("IGSR", cache, chromosomes = 1L), input)
  expect_length(downloads, 5L)
  receipt <- prepareIGSRMarkers(1L, input$panel_dir, input$fasta_dir, input$rsid_dir, dirs[4L])
  expect_identical(receipt$sources[[3L]]$upstream$source_md5, md5)
  expect_length(list.files(cache, pattern = "\\.vcf\\.bgz$", recursive = TRUE), 0L)
  table <- file.path(dirs[4L], "chr1.tsv.gz")
  markers <- data.table::fread(table)
  expect_identical(markers$rsids, c("rs1", "rs2", "rs8", "rs88", "rs10"))
  expect_equal(markers$af, c(.02, .12, .21, .22, .02))
  expect_equal(unname(unlist(receipt$counts)), c(11, 2, 1, 1, 2, 5))
  expect_identical(receipt$output$md5, unname(tools::md5sum(table)))
  expect_length(jsonlite::read_json(file.path(dirs[4L], "chr1.json"))$sources, 3L)
  expect_error(prepareIGSRMarkers(1L, dirs[1L], dirs[2L], dirs[3L], dirs[4L]),
               "already exists")
  write_gz(c(">chr1", "GGGGGGGGGG"), fasta)
  expect_error(prepareIGSRMarkers(1L, dirs[1L], dirs[2L], dirs[3L], dirs[4L]),
               "checksum mismatch")
  expect_error(downloadReferenceInputs("IGSR", file.path(root, "bad"), chromosomes = 1L),
               "Release checksum mismatch for chr1")
  expect_false(file.exists(file.path(root, "bad", "ucsc-hg38", "chr1.fa.gz")))
  expect_error(prepareIGSRMarkers(1L, dirs[1L], dirs[2L], file.path(root, "missing"), dirs[4L]),
               "gnomad_rsid_chr1\\.tsv\\.bgz")
})
