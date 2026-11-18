test_that("IGSR preparation preserves allele and rsID rules and writes verified receipts", {
  skip_if_not_installed("gsim")
  skip_if_not("gsim_scan_vcf" %in% getNamespaceExports("gsim"))
  root <- tempfile("igsr-")
  on.exit(unlink(root, recursive = TRUE))
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
  testthat::local_mocked_bindings(curl_download = function(url, destfile, quiet, mode) {
    expect_false(quiet)
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
  messages <- capture.output(input <- downloadReferenceInputs("IGSR", cache, chromosomes = 1L),
                             type = "message")
  expect_length(grep("^gnomAD chr1 (download|checksum|scan and lookup write): [0-9.]+ s$",
                     messages), 3L)
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

  bundle <- packageGnomadLookups(cache, file.path(root, "bundle"), chromosomes = 1L)
  fresh <- file.path(root, "fresh")
  expect_error(installGnomadLookups(bundle$archive, fresh, checksum_file = fasta),
               "archive SHA-256 mismatch")
  expect_false(dir.exists(fresh))
  installed <- installGnomadLookups(bundle$archive, fresh)
  expect_identical(installGnomadLookups(bundle$archive, fresh), installed)
  members <- list.files(input$rsid_dir)
  expect_identical(unname(tools::md5sum(file.path(installed, members))),
                   unname(tools::md5sum(file.path(input$rsid_dir, members))))
  downloads <- character()
  reuse <- downloadReferenceInputs("IGSR", fresh, chromosomes = 1L)
  expect_length(downloads, 4L)
  expect_false(basename(gnomad) %in% downloads)
  output <- file.path(root, "prepared-from-bundle")
  reused <- prepareIGSRMarkers(1L, reuse$panel_dir, reuse$fasta_dir, reuse$rsid_dir, output)
  expect_identical(reused$sources, receipt$sources)
  expect_identical(data.table::fread(file.path(output, "chr1.tsv.gz")), markers)
  writeLines("changed", file.path(installed, basename(rsid)))
  expect_error(installGnomadLookups(bundle$archive, fresh), "destination differs")
  expect_identical(readLines(file.path(installed, basename(rsid))), "changed")
  expect_false(dir.exists(file.path(fresh, ".gact-igsr-download")))

  altered <- file.path(root, "altered")
  dir.create(altered)
  utils::untar(bundle$archive, exdir = altered, tar = "internal")
  writeLines("changed", file.path(altered, basename(rsid)))
  bad_archive <- file.path(root, "altered.tar.gz")
  local({
    previous <- setwd(altered)
    on.exit(setwd(previous))
    expect_equal(utils::tar(bad_archive, list.files(),
                            compression = "gzip", tar = "internal"), 0L)
  })
  writeLines(gact:::.gact_sha256_file(bad_archive), paste0(bad_archive, ".sha256"))
  expect_error(installGnomadLookups(bad_archive, file.path(root, "invalid")),
               "member checksum mismatch")
  expect_false(dir.exists(file.path(root, "invalid", "gnomad-rsids")))

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

test_that("parallel gnomAD workers reuse lookups and identify a failing chromosome", {
  root <- tempfile("gnomad-parallel-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  stage <- file.path(root, "stage")
  dir.create(stage)
  for (chr in 1:2) {
    paths <- file.path(root, gact:::.gact_gnomad_lookup_files(chr))
    con <- gzfile(paths[1L], "wt")
    writeLines(c("#CHROM\tPOS\tID\tREF\tALT", paste0("chr", chr, "\t1\trs", chr, "\tA\tG")), con)
    close(con)
    jsonlite::write_json(list(version = "4.1", chromosome = chr,
      url = paste0("https://example.org/chr", chr, ".vcf.bgz"),
      source_md5 = strrep("a", 32L), output_md5 = unname(tools::md5sum(paths[1L]))),
      paths[2L], auto_unbox = TRUE)
  }
  messages <- capture.output(gact:::.gact_prepare_gnomad_lookups(1:2, root, stage, 2L),
                             type = "message")
  expect_true("Preparing or verifying gnomAD lookups with 2 workers." %in% messages)
  expect_true(all(paste0("Reusing verified gnomAD rsID lookup for chr", 1:2) %in% messages))
  expect_length(list.files(stage, all.files = TRUE, no.. = TRUE), 0L)

  first <- file.path(root, gact:::.gact_gnomad_lookup_files(1L))
  hashes <- tools::md5sum(first)
  writeLines("changed", file.path(root, "gnomad_rsid_chr2.tsv.bgz"))
  expect_error(gact:::.gact_prepare_gnomad_lookups(1:2, root, stage, 2L),
               "failed for chr2: Invalid gnomAD lookup or provenance")
  expect_identical(tools::md5sum(first), hashes)
  expect_length(list.files(stage, all.files = TRUE, no.. = TRUE), 0L)
})
