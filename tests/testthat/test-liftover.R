test_that("Pan-UKB liftover verifies unique positions and keeps the effect on its allele", {
  root <- tempfile("liftover-grch37-grch38-")
  ga <- createDB("hsa.GRCh38.0.0.1", root)
  on.exit(unlink(root, recursive = TRUE))
  chain <- file.path(root, "example.chain.gz")
  con <- gzfile(chain, "wt")
  writeLines(c(
    "chain 100 chr1 100 + 0 7 chr1 100 + 10 18 1", "3 1 2", "3", "",
    "chain 90 chr1 100 + 20 23 chr1 100 - 10 13 2", "3", "",
    "chain 80 chr1 100 + 1 2 chr1 100 + 41 42 3", "1", "",
    "chain 70 chr1 100 + 30 31 chr2 100 + 10 11 4", "1", "",
    "chain 60 chr1 100 + 41 42 chr1 100 + 16 17 5", "1"
  ), con)
  close(con)
  catalog <- file.path(ga$dirs[["marker"]], "markers.txt.gz")
  data.table::fwrite(data.frame(rsids = paste0("rs", 1:5), chr = 1L,
    pos = c(11L, 13L, 17L, 18L, 90L), nea = c("A", "A", "A", "A", "T"),
    ea = "C"), catalog, sep = "\t")
  jsonlite::write_json(list(build = "GRCh38", effect_allele = "ALT",
    checksums = list("markers.txt.gz" = unname(tools::md5sum(catalog)))),
    file.path(dirname(catalog), "marker_manifest.json"), auto_unbox = TRUE)
  stat <- data.frame(chr = c(rep("chr1", 12L), "chrX"),
    pos = c(21L, 1:7, 7L, 31L, 42L, 51L, 1L),
    ref = c("A", "A", "A", "G", rep("A", 9L)),
    alt = c("G", rep("C", 10L), "AA", "C"), beta_EUR = seq_len(13L) / 10)
  variants <- file.path(root, "variants.tsv.gz")
  provider <- stat
  names(provider)[1L] <- "chrom"
  data.table::fwrite(provider, variants, sep = "\t")
  sources <- list(chain = list(url = "https://example.org/ucsc.chain.gz",
                               md5 = unname(tools::md5sum(chain))),
                  variant_manifest = list(url = "https://example.org/variants.tsv.gz"))
  local_mocked_bindings(.gact_liftover_sources = function() sources, .package = "gact")
  requests <- character()
  local_mocked_bindings(curl_download = function(url, destfile, ...) {
    requests <<- c(requests, url)
    input <- if (url == sources$chain$url) chain else variants
    stopifnot(file.copy(input, destfile))
  }, .package = "curl")
  chain_md5 <- sources$chain$md5
  sources$chain$md5 <- strrep("0", 32)
  expect_error(buildLiftoverMap(ga), "input checksum mismatch")
  expect_identical(requests, sources$chain$url)
  sources$chain$md5 <- chain_md5
  requests <- character()
  fwrite <- data.table::fwrite
  local_mocked_bindings(fwrite = function(x, file, ...) {
    if (grepl("\\.(gz|bgz)$", file, ignore.case = TRUE))
      stop("Synthetic data.table build has no gzip support")
    fwrite(x, file, ...)
  }, .package = "data.table")
  local({
    local_mocked_bindings(gzip = function(...) stop("Simulated compression failure"),
                          .package = "R.utils")
    expect_error(buildLiftoverMap(ga), "Simulated compression failure")
  })
  expect_false(dir.exists(gact:::.gact_liftover_dir(ga)))
  expect_length(list.files(ga$dirs[["marker"]], pattern = "^\\.gact-liftover-"), 0L)
  map <- buildLiftoverMap(ga)
  expect_identical(dirname(map), unname(ga$dirs[["marker"]]))
  expect_true(file.exists(file.path(map, "index.bin")))
  expect_identical(buildLiftoverMap(ga), map)
  expect_identical(requests, c(sources$chain$url, sources$variant_manifest$url))
  expect_identical(jsonlite::read_json(file.path(map, "manifest.json"))$inputs$variant_manifest$url,
                   sources$variant_manifest$url)
  expect_true(file.exists(file.path(ga$dirs[["reference"]], basename(map), "ucsc.chain.gz")))
  expect_identical(gact:::.gact_cache_liftover_source(ga, sources$chain),
                   file.path(ga$dirs[["reference"]], basename(map), "ucsc.chain.gz"))
  expect_length(requests, 2L)
  counts <- jsonlite::read_json(file.path(map, "manifest.json"))
  expect_equal(counts$non_autosomal, 1)
  counts <- counts$chromosomes[[1L]]$counts
  expect_equal(unlist(counts), c(input = 12, unsupported = 1, duplicate_source = 2,
    unmapped = 1, ambiguous = 1, different_chromosome = 1, outside_backbone = 1,
    allele_mismatch = 1, duplicate_target = 2, retained = 2))
  lifted <- liftGWASStats(stat, map)
  expect_identical(lifted$pos, c(90L, 11L))
  expect_identical(lifted$ref, c("T", "A"))
  expect_identical(lifted$alt, c("C", "C"))
  expect_identical(lifted$marker, c("rs5", "rs1"))
  expect_identical(lifted$beta_EUR, stat$beta_EUR[1:2])
  expect_identical(attr(lifted, "liftover_qc"),
                   list(n_input = 13L, n_lifted = 2L, n_dropped = 11L))
  expect_equal(nrow(liftGWASStats(stat[13L, , drop = FALSE], map)), 0)
  expect_identical(verifyGWASBuild(stat, map, min_informative = 2L)$build, "GRCh37")
  expect_identical(verifyGWASBuild(lifted, map, min_informative = 2L)$build, "GRCh38")
  expect_identical(verifyGWASBuild(rbind(stat[1L, ], lifted[1L, names(stat)]),
                                   map, min_informative = 1L)$build, NA_character_)
  expect_identical(verifyGWASBuild(stat[rep(1L, 20L), ], map)$grch37_only, 1L)
  expect_identical(verifyGWASBuild(stat, map)$build, NA_character_)
  expect_error(buildLiftoverMap(ga, variants, chain), "already exists")
  malformed <- file.path(root, "malformed.tsv")
  data.table::fwrite(data.frame(chr = 1L, chrom = 1L, pos = 1L, ref = "A", alt = "C"),
                     malformed, sep = "\t")
  expect_error(buildLiftoverMap(ga, malformed, chain, file.path(root, "bad-map")),
               "exactly one chromosome column")
  data.table::fwrite(data.frame(chrom = 1L, pos = 1L, ref = "A"),
                     malformed, sep = "\t")
  expect_error(buildLiftoverMap(ga, malformed, chain, file.path(root, "bad-map")),
               "plus pos, ref and alt")

  bundle <- packageLiftoverMap(ga, file.path(root, "bundle"))
  other <- createDB(ga$version, tempfile(tmpdir = root))
  expect_error(downloadDB(other, what = "liftover"), "No liftover bundle is configured")
  expect_true(all(file.copy(c(catalog, file.path(dirname(catalog), "marker_manifest.json")),
                           other$dirs[["marker"]])))
  expect_error(downloadDB(other, what = "liftover", archive = bundle$archive,
                          sha256 = strrep("0", 64)), "archive SHA-256 mismatch")
  other <- downloadDB(other, what = "liftover", archive = bundle$archive)
  installed <- gact:::.gact_liftover_dir(other)
  expect_identical(liftGWASStats(stat, installed), lifted)
  expect_identical(tail(downloadDB(other, "liftover")$downloads$dest_file, 1L), installed)

  remote_db <- createDB(ga$version, tempfile(tmpdir = root))
  remote_catalog <- file.path(remote_db$dirs[["marker"]], basename(catalog))
  writeLines("mismatched catalog", remote_catalog)
  expect_error(downloadDB(remote_db, "liftover", archive = bundle$archive),
               "different GRCh38 marker catalog")
  expect_false(dir.exists(gact:::.gact_liftover_dir(remote_db)))
  expect_true(file.copy(catalog, remote_catalog, overwrite = TRUE))
  local_mocked_bindings(curl_download = function(url, destfile, ...) {
    expect_identical(url, "https://example.org/liftover.tar.gz")
    stopifnot(file.copy(bundle$archive, destfile))
  }, .package = "curl")
  versions <- gact:::gact_versions
  versions[[ga$version]]$sources$liftover <- "https://example.org/liftover.tar.gz"
  versions[[ga$version]]$liftover_sha256 <- bundle$sha256
  local_mocked_bindings(.gact_version_metadata = function(version) versions[[version]],
                        .package = "gact")
  remote_db <- downloadDB(remote_db, "liftover")
  expect_true(file.exists(file.path(gact:::.gact_liftover_dir(remote_db), "index.bin")))
  writeLines("altered", file.path(map, "chr1.tsv.gz"))
  expect_error(liftGWASStats(stat, map), "checksum mismatch")
  expect_error(buildLiftoverMap(ga), "checksum mismatch")
})

test_that("native manifest aligns verified GRCh37 and GRCh38 in one batch", {
  root <- tempfile("panukb-native-")
  ga <- createDB("hsa.GRCh38.0.0.1", root)
  on.exit(unlink(root, recursive = TRUE))
  descriptor <- file.path(root, "GAlist.rds")
  catalog <- file.path(ga$dirs[["marker"]], "markers.txt.gz")
  sites <- data.frame(rsids = paste0("rs", 1:13), chr = 1L,
    pos = c(101:111, 180L, 15L), ea = c(rep("C", 11L), "C", "C"),
    nea = c(rep("A", 11L), "T", "A"), eaf = .2)
  ga$rsids <- sites$rsids
  saveRDS(ga, descriptor)
  data.table::fwrite(sites, catalog)
  data.table::fwrite(data.frame(cpra = with(sites, paste(chr, pos, ea, nea, sep = "_"))),
    file.path(dirname(catalog), "cpra.txt.gz"), col.names = FALSE)
  jsonlite::write_json(list(build = "GRCh38", effect_allele = "ALT",
    checksums = list("markers.txt.gz" = unname(tools::md5sum(catalog)))),
    file.path(dirname(catalog), "marker_manifest.json"), auto_unbox = TRUE)
  raw <- data.frame(chr = 1L, pos = c(1:11, 20L, 15L), ref = "A",
    alt = c(rep("C", 11L), "G", "C"), beta_EUR = .25, se_EUR = .1,
    neglog10_pval_EUR = 3, af_EUR = .2)
  chain <- file.path(root, "example.chain")
  writeLines(c("chain 100 chr1 200 + 0 11 chr1 200 + 100 111 1", "11", "",
    "chain 90 chr1 200 + 19 20 chr1 200 - 20 21 2", "1"), chain)
  variants <- file.path(root, "variants.tsv")
  data.table::fwrite(raw, variants, sep = "\t")
  map <- buildLiftoverMap(ga, variants, chain)
  index <- file.path(map, "index.bin")
  original_index <- file.info(index)[c("size", "mtime")]
  target <- transform(raw, pos = sites$pos, ref = sites$nea, alt = sites$ea)
  filenames <- c("original.tsv.gz", "target.tsv.gz", "mismatch.tsv.gz")
  files <- file.path(ga$dirs[["download"]], filenames)
  source <- raw
  source[2L, c("ref", "alt", "beta_EUR", "af_EUR")] <- list("C", "A", -.25, .8)
  data.table::fwrite(source, files[1L], sep = "\t")
  data.table::fwrite(target, files[2L], sep = "\t")
  data.table::fwrite(transform(raw, beta_EUR = .5), files[3L], sep = "\t")
  manifest <- data.frame(accession = c("PAN37", "PAN38", "MISMATCH"),
    filename = filenames, url = "", md5 = unname(tools::md5sum(files)),
    source = "test Pan-UKB", trait = "Test", type = "quantitative", ancestry = "EUR",
    strand = "forward", build = c("GRCh37", "GRCh38", "GRCh38"))
  mapping <- setNames(as.list(names(raw)), c("chr", "pos", "nea", "ea", "b", "seb", "p", "eaf"))
  profile <- list(source = "test Pan-UKB", headers = names(raw), mapping = mapping,
    effect_scale = "beta", p_scale = "neglog10p", verified = TRUE, evidence = "synthetic")
  knowledge <- file.path(root, "formats.json")
  jsonlite::write_json(list(formats = list(profile)), knowledge, auto_unbox = TRUE)
  run <- function(rows, ...) runGWASManifest(manifest[rows, , drop = FALSE], descriptor,
    knowledge = knowledge, engine = "cpp", ...)
  expect_identical(run(1L, import = FALSE)$status, "ready")
  expect_match(run(1L, import = FALSE, liftover_map = map, preview_rows = 3L)$issues,
    "Preview does not verify")
  result <- run(1:3, liftover_map = map)
  expect_identical(result$status, c("imported", "imported", "needs_review"))
  expect_match(result$issues[3L], "12 GRCh37-only")
  expect_equal(result$n_final[1:2], c(12, 13))
  expect_equal(result$n_flipped[1:2], c(1, 0))
  expect_identical(file.info(index)[c("size", "mtime")], original_index)
  saved <- readRDS(descriptor)
  expect_identical(saved$study$build, rep("GRCh38", 2L))
  proof <- jsonlite::read_json(result$recipe[1L])
  expect_identical(proof$build_check$build, "GRCh37")
  qc <- jsonlite::read_json(file.path(ga$dirs[["gstat"]], saved$study$qc_file[1L]))
  expect_equal(qc$n_lifted, 12)
  expect_equal(qc$dropped_by_reason$liftover_unmapped, 1)
  expect_equal(qc$n_input, qc$n_final + sum(unlist(qc$dropped_by_reason)))
  out <- getGWAS(saved, studyID = "GWAS1", rsids = sites$rsids, format = "data.frame")
  expect_equal(out$pos, sites$pos[1:12])
  expect_equal(out$b, rep(.25, 12))
  expect_equal(out$p, rep(.001, 12))
  same_build_qc <- jsonlite::read_json(file.path(ga$dirs[["gstat"]], saved$study$qc_file[2L]))
  expect_null(same_build_qc$n_lifted)
  writeLines("altered", index)
  expect_match(run(3L, import = FALSE, liftover_map = map)$issues, "index checksum mismatch")
})
