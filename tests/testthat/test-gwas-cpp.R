cpp_gwas_fixture <- function() {
  root <- tempfile("gact-cpp-")
  dir.create(root)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  dirs <- setNames(file.path(root, c("marker", "gstat", "download")),
                   c("marker", "gstat", "download"))
  for (dir in dirs) dir.create(dir)
  marker <- data.frame(rsids = c("rs1", "rs2", "rs3"), chr = 1L, pos = 101:103,
                       ea = c("A", "C", "A"), nea = c("C", "T", "T"), eaf = c(.2, .3, .4))
  data.table::fwrite(marker, file.path(dirs[["marker"]], "markers.txt.gz"))
  data.table::fwrite(data.frame(cpra = with(marker, paste(chr, pos, ea, nea, sep = "_"))),
    file.path(dirs[["marker"]], "cpra.txt.gz"), col.names = FALSE)
  ga <- list(dirs = dirs, rsids = marker$rsids, study = list(id = character()),
             metadata = list(build = "GRCh38"))
  descriptor <- file.path(root, "GAlist.rds")
  saveRDS(ga, descriptor)
  source <- data.frame(chromosome = 1L, base_pair_location = 101:103, effect_allele = c("A", "T", "A"),
    other_allele = c("C", "C", "T"), beta = c(.5, .2, .3), standard_error = .1,
    effect_allele_frequency = c(.2, .7, .4), p_value = .05, rsid = marker$rsids)
  filenames <- c("a.h.tsv.gz", "b.h.tsv.gz")
  files <- file.path(dirs[["download"]], filenames)
  data.table::fwrite(source, files[1L], sep = "\t")
  data.table::fwrite(transform(source, beta = beta * 2), files[2L], sep = "\t")
  mapping <- setNames(as.list(names(source)), c("chr", "pos", "ea", "nea", "b", "seb", "eaf", "p", "marker"))
  knowledge <- file.path(root, "formats.json")
  formats <- list(list(source = "test SSF", headers = names(source), mapping = mapping,
    effect_scale = "beta", p_scale = "p", verified = TRUE, evidence = "synthetic"))
  jsonlite::write_json(list(formats = formats), knowledge, auto_unbox = TRUE)
  manifest <- data.frame(accession = c("GCST123", "GCST456"), filename = filenames, url = "",
    md5 = unname(tools::md5sum(files)), source = "test SSF", trait = c("A", "B"),
    type = "quantitative", build = "GRCh38", ancestry = "Mixed", strand = "forward",
    study_n_total = c(1200L, 1000L))
  list(root = root, ga = ga, descriptor = descriptor, manifest = manifest, knowledge = knowledge,
       files = files, mapping = mapping)
}

test_that("native staging publishes verified recipes and reads both studies", {
  x <- cpp_gwas_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  scratch <- file.path(x$root, "scratch")
  dir.create(scratch)
  writeLines("owned by the caller", file.path(scratch, "keep.me"))
  inputs <- x$files
  inputs[1L] <- file.path(x$root, "upstream.tsv.gz")
  expect_true(file.rename(x$files[1L], inputs[1L]))
  x$manifest$url[1L] <- "https://example.org/a.h.tsv.gz"
  testthat::local_mocked_bindings(curl_download = function(url, destfile, mode) {
    expect_identical(url, x$manifest$url[1L])
    expect_true(file.copy(inputs[1L], destfile))
  }, .package = "curl")
  run <- function(name, import = TRUE) runGWASManifest(x$manifest, x$descriptor,
    out_dir = if (is.null(name)) NULL else file.path(x$root, name),
    knowledge = x$knowledge, engine = "cpp", workers = 2L, import = import,
    scratch_dir = scratch)
  original <- gact:::.publish_gwas_import_file
  fail <- TRUE
  writes <- 0L
  testthat::local_mocked_bindings(.publish_gwas_import_file = function(from, to, ...) {
    if (identical(to, x$descriptor)) {
      work <- setdiff(list.files(scratch, full.names = TRUE), file.path(scratch, "keep.me"))
      expect_length(work, 1L)
      expect_true(file.exists(file.path(work, "marker-index.bin")))
      expect_length(list.files(work, pattern = "\\.txt\\.gz$"), 2L)
      writes <<- writes + 1L
      if (fail) stop("Simulated batch descriptor failure")
    }
    original(from, to, ...)
  }, .package = "gact")
  expect_identical(run("preview", FALSE)$status, rep("ready", 2))
  expect_false(dir.exists(file.path(x$root, "recipes")))
  expect_identical(list.files(scratch), "keep.me")
  expect_true(all(file.exists(x$files)))
  failed <- run("failure")
  expect_identical(failed$status, rep("needs_review", 2))
  expect_match(failed$issues, "Simulated batch descriptor failure", all = TRUE)
  expect_length(list.files(x$ga$dirs[["gstat"]], all.files = TRUE, no.. = TRUE), 0L)
  expect_length(readRDS(x$descriptor)$study$id, 0L)
  expect_identical(list.files(scratch), "keep.me")
  expect_true(all(file.exists(x$files)))
  fail <- FALSE
  writes <- 0L
  result <- run(NULL)
  expect_identical(result$status, rep("imported", 2))
  expect_equal(result$n_flipped, c(1, 1))
  expect_equal(writes, 1L)
  expect_true(all(file.exists(result$recipe)))
  reviews <- normalizePath(file.path(x$root, "recipes", "reviews"), winslash = "/")
  reports <- attr(result, "out_dir")
  expect_identical(dirname(reports), reviews)
  expect_true(all(file.exists(file.path(reports, c("manifest.tsv", "manifest_review.tsv")))))
  expect_length(list.files(reports, recursive = TRUE, pattern = "\\.(gz|bin)$"), 0L)
  expect_identical(list.files(scratch), "keep.me")
  expect_identical(result$cleanup, c("removed", NA_character_))
  expect_identical(file.exists(x$files), c(FALSE, TRUE))
  ga <- readRDS(x$descriptor)
  expect_identical(ga$study$id, c("GWAS1", "GWAS2"))
  expect_equal(ga$study$study_n_total, c(1200, 1000))
  for (i in seq_along(inputs)) {
    source <- data.table::fread(inputs[i], data.table = FALSE)
    normalized <- normalizeStatSchema(source, x$mapping, effect_scale = "beta", p_scale = "p")
    expected <- gact:::qcStatDB(x$ga, normalized, excludeMAFDIFF = .5, strand = "forward")
    actual <- getGWAS(ga, ga$study$id[i], rsids = x$ga$rsids, format = "data.frame")
    expect_equal(lapply(names(expected), function(field) actual[[field]]),
                 lapply(names(expected), function(field) expected[[field]]))
    expect_true(all(is.na(actual$n)))
  }
  expect_identical(run("repeat")$status, rep("already_imported", 2))
})
