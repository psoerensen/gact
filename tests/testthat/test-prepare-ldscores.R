test_that("prepared scores match dense LD and the existing sparse workflow", {
  for (build in c("GRCh37", "GRCh38")) {
    x <- ldscore_fixture(build)
    tryCatch({
      inputs <- list.files(x$root, full.names = TRUE)
      before <- tools::md5sum(inputs)
      result <- prepareLDscoresDB(x$ga, progress = FALSE)
      expect_identical(tools::md5sum(inputs), before)
      expect_false(any(grepl("\\.lock$", list.files(x$root))))
      scores <- getLDscoresDB(result, ancestry = "EUR", version = "1000G")
      gl <- gbase::gprep(bedfiles = paste0(x$prefix, ".bed"), rsids = x$ga$rsids)
      sparse <- gbase::gprep(gl, task = "sparseld", msize = 1000L,
        out_prefix = file.path(x$root, "oracle"), progress = FALSE)
      expect_equal(scores, gbase::getLDscores(sparse), tolerance = 1e-6)
      # rs1 and rs2 are identical despite being >1 Mb apart. rs3 is
      # monomorphic, rs4 is rs1 with a missing first genotype.
      genotype <- cbind(c(2, 1, 0, 2), c(2, 1, 0, 2), rep(0, 4), c(NA, 1, 0, 2))
      centered <- sweep(genotype, 2, colMeans(genotype, na.rm = TRUE))
      centered[is.na(centered)] <- 0
      scale <- sqrt(colSums(centered^2))
      centered[, scale > 0] <- sweep(centered[, scale > 0], 2, scale[scale > 0], "/")
      expected <- rowSums(crossprod(centered)^2)
      expected[3] <- 1
      expect_equal(unname(scores), c(expected, 1), tolerance = 1e-6)
      expect_identical(names(scores), paste0("rs", 1:5))
      expect_equal(getLDscoresDB(result, chr = 2, ancestry = "EUR", version = "1000G"), c(rs5 = 1))
      expect_equal(getLDscoresDB(result, ancestry = "EUR", version = "1000G", rsids = "rs2"), scores["rs2"])
      receipt <- result$ldscorePanels$EUR$manifest
      meta <- jsonlite::read_json(receipt, simplifyVector = TRUE)
      expect_identical(meta$build, build)
      expect_equal(meta$n_samples, 4)
      expect_equal(meta$n_markers, 5)
      expect_equal(meta$settings$max_distance_bp, 0)
      expect_false(meta$settings$corrected)
      expect_false(any(grepl("g1000_eur.*row_ptr", list.files(x$root))))
      hashes <- tools::md5sum(list.files(x$root, full.names = TRUE))
      expect_message(again <- prepareLDscoresDB(result), "reusing verified")
      expect_identical(again, result)
      expect_identical(tools::md5sum(names(hashes)), hashes)
      subset <- result
      subset$rsids <- subset$rsids[-1]
      expect_error(prepareLDscoresDB(subset, progress = FALSE), "changed")
      wrong <- result
      wrong$metadata$build <- if (build == "GRCh38") "GRCh37" else "GRCh38"
      expect_error(getLDscoresDB(wrong, ancestry = "EUR", version = "1000G"), "build or ancestry")
      cat("corrupt", file = result$ldscorePanels$EUR$file, append = TRUE)
      expect_error(getLDscoresDB(result, ancestry = "EUR", version = "1000G"), "checksum")
      expect_error(prepareLDscoresDB(result, progress = FALSE), "checksum")
      expect_silent(prepareLDscoresDB(result, overwrite = TRUE, progress = FALSE))
    }, finally = unlink(x$root, recursive = TRUE))
  }
})

test_that("the legacy window includes distance 1000 but excludes distance 1001", {
  x <- ldscore_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  count <- 1002L
  markers <- data.frame(rsids = paste0("rs", seq_len(count)), chr = 1L,
    pos = seq_len(count) * 100L, ea = "A", nea = "G")
  writeBin(as.raw(c(0x6c, 0x1b, 0x01, rep(0x38, count))), paste0(x$prefix, ".bed"))
  writeLines(with(markers, paste(chr, rsids, 0, pos, ea, nea)), paste0(x$prefix, ".bim"))
  data.table::fwrite(markers, file.path(x$root, "markers.txt.gz"))
  x$ga$rsids <- markers$rsids
  ga <- prepareLDscoresDB(x$ga, progress = FALSE)
  scores <- getLDscoresDB(ga, ancestry = "EUR", version = "1000G")
  expect_equal(unname(scores), c(1001, rep(1002, 1000), 1001), tolerance = 1e-6)
})

test_that("legacy GRCh37 tables remain readable and GRCh38 requires provenance", {
  x <- ldscore_fixture("GRCh37")
  on.exit(unlink(x$root, recursive = TRUE))
  path <- file.path(x$root, "markers_1000G_eur_filtered.txt.gz")
  data.table::fwrite(data.frame(rsids = "rs1", chr = 1, ldscores = 1), path)
  expect_equal(getLDscoresDB(x$ga, ancestry = "EUR", version = "1000G"), c(rs1 = 1))
  x$ga$metadata$build <- "GRCh38"
  expect_error(getLDscoresDB(x$ga, ancestry = "EUR", version = "1000G"), "require provenance")
  expect_error(prepareLDscoresDB(x$ga, progress = FALSE), "no provenance")
  unlink(path)
  expect_error(getLDscoresDB(x$ga, ancestry = "EUR", version = "1000G"), "prepareLDscoresDB")
})

test_that("preparation validates panels and supports all five reference ancestries", {
  x <- ldscore_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  expect_error(prepareLDscoresDB(x$ga, ancestry = "MID"), "ancestries")
  expect_error(prepareLDscoresDB(x$ga, nthreads = 0), "positive integer")
  expect_error(prepareLDscoresDB(x$ga, ancestry = "AFR"), "Missing reference")
  for (pop in c("afr", "eas", "sas", "amr"))
    for (suffix in c(".bed", ".bim", ".fam"))
      file.copy(paste0(x$prefix, suffix), file.path(x$root, paste0("g1000_", pop, suffix)))
  ga <- prepareLDscoresDB(x$ga, ancestry = c("EUR", "AFR", "EAS", "SAS", "AMR"), progress = FALSE)
  for (pop in names(ga$ldscorePanels))
    expect_equal(getLDscoresDB(ga, ancestry = pop, version = "1000G"),
      getLDscoresDB(ga, ancestry = "EUR", version = "1000G"))
  bim <- paste0(x$prefix, ".bim")
  original <- readLines(bim)
  writeLines(sub(" 10 ", " 11 ", original), bim)
  expect_error(prepareLDscoresDB(x$ga, overwrite = TRUE, progress = FALSE), "coordinates or alleles")
  # Restore the original positions, then verify that allele swaps are allowed.
  writeLines(sub(" A G$", " G A", original), bim)
  expect_silent(prepareLDscoresDB(x$ga, overwrite = TRUE, progress = FALSE))
})
