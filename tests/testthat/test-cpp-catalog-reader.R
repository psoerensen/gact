catalog_reader_fixture <- function(writer = "fwrite") {
  root <- tempfile("catalog-reader-")
  dir.create(root)
  markers <- data.frame(rsids = c("rs1", "rs2", "rs3"), chr = 1L,
    pos = c(1000000, 12300000, 20000), ea = c("A", "G", "A"),
    nea = c("C", "T", "T"), eaf = c(.2, .3, .4))
  catalog <- file.path(root, "markers.txt.gz")
  if (writer == "write.csv") {
    con <- gzfile(catalog, "wt")
    write.csv(markers, con, row.names = FALSE)
    close(con)
  } else data.table::fwrite(markers, catalog, sep = if (writer == "tsv") "\t" else ",")
  # Numeric spellings need not match byte for byte between CSV and paste().
  cpra <- data.frame(cpra = with(markers, paste(chr,
    format(pos, scientific = FALSE, trim = TRUE), ea, nea, sep = "_")))
  keys <- file.path(root, "cpra.txt.gz")
  if (writer == "write.csv") {
    con <- gzfile(keys, "wt")
    write.csv(cpra, con, row.names = FALSE) # Quoted header and quoted keys.
    close(con)
  } else data.table::fwrite(cpra, keys, col.names = FALSE)
  source <- data.frame(chromosome = 1L, base_pair_location = markers$pos,
    effect_allele = c("A", "T", "A"), other_allele = c("C", "G", "T"),
    beta = c(.5, .2, .3), standard_error = .1, p_value = .05,
    effect_allele_frequency = c(.2, .7, .4))
  input <- file.path(root, "study.tsv.gz")
  data.table::fwrite(source, input, sep = "\t")
  list(root = root, catalog = catalog, keys = keys, markers = markers,
       source = source, input = input, index = file.path(root, "index.bin"))
}

test_that("native indexing accepts R writer formats and equivalent CPRA coordinates", {
  hashes <- character()
  for (writer in c("fwrite", "write.csv", "tsv")) {
    x <- catalog_reader_fixture(writer)
    if (writer == "fwrite") {
      con <- gzfile(x$catalog)
      text <- readLines(con)
      close(con)
      expect_match(text[2L], "1e\\+06") # Reproduces the old integer-parser failure.
    }
    index <- get("C_gact_gwas_index", asNamespace("gact"))
    expect_null(.Call(index, x$catalog, x$keys, x$index))
    hashes <- c(hashes, unname(tools::md5sum(x$index)))
    profile <- gact:::.gact_catalog_preview_profile(x$source, "Synthetic Catalog")
    recipe <- list(mapping = profile$mapping, inspection = list(headers = names(x$source)),
                   effect_scale = "beta", p_scale = "p")
    output <- file.path(x$root, "aligned.txt.gz")
    qc_path <- file.path(x$root, "qc.json")
    expect_null(.Call(get("C_gact_gwas_stage", asNamespace("gact")), x$index, x$input,
      output, qc_path, names(x$source), gact:::.gact_cpp_gwas_mapping(recipe), TRUE, NULL, FALSE))
    actual <- data.table::fread(output)
    expect_identical(actual$rsids, x$markers$rsids)
    expect_equal(actual$pos, x$markers$pos)
    expect_equal(actual$b, c(.5, -.2, .3))
    expect_equal(actual$eaf, x$markers$eaf)
    expect_true(all(is.na(actual$n)))
    qc <- jsonlite::read_json(qc_path)
    expect_equal(c(qc$n_input, qc$n_final, qc$n_dropped, qc$n_flipped), c(3, 3, 0, 1))
    unlink(x$root, recursive = TRUE)
  }
  expect_length(unique(hashes), 1L)
})

test_that("invalid catalog values and CPRA mismatches still fail with row diagnostics", {
  x <- catalog_reader_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  index <- get("C_gact_gwas_index", asNamespace("gact"))
  bad_values <- list(pos = c("1000000.5", "0", "-1", "2147483648", "Inf", "NaN", "1e999", "12junk"),
                     chr = c("26", "1.5"), ea = c("AC", "N", "T"))
  for (field in names(bad_values)) for (value in bad_values[[field]]) {
    markers <- x$markers
    markers[[field]] <- as.character(markers[[field]])
    markers[[field]][2L] <- value
    data.table::fwrite(markers, x$catalog)
    expect_error(.Call(index, x$catalog, x$keys, x$index), "Marker catalog row 2: Invalid catalog")
    expect_false(file.exists(x$index))
  }
  data.table::fwrite(x$markers, x$catalog)
  keys <- with(x$markers, paste(chr, format(pos, scientific = FALSE, trim = TRUE), ea, nea, sep = "_"))
  keys[2L] <- "1_12300001_G_T"
  data.table::fwrite(data.frame(cpra = keys), x$keys, col.names = FALSE)
  expect_error(.Call(index, x$catalog, x$keys, x$index), "row 2: CPRA and marker catalog disagree")
  expect_false(file.exists(x$index))
})
