test_that("unused catalog columns and row order do not affect either LD route", {
  x <- ldscore_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  ga <- prepareLDscoresDB(x$ga, progress = FALSE)
  expected <- getLDscoresDB(ga, ancestry = "EUR", version = "1000G")
  file <- file.path(x$root, "markers.txt.gz")
  markers <- data.table::fread(file)
  markers$unused_text <- paste0("ignored-", markers$rsids)
  markers$unused_missing <- NA_real_
  data.table::fwrite(markers[rev(seq_len(nrow(markers)))], file)
  ga <- prepareLDscoresDB(x$ga, overwrite = TRUE, progress = FALSE)
  expect_identical(getLDscoresDB(ga, ancestry = "EUR", version = "1000G"), expected)
  ga <- prepareLDDB(ga, progress = FALSE)
  expect_equal(getLDscoresDB(ga, ancestry = "EUR", version = "1000G"), expected, tolerance = 1e-6)
  gl <- readRDS(ga$ldPanels$EUR$descriptor)
  expect_silent(gbase::verify_ld(gl))
})

test_that("chunked alignment checks cover the final partial chunk", {
  root <- tempfile("ld-alignment-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  count <- 100003L
  markers <- data.table::data.table(rsids = paste0("rs", seq_len(count)),
    chr = "1", pos = seq_len(count), ea = "A", nea = "G")
  catalog <- file.path(root, "markers.tsv")
  data.table::fwrite(markers[rev(seq_len(count))], catalog)
  panel <- data.frame(chr = markers$chr, id = markers$rsids, map = 0,
    pos = markers$pos, a1 = "A", a2 = "G")
  reference <- function() list(rsids = list(panel$id), chr = list(panel$chr),
    pos = list(panel$pos), a1 = list(panel$a1), a2 = list(panel$a2))
  select <- gact:::.gact_ld_select_markers
  expect_identical(select(catalog, reference(), rev(markers$rsids), TRUE), markers$rsids)
  panel$pos[count] <- panel$pos[count] + 1L
  expect_error(select(catalog, reference(), markers$rsids, TRUE), "coordinates or alleles")
  panel$pos[count] <- markers$pos[count]
  panel$a1[count] <- "G"; panel$a2[count] <- "A"
  expect_identical(select(catalog, reference(), markers$rsids, FALSE), markers$rsids)
  expect_error(select(catalog, reference(), markers$rsids, TRUE), "BIM A1/A2")
})
