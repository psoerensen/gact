mining_fixture <- function() {
  root <- tempfile("mine-")
  dir.create(root)
  data <- file.path(root, "data")
  dir.create(data)
  stat <- data.frame(chr = 1L, pos = 1001:1020, ea = "A", nea = "G",
                     b = 0.2, seb = 0.1, p = 0.05)
  con <- gzfile(file.path(data, "trait.tsv.gz"), "wt")
  write.table(stat, con, sep = "\t", row.names = FALSE, quote = FALSE)
  close(con)
  meta <- data.frame(filename = "trait.tsv.gz", source = "example-release",
    trait = "example", type = "quantitative", build = "GRCh38", ancestry = "EUR")
  list(root = root, data = data, meta = meta)
}

test_that("mining reads bounded compressed previews and prepares scripts", {
  x <- mining_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  original <- tools::md5sum(list.files(x$data, full.names = TRUE))
  out <- file.path(x$root, "output")
  result <- inspectGWAS(x$data, x$meta, out, preview_rows = 3)
  expect_equal(result$review$preview_rows, 3)
  expect_true(result$review$preview_ok)
  recipe <- jsonlite::read_json(result$recipes[[1]])
  expect_identical(recipe$mapping$b, "b")
  expect_identical(recipe$build, "GRCh38")
  expect_equal(recipe$inspection$mtime, as.numeric(file.info(recipe$file_name)$mtime))
  expect_silent(parse(result$script))
  expect_identical(tools::md5sum(list.files(x$data, full.names = TRUE)), original)
  expect_error(inspectGWAS(x$data, x$meta, out), "new or empty")
})

test_that("verified source-specific knowledge is reused", {
  x <- mining_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  first <- file.path(x$root, "first")
  inspectGWAS(x$data, x$meta, first)
  path <- file.path(first, "formats.json")
  book <- jsonlite::read_json(path)
  book$formats[[1]]$verified <- TRUE
  book$formats[[1]]$evidence <- "Confirmed in release documentation"
  jsonlite::write_json(book, path, auto_unbox = TRUE)
  second <- inspectGWAS(x$data, x$meta, file.path(x$root, "second"), knowledge = path)
  expect_true(jsonlite::read_json(second$recipes[[1]])$inspection$profile_reused)
  x$meta$source <- "different-release"
  third <- inspectGWAS(x$data, x$meta, file.path(x$root, "third"), knowledge = path)
  expect_false(jsonlite::read_json(third$recipes[[1]])$inspection$profile_reused)
})

test_that("ambiguous matches and malformed files remain visible", {
  x <- mining_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  writeLines("empty", file.path(x$data, "bad.tsv"))
  result <- inspectGWAS(x$data, rbind(x$meta, x$meta), file.path(x$root, "out"))
  expect_equal(nrow(result$review), 2L)
  expect_match(result$review$issues[1], "no data rows")
  expect_match(result$review$issues[2], "Metadata matches: 2")
  expect_equal(nrow(result$unmatched_metadata), 2L)
  expect_error(inspectGWAS(x$data, x$meta, file.path(x$data, "out")), "outside")
})

test_that("generated script passes full data and recipe to the ingestion API", {
  x <- mining_fixture()
  on.exit(unlink(x$root, recursive = TRUE))
  result <- inspectGWAS(x$data, x$meta, file.path(x$root, "out"), preview_rows = 2)
  env <- new.env(parent = globalenv())
  env$library <- function(...) invisible(NULL)
  env$commandArgs <- function(...) file.path(x$root, "database.rds")
  env$readRDS <- function(...) list()
  env$saveRDS <- function(...) invisible(NULL)
  calls <- list()
  env$ingestStatDB <- function(GAlist, stat, ...) {
    args <- list(...)
    normalized <- normalizeStatSchema(stat, args$schema,
      effect_scale = args$effect_scale, p_scale = args$p_scale)
    expect_true(validateStatSchema(normalized)$ok)
    calls[[length(calls) + 1L]] <<- list(rows = nrow(stat), args = args)
    list(GAlist = GAlist)
  }
  eval(parse(result$script), env)
  expect_length(calls, 1L)
  expect_equal(calls[[1]]$rows, 20L)
  expect_identical(calls[[1]]$args$build, "GRCh38")
  con <- gzfile(file.path(x$data, "trait.tsv.gz"), "at")
  writeLines("1\t99999\tA\tG\t0.2\t0.1\t0.05", con)
  close(con)
  expect_error(eval(parse(result$script), env), "File changed")
  expect_length(calls, 1L)
})
