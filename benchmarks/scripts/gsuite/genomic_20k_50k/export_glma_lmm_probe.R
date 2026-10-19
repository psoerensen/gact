# Export the frozen input order for bounded native diagnostics; no fitting.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
x <- readRDS(args[[1L]])
g <- x$Glist
stopifnot(length(x$y) == 5000L, length(g$bedfiles) == 1L,
  identical(as.integer(g$sample_source_indices_0based), 0:4999))
dir.create(args[[2L]], recursive = TRUE, showWarnings = FALSE)
writeLines(format(x$y, digits = 17, scientific = TRUE),
  file.path(args[[2L]], "frozen.phen"))
writeLines(unlist(g$rsids, use.names = FALSE),
  file.path(args[[2L]], "frozen.markers"))
writeLines(c(g$bedfiles, g$bimfiles, g$famfiles),
  file.path(args[[2L]], "frozen.paths"))
cat("Frozen phenotype and 37991-marker order exported\n")
