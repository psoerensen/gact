#' Align a summary-statistics file to the gact marker catalog
#'
#' Match by rsid, chromosome, position and allele pair. Swap the allele labels,
#' reverse beta and Z signs, and complement effect-allele frequency when needed.
#' Repair older Z signs that equal the negative of beta divided by SE.
#' @param file An existing gact-format CSV or TSV file, plain or gzipped.
#' @param GAlist A gact descriptor pointing to the marker catalog.
#' @param out_file Output path. Defaults to replacing the source after validation.
#' @return Invisibly, the output path, recode counts and the file's MD5.
#' @export
recodeSumstat <- function(file, GAlist, out_file = file) {
  one_path <- function(x) is.character(x) && length(x) == 1L &&
    !is.na(x) && nzchar(x)
  if (!one_path(file) || !one_path(out_file))
    stop("file and out_file must each name one file.", call. = FALSE)
  file <- normalizePath(file, winslash = "/", mustWork = TRUE)
  out_file <- path.expand(out_file)
  parent <- normalizePath(dirname(out_file), winslash = "/", mustWork = TRUE)
  out_file <- file.path(parent, basename(out_file))
  if (file.exists(out_file))
    out_file <- normalizePath(out_file, winslash = "/", mustWork = TRUE)
  if (dir.exists(file) || dir.exists(out_file) ||
      (file != out_file && file.exists(out_file)))
    stop("Provide a study file and either the same file or a new output path.", call. = FALSE)
  marker_file <- file.path(GAlist$dirs[["marker"]], "markers.txt.gz")
  if (length(marker_file) != 1L || !file.exists(marker_file))
    stop("GAlist must point to an installed markers.txt.gz.", call. = FALSE)
  marker_file <- normalizePath(marker_file, winslash = "/", mustWork = TRUE)
  if (identical(out_file, marker_file))
    stop("Cannot overwrite the marker catalog.", call. = FALSE)
  stage <- tempfile(".gact-recode-", tmpdir = parent)
  if (!dir.create(stage)) stop("Cannot create a staging directory.", call. = FALSE)
  on.exit(unlink(stage, recursive = TRUE))
  counts <- .Call(C_gact_recode_sumstat, enc2utf8(file), enc2utf8(marker_file),
    enc2utf8(file.path(stage, "output")), enc2utf8(out_file),
    grepl("\\.gz$", out_file, ignore.case = TRUE))
  invisible(list(file = out_file, reversed = counts[1L], repaired_z = counts[2L],
    checksum = unname(tools::md5sum(out_file))))
}
