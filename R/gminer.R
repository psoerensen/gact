#' Run gminer on a prepared association matrix
#' @param GAlist A gact descriptor.
#' @param task Native task. Currently only \code{"impute"} is supported.
#' @param matrix_dir Directory written by \code{\link{assocMatrix}} with LD
#'   maps. Defaults to \code{gwas_matrix} next to the gact marker directory.
#' @param out_dir New directory for output. Defaults to \code{imputed}
#'   within \code{matrix_dir}.
#' @return Invisibly, paths to the filled Z matrix, mask, predicted r-squared
#'   matrix, per-study counts and task metadata.
#' @export
gminer <- function(GAlist, task = "impute", matrix_dir = NULL, out_dir = NULL) {
  if (!identical(task, "impute"))
    stop("The supported gminer task is \"impute\".", call. = FALSE)
  if (is.null(matrix_dir)) {
    marker_dir <- GAlist$dirs[["marker"]]
    if (length(marker_dir) != 1L || is.na(marker_dir) || !nzchar(marker_dir))
      stop("Supply a GAlist marker directory or matrix_dir.", call. = FALSE)
    matrix_dir <- file.path(dirname(marker_dir), "gwas_matrix")
  }
  matrix_dir <- normalizePath(matrix_dir, winslash = "/", mustWork = TRUE)
  if (is.null(out_dir)) out_dir <- file.path(matrix_dir, "imputed")
  out_dir <- path.expand(out_dir)
  out_dir <- file.path(normalizePath(dirname(out_dir), winslash = "/", mustWork = TRUE),
                        basename(out_dir))
  exe <- Sys.which("gminer_impute")
  if (!nzchar(exe))
    stop("Put the gminer_impute executable on PATH.", call. = FALSE)
  status <- system2(exe, shQuote(c("--assoc", matrix_dir, out_dir)))
  if (!identical(status, 0L))
    stop("gminer imputation failed (exit status ", status, ").", call. = FALSE)
  files <- c(z = "z.bin", imputed = "imputed.bin", r2 = "r2.bin",
             counts = "counts.tsv", task = "task.tsv")
  invisible(as.list(setNames(file.path(out_dir, files), names(files))))
}
