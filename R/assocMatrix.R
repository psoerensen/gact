#' Export a GWAS association matrix
#'
#' Writes little-endian IEEE FP64 values in column-major order to `z.f64.bin`
#' and creates `markers.tsv` and `studies.tsv` in the same directory.
#' Rows follow `GAlist$rsids` and columns follow the study inventory.
#' The writer aligns study scores to the gact marker catalog. It computes Z
#' from finite beta and positive finite SE when available, with stored Z as a
#' fallback. Missing and unusable scores are stored as NaN.
#' All studies use the gact catalog effect allele. Reference panels must use
#' that same dosage allele. Each panel's binary map indexes Z rows in LD order.
#' @param GAlist A gact database descriptor with marker and GWAS data installed.
#' @param out_dir A new output directory. Its parent must already exist.
#' @param Glist Optional list of Glist objects named by ancestry, for example
#'   list(EUR = eur, EAS = eas, SAS = sas), with gmat sparse LD prepared.
#'   Panels must use the same genome build, strand and allele order as the
#'   gact catalog for shared markers.
#' @param stat Statistic to export. The default, `"z"`, exports Z scores and
#'   is the only supported value.
#' @return Invisibly, a list with the matrix and catalog paths, matrix dimensions,
#'   and any LD row map paths.
#' @export
assocMatrix <- function(GAlist, out_dir, Glist = NULL, stat = "z") {
  if (!identical(stat, "z"))
    stop("stat must be \"z\".", call. = FALSE)

  valid_ids <- function(x) is.character(x) && length(x) > 0L &&
    !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
  rsids <- GAlist$rsids
  if (!is.character(rsids) || !length(rsids))
    stop("GAlist$rsids must contain the marker backbone.", call. = FALSE)
  studies <- GAlist$studies
  if (is.null(studies$id)) studies <- GAlist$study
  ids <- as.character(studies$id)
  if (!valid_ids(ids))
    stop("GAlist must contain unique, nonempty study IDs.", call. = FALSE)

  files <- GAlist$studyfiles
  at <- match(ids, names(files))
  if (anyNA(at) || anyDuplicated(names(files)[names(files) %in% ids]))
    stop("GAlist$studyfiles must identify one file per study ID.", call. = FALSE)
  files <- unname(files[at])
  if (anyNA(files) || any(!file.exists(files) | dir.exists(files)))
    stop("Missing aligned study file for ", ids[which(
      is.na(files) | !file.exists(files) | dir.exists(files))[1L]], ".", call. = FALSE)
  marker_dir <- GAlist$dirs[["marker"]]
  marker_file <- file.path(marker_dir, "markers.txt.gz")
  if (length(marker_file) != 1L || !file.exists(marker_file))
    stop("GAlist must point to an installed markers.txt.gz.", call. = FALSE)
  ancestry <- as.character(studies$ancestry)
  if (!length(ancestry)) ancestry <- rep(NA_character_, length(ids))
  if (length(ancestry) != length(ids))
    stop("The study ancestry column must match the study IDs.", call. = FALSE)
  ancestry <- .gact_ancestry_codes(ancestry)
  ancestry[is.na(ancestry) | !nzchar(ancestry)] <- "unknown"
  if (any(grepl("[\t\r\n]", c(ids, ancestry))))
    stop("Study IDs and ancestries must be single TSV fields.", call. = FALSE)

  panels <- panel_files <- character()
  if (length(Glist)) {
    panels <- .gact_ancestry_codes(names(Glist))
    if (!is.list(Glist) || length(panels) != length(Glist) ||
        !valid_ids(panels) || any(!grepl("^[A-Z][A-Z0-9_-]*$", panels)))
      stop("Name each Glist by its ancestry, for example EUR.", call. = FALSE)
    panel_files <- vapply(seq_along(Glist), function(i) {
      panel <- Glist[[i]]
      ld_files <- if (is.list(panel)) panel$ldfiles
      path <- if ("markers" %in% names(ld_files)) ld_files[["markers"]]
      if (!is.character(path) || length(path) != 1L || is.na(path) ||
          !file.exists(path) || dir.exists(path))
        stop("Missing gmat sparse LD marker file for ", panels[i], ".", call. = FALSE)
      normalizePath(path, winslash = "/", mustWork = TRUE)
    }, character(1L))
    if (any(grepl("[\t\r\n]", panel_files)))
      stop("LD marker paths must be single TSV fields.", call. = FALSE)
  }

  if (!is.character(out_dir) || length(out_dir) != 1L ||
      is.na(out_dir) || !nzchar(out_dir))
    stop("out_dir must be a single directory path.", call. = FALSE)
  out_dir <- path.expand(out_dir)
  parent <- normalizePath(dirname(out_dir), winslash = "/", mustWork = TRUE)
  out_dir <- file.path(parent, basename(out_dir))
  if (file.exists(out_dir) || dir.exists(out_dir))
    stop("The Z matrix output directory already exists.", call. = FALSE)
  stage <- tempfile(".gact-z-", tmpdir = parent)
  if (!dir.create(stage)) stop("Cannot create the Z matrix staging directory.", call. = FALSE)
  on.exit(unlink(stage, recursive = TRUE))
  .Call(C_gact_write_z, enc2utf8(rsids),
        enc2utf8(normalizePath(marker_file, winslash = "/", mustWork = TRUE)),
        enc2utf8(normalizePath(files, winslash = "/", mustWork = TRUE)),
        enc2utf8(ids), enc2utf8(ancestry), enc2utf8(stage),
        enc2utf8(panels), enc2utf8(panel_files))
  if (file.exists(out_dir) || !file.rename(stage, out_dir))
    stop("Cannot publish the Z matrix directory: ", out_dir, call. = FALSE)
  result <- list(z = file.path(out_dir, "z.f64.bin"),
    markers = file.path(out_dir, "markers.tsv"), studies = file.path(out_dir, "studies.tsv"),
    dim = c(length(rsids), length(ids)))
  if (length(panels)) {
    result$ld <- file.path(out_dir, "ld.tsv")
    result$ld_maps <- setNames(file.path(out_dir, paste0("ld-", panels, ".u32le.bin")), panels)
  }
  invisible(result)
}

#' @rdname assocMatrix
#' @export
writeZMatrix <- function(GAlist, out_dir, Glist = NULL) {
  assocMatrix(GAlist, out_dir, Glist, stat = "z")
}
