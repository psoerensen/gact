#' Register custom sparse LD panels for marker statistics
#'
#' Records the path to each panel's gmat sparse LD marker file and the panel's
#' assembly when supplied. During [writeMarkerStat()], each registered panel
#' takes precedence over the conventional filename for that ancestry.
#' @param GAlist A gact database descriptor.
#' @param Glist An ancestry-named list of aligned Glist objects with sparse LD,
#'   for example \code{list(EUR = eur, EAS = eas, SAS = sas)}.
#' @return The updated descriptor with reference paths in GAlist$ldPanels.
#' @export
addLDPanels <- function(GAlist, Glist) {
  if (!length(Glist)) {
    GAlist$ldPanels <- list()
    return(GAlist)
  }
  panels <- .gact_ancestry_codes(names(Glist))
  if (!is.list(Glist) || length(panels) != length(Glist) || anyNA(panels) ||
      anyDuplicated(panels) || any(!grepl("^[A-Z][A-Z0-9_-]*$", panels)))
    stop("Name each Glist by its ancestry, for example EUR.", call. = FALSE)
  references <- lapply(seq_along(Glist), function(i) {
    panel <- Glist[[i]]
    ld_files <- if (is.list(panel)) panel$ldfiles
    path <- if ("markers" %in% names(ld_files)) ld_files[["markers"]]
    if (!is.character(path) || length(path) != 1L || is.na(path) ||
        !file.exists(path) || dir.exists(path))
      stop("Missing gmat sparse LD marker file for ", panels[i], ".", call. = FALSE)
    list(markers = normalizePath(path, winslash = "/", mustWork = TRUE),
         assembly = panel$assembly)
  })
  GAlist$ldPanels <- setNames(references, panels)
  GAlist
}

.find_aligned_ld <- function(marker_dir, study_ancestry, panels) {
  ancestries <- setdiff(intersect(c("EUR", "EAS", "SAS"), study_ancestry), names(panels))
  suffixes <- c(".markers.tsv", ".meta.txt", ".row_ptr.u64.bin",
                ".col_idx.u32.0based.bin", ".values.f32.bin")
  for (ancestry in ancestries) {
    stem <- paste0("g1000_", tolower(ancestry), ".subset.aligned")
    paths <- file.path(marker_dir, paste0(stem, suffixes))
    present <- file.exists(paths)
    if (!any(present)) next
    missing <- !present | dir.exists(paths)
    if (any(missing))
      stop("Incomplete sparse LD for ", ancestry, ". Missing files: ",
           paste(basename(paths[missing]), collapse = ", "), call. = FALSE)
    panels[[ancestry]] <- list(markers = normalizePath(paths[1L], winslash = "/",
                                                        mustWork = TRUE))
  }
  panels
}
