#' Register Berisa–Pickrell LD regions
#'
#' Installs the published EUR and EAS GRCh37 interval files in the marker
#' directory. Prefer `downloadDB(GAlist, what = "ldetect")` for database setup.
#' EAS is labelled ASN upstream. Starts are inclusive and stops are exclusive.
#' @param GAlist A gact database descriptor.
#' @param assembly Genome assembly of the marker catalog.
#' @return The updated descriptor with file paths and provenance in
#'   GAlist$ldRegions.
#' @export
addLDRegions <- function(GAlist, assembly = "GRCh37") {
  if (!identical(assembly, "GRCh37"))
    stop("Bundled Berisa–Pickrell regions are available for GRCh37.", call. = FALSE)
  if (!is.null(GAlist$assembly) && !identical(GAlist$assembly, assembly))
    stop("LD region assembly differs from GAlist$assembly.", call. = FALSE)
  marker_dir <- GAlist$dirs[["marker"]]
  if (length(marker_dir) != 1L || is.na(marker_dir) || !dir.exists(marker_dir))
    stop("GAlist must point to an existing marker directory.", call. = FALSE)
  marker_dir <- normalizePath(marker_dir, winslash = "/", mustWork = TRUE)
  revision <- "ac125e47bf7ff3e90be31f278a7b6a61daaba0dc"
  populations <- c(EUR = "EUR", EAS = "ASN")
  entries <- lapply(names(populations), function(ancestry) {
    name <- paste0(assembly, "-", ancestry, ".bed")
    bundled <- system.file("extdata", "ld-regions", name, package = "gact", mustWork = TRUE)
    path <- file.path(marker_dir, paste0("ldetect-", name))
    md5 <- unname(tools::md5sum(bundled))
    if (!file.exists(path)) {
      staged <- tempfile(".ldetect-", tmpdir = marker_dir)
      on.exit(unlink(staged))
      if (!file.copy(bundled, staged) || !file.rename(staged, path))
        stop("Cannot install LD region file: ", path, call. = FALSE)
    }
    if (!identical(unname(tools::md5sum(path)), md5))
      stop("Installed LD region file differs from the bundled file: ", path, call. = FALSE)
    list(file = path,
         source = sprintf("https://bitbucket.org/nygcresearch/ldetect-data/raw/%s/%s/fourier_ls-all.bed",
                          revision, populations[[ancestry]]),
         revision = revision, md5 = md5)
  })
  GAlist$assembly <- assembly
  GAlist$ldRegions[[assembly]] <- setNames(entries, names(populations))
  GAlist
}

.stage_ld_regions <- function(GAlist, stage) {
  panels <- names(GAlist$ldPanels)
  files <- rep("", length(panels))
  assembly <- GAlist$assembly
  if (!length(GAlist$ldRegions) || !length(panels)) return(files)
  if (!is.character(assembly) || length(assembly) != 1L || is.na(assembly) ||
      !grepl("^[A-Za-z0-9]+$", assembly))
    stop("Set GAlist$assembly to select the LD regions.", call. = FALSE)
  regions <- GAlist$ldRegions[[assembly]]
  if (is.null(regions)) stop("No LD regions for GAlist$assembly.", call. = FALSE)
  for (i in seq_along(panels)) {
    entry <- regions[[panels[i]]]
    if (is.null(entry)) next
    panel_assembly <- GAlist$ldPanels[[i]]$assembly
    if (!is.null(panel_assembly) && !identical(panel_assembly, assembly))
      stop("LD assembly differs from GAlist$assembly for ", panels[i], ".", call. = FALSE)
    files[i] <- paste0("ld-regions-", assembly, "-", panels[i], ".bed")
    intervals <- if (is.null(entry$file)) entry$intervals else
      utils::read.table(entry$file, header = TRUE,
                        colClasses = c("character", "integer", "integer"))
    if (!is.data.frame(intervals) || !identical(names(intervals), c("chr", "start", "stop")))
      stop("LD intervals must have chr, start and stop columns.", call. = FALSE)
    utils::write.table(intervals, file.path(stage, files[i]), sep = "\t",
                       quote = FALSE, row.names = FALSE)
  }
  files
}
