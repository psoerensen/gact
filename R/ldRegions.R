#' Register Berisa–Pickrell LD regions
#'
#' Adds the published EUR and EAS GRCh37 interval tables. The EAS source is
#' labelled ASN upstream. Intervals use inclusive starts and exclusive stops.
#' @param GAlist A gact database descriptor.
#' @param assembly Genome assembly of the marker catalog.
#' @return The updated descriptor with intervals in GAlist$ldRegions.
#' @export
addLDRegions <- function(GAlist, assembly = "GRCh37") {
  if (!identical(assembly, "GRCh37"))
    stop("Bundled Berisa–Pickrell regions are available for GRCh37.", call. = FALSE)
  if (!is.null(GAlist$assembly) && !identical(GAlist$assembly, assembly))
    stop("LD region assembly differs from GAlist$assembly.", call. = FALSE)
  revision <- "ac125e47bf7ff3e90be31f278a7b6a61daaba0dc"
  populations <- c(EUR = "EUR", EAS = "ASN")
  entries <- lapply(names(populations), function(ancestry) {
    name <- paste0(assembly, "-", ancestry, ".bed")
    path <- system.file("extdata", "ld-regions", name, package = "gact", mustWork = TRUE)
    list(intervals = utils::read.table(path, header = TRUE,
           colClasses = c("character", "integer", "integer")),
         source = sprintf("https://bitbucket.org/nygcresearch/ldetect-data/raw/%s/%s/fourier_ls-all.bed",
                          revision, populations[[ancestry]]),
         revision = revision, md5 = unname(tools::md5sum(path)))
  })
  GAlist$assembly <- assembly
  GAlist$ldRegions[[assembly]] <- setNames(entries, names(populations))
  GAlist
}

.stage_ld_regions <- function(GAlist, panels, Glists, stage) {
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
    panel_assembly <- Glists[[i]]$assembly
    if (!is.null(panel_assembly) && !identical(panel_assembly, assembly))
      stop("LD assembly differs from GAlist$assembly for ", panels[i], ".", call. = FALSE)
    files[i] <- paste0("ld-regions-", assembly, "-", panels[i], ".bed")
    intervals <- entry$intervals
    if (!is.data.frame(intervals) || !identical(names(intervals), c("chr", "start", "stop")))
      stop("LD intervals must have chr, start and stop columns.", call. = FALSE)
    utils::write.table(intervals, file.path(stage, files[i]), sep = "\t",
                       quote = FALSE, row.names = FALSE)
  }
  files
}
