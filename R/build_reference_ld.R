# LD panel registration and matching --------

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

# LD region registration and staging --------

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

# LD score retrieval --------

#' Retrieve Linkage Disequilibrium (LD) Scores for Specific Genetic Markers
#'
#' This function retrieves LD scores for specified genetic markers based on the population ancestry and version of the reference panel. The function reads in marker files and extracts LD scores, allowing for optional filtering by specific rsIDs.
#'
#' @param GAlist A list containing directory paths where marker data is stored. It is assumed that `GAlist$dirs["marker"]` points to the correct directory.
#' @param ancestry A character string specifying the ancestry population. The options are "EUR" (European), "EAS" (East Asian), or "SAS" (South Asian). Default is "EUR".
#' @param version A character string specifying the version of the reference panel. The options are "HapMap3", "1000G", or "Original". Default is "HapMap3".
#' @param rsids An optional character vector of rsIDs to filter the LD scores. If not specified, all LD scores will be returned.
#'
#' @return A named vector of LD scores, where the names are rsIDs.
#'
#' @details
#' The function reads in precomputed LD score files for the selected ancestry and version, extracts the LD scores, and returns them as a named vector. If the `rsids` argument is provided, the returned vector will only include LD scores for the specified rsIDs.
#'
#' @examples
#' # Example usage:
#' GAlist <- list(dirs = c(marker = "/path/to/marker/files"))
#' ldscores <- getLDscoresDB(GAlist, ancestry = "EUR", version = "HapMap3", rsids = c("rs123", "rs456"))
#'
#' @export
#'
getLDscoresDB <- function(GAlist=NULL, chr=NULL, ancestry="EUR", version="HapMap3", rsids=NULL) {
 # Check if GAlist is provided
 if (is.null(GAlist) || !("marker" %in% names(GAlist$dirs))) {
  stop("Error: GAlist or GAlist$dirs['marker'] is missing.")
 }

 # Check if ancestry and version are valid
 valid_ancestries <- c("EUR", "EAS", "SAS")
 valid_versions <- c("HapMap3", "1000G", "Original")

 if (!(ancestry %in% valid_ancestries)) {
  stop("Error: Invalid ancestry. Choose from 'EUR', 'EAS', or 'SAS'.")
 }

 if (!(version %in% valid_versions)) {
  stop("Error: Invalid version. Choose from 'HapMap3', '1000G', or 'Original'.")
 }

 # Define the file path based on ancestry and version
 file_name <- switch(
  paste(ancestry, version, sep = "_"),
  "EUR_HapMap3" = "markers_1000G_eur_hm3.txt.gz",
  "EAS_HapMap3" = "markers_1000G_eas_hm3.txt.gz",
  "SAS_HapMap3" = "markers_1000G_sas_hm3.txt.gz",
  "EUR_1000G" = "markers_1000G_eur_filtered.txt.gz",
  "EAS_1000G" = "markers_1000G_eas_filtered.txt.gz",
  "SAS_1000G" = "markers_1000G_sas_filtered.txt.gz",
  "EUR_Original" = "markers_1000G_eur_w_ld.txt.txt.gz",
  stop("Error: No marker file found for the selected ancestry and version.")
 )

 # Read the marker file
 marker_file_path <- file.path(GAlist$dirs["marker"], file_name)

 if (!file.exists(marker_file_path)) {
  stop("Error: The marker file does not exist at the specified location.")
 }

 marker <- fread(marker_file_path, data.table = FALSE)
 if(!is.null(chr)) marker <- marker[marker$chr==chr,]
 # Check if 'ldscores' and 'rsids' columns exist
 if (!all(c("ldscores", "rsids") %in% colnames(marker))) {
  stop("Error: The marker file does not contain the required columns 'ldscores' and 'rsids'.")
 }

 # Extract ldscores and assign rsids as names
 ldscores <- marker[,"ldscores"]
 names(ldscores) <- marker$rsids

 # Filter by rsids if provided
 if (!is.null(rsids)) {
  ldscores <- ldscores[names(ldscores) %in% rsids]
 }

 return(ldscores)
}

# Association matrices and LD row maps --------

#' Export a GWAS association matrix
#'
#' Writes little-endian IEEE FP64 values in column-major order to `z.f64.bin`
#' and creates `markers.tsv` and `studies.tsv` in the same directory.
#' Rows follow `GAlist$rsids`. Columns follow the full study inventory in
#' GAlist order. The inventory must contain at least two studies.
#' The writer aligns study scores to the gact marker catalog. It computes Z
#' from finite beta and positive finite SE when available, with stored Z as a
#' fallback. Missing and unusable scores are stored as NaN.
#' All studies use the gact catalog effect allele. Reference panels must use
#' that same dosage allele. Each panel's binary map indexes Z rows in LD order.
#' For EUR, EAS and SAS studies, the writer looks for sparse LD with prefix
#' `g1000_<ancestry>.subset.aligned` in the marker directory, using lowercase
#' ancestry codes. Use [addLDPanels()] for custom paths.
#' With LD regions installed by `downloadDB(GAlist, what = "ldetect")`, the
#' export also writes zero-based Z and LD row pairs for each matching region.
#' @param GAlist A gact database descriptor with marker and GWAS data installed.
#' @param out_dir A new output directory. Defaults to `matrix` inside the
#'   database's `gstat` directory. Its parent must already exist.
#' @param stat Statistic to export. The default, `"z"`, exports Z scores and
#'   is the only supported value.
#' @param progress Show progress while writing the marker catalog, LD maps
#'   and study columns. The default is `TRUE`.
#' @return Invisibly, a list with the matrix and catalog paths, matrix dimensions,
#'   any LD row map paths and the LD region catalog path when present.
#' @export
writeMarkerStat <- function(GAlist, out_dir = NULL, stat = "z", progress = TRUE) {
  if (!identical(stat, "z"))
    stop("stat must be \"z\".", call. = FALSE)
  if (!identical(progress, TRUE) && !identical(progress, FALSE))
    stop("progress must be TRUE or FALSE.", call. = FALSE)

  valid_ids <- function(x) is.character(x) && length(x) > 0L &&
    !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
  rsids <- GAlist$rsids
  if (!is.character(rsids) || !length(rsids))
    stop("GAlist$rsids must contain the marker backbone.", call. = FALSE)
  studies <- GAlist$studies
  if (is.null(studies$id)) studies <- GAlist$study
  ids <- as.character(studies$id)
  if (length(ids) < 2L)
    stop("The association matrix requires at least two studies.", call. = FALSE)
  if (!valid_ids(ids))
    stop("GAlist must contain unique, nonempty study IDs.", call. = FALSE)
  ancestry <- as.character(studies$ancestry)
  if (!length(ancestry)) ancestry <- rep(NA_character_, length(ids))
  if (length(ancestry) != length(ids))
    stop("The study ancestry column must match the study IDs.", call. = FALSE)
  ancestry <- .gact_ancestry_codes(ancestry)
  ancestry[is.na(ancestry) | !nzchar(ancestry)] <- "unknown"
  if (any(grepl("[\t\r\n]", c(ids, ancestry))))
    stop("Study IDs and ancestries must be single TSV fields.", call. = FALSE)

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

  GAlist$ldPanels <- .find_aligned_ld(marker_dir, ancestry, GAlist$ldPanels)
  panels <- as.character(names(GAlist$ldPanels))
  panel_files <- vapply(GAlist$ldPanels, function(panel) panel$markers, character(1L))
  if (any(grepl("[\t\r\n]", panel_files)))
    stop("LD marker paths must be single TSV fields.", call. = FALSE)

  if (is.null(out_dir)) {
    gstat_dir <- if ("gstat" %in% names(GAlist$dirs)) GAlist$dirs[["gstat"]]
    if (!is.character(gstat_dir) || length(gstat_dir) != 1L ||
        is.na(gstat_dir) || !nzchar(gstat_dir))
      stop("Supply a GAlist gstat directory or out_dir.", call. = FALSE)
    out_dir <- file.path(gstat_dir, "matrix")
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
  region_files <- .stage_ld_regions(GAlist, stage)
  assembly <- if (any(nzchar(region_files))) GAlist$assembly else ""
  without_ld <- !ancestry %in% panels
  if (progress && any(without_ld))
    message("[writeMarkerStat] No matching LD for ",
            paste(paste0(ids[without_ld], " (", ancestry[without_ld], ")"), collapse = ", "), ".")
  .Call(C_gact_write_z, enc2utf8(rsids),
        enc2utf8(normalizePath(marker_file, winslash = "/", mustWork = TRUE)),
        enc2utf8(normalizePath(files, winslash = "/", mustWork = TRUE)),
        enc2utf8(ids), enc2utf8(ancestry), enc2utf8(stage),
        enc2utf8(panels), enc2utf8(panel_files), enc2utf8(region_files), enc2utf8(assembly),
        progress)
  if (file.exists(out_dir) || !file.rename(stage, out_dir))
    stop("Cannot publish the Z matrix directory: ", out_dir, call. = FALSE)
  result <- list(z = file.path(out_dir, "z.f64.bin"),
    markers = file.path(out_dir, "markers.tsv"), studies = file.path(out_dir, "studies.tsv"),
    dim = c(length(rsids), length(ids)))
  if (length(panels)) {
    result$ld <- file.path(out_dir, "ld.tsv")
    result$ld_maps <- setNames(file.path(out_dir, paste0("ld-", panels, ".u32le.bin")), panels)
  }
  if (nzchar(assembly)) result$ld_regions <- file.path(out_dir, "ld_regions.tsv")
  if (progress) message("[writeMarkerStat] Matrix written to ", out_dir)
  invisible(result)
}
