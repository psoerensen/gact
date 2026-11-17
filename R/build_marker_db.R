#' Build the GRCh38 human GWAS database marker catalog
#'
#' Publishes a new GRCh38 database directory with a marker catalog, a
#' chromosome-position-allele index, provenance and a saved GAlist. The
#' catalog uses ALT as its effect allele and the frequency of ALT as EAF.
#' Rows at repeated chromosome-position pairs or with repeated rsIDs are
#' omitted. All 22 autosomes must have retained markers. This version
#' provides resources for GWAS ingestion and retrieval.
#' The build reports progress as it verifies and reads each chromosome.
#'
#' @param input_dir Directory containing the checksum-verified `chr1.tsv.gz`
#'   to `chr22.tsv.gz` tables and their `chr1.json` to `chr22.json` receipts
#'   from [prepareIGSRMarkers()].
#' @param dbdir Parent directory for the new database directory.
#'   The build stops if that database directory already exists.
#' @return Invisibly, the saved GAlist.
#' @export
buildMarkerDB <- function(input_dir, dbdir = getwd()) {
  version <- "hsa.GRCh38.0.0.1"
  dbdir <- normalizePath(dbdir, winslash = "/", mustWork = TRUE)
  target <- file.path(dbdir, version)
  if (file.exists(target)) stop("Database already exists: ", target)
  autosomes <- seq_len(22L)
  files <- file.path(input_dir, paste0("chr", autosomes, ".tsv.gz"))
  receipts <- file.path(input_dir, paste0("chr", autosomes, ".json"))
  if (!all(file.exists(c(files, receipts)))) {
    stop("Supply prepared marker tables and receipts for all 22 autosomes.")
  }
  message("Preparing the GRCh38 marker backbone from ", input_dir)
  provenance <- lapply(receipts, jsonlite::read_json)
  expected <- c("rsids", "chr", "pos", "ref", "alt", "af")
  parts <- lapply(autosomes, function(i) {
    message("Chromosome ", i, "/", length(autosomes), ": checking and reading ", basename(files[i]))
    receipt <- provenance[[i]]
    checksum <- unname(tools::md5sum(files[i]))
    if (!identical(receipt$build, "GRCh38") ||
        !isTRUE(receipt$chromosome == i) ||
        !identical(receipt$output$file, basename(files[i])) ||
        is.na(checksum) || !identical(receipt$output$md5, checksum)) {
      stop("Invalid GRCh38 receipt or checksum for chromosome ", i)
    }
    markers <- data.table::fread(files[i],
      colClasses = c("character", "integer", "integer", "character", "character", "numeric"))
    if (!identical(names(markers), expected) ||
        !nrow(markers) || anyNA(markers) ||
        !isTRUE(receipt$counts$retained == nrow(markers)) ||
        any(markers$chr != i | markers$pos < 1L |
            markers$ref == markers$alt | !markers$ref %in% c("A", "C", "G", "T") |
            !markers$alt %in% c("A", "C", "G", "T") |
            !is.finite(markers$af) | markers$af <= 0 | markers$af >= 1 |
            !grepl("^rs[0-9]+$", markers$rsids))) {
      stop("Invalid GRCh38 marker data for chromosome ", i)
    }
    markers
  })
  message("Combining the chromosome tables and removing duplicate markers.")
  markers <- data.table::rbindlist(parts)
  rm(parts)
  n_input <- nrow(markers)
  position_cols <- c("chr", "pos")
  repeated_position <- duplicated(markers, by = position_cols) |
    duplicated(markers, by = position_cols, fromLast = TRUE)
  repeated_rsid <- duplicated(markers$rsids) | duplicated(markers$rsids, fromLast = TRUE)
  markers <- markers[!(repeated_position | repeated_rsid)]
  if (!setequal(markers$chr, autosomes)) stop("Catalog must cover all 22 autosomes.")
  data.table::setorderv(markers, position_cols)
  data.table::setnames(markers, c("ref", "alt", "af"), c("nea", "ea", "eaf"))
  data.table::setcolorder(markers, c("rsids", "chr", "pos", "ea", "nea", "eaf"))

  message("Writing the marker catalog and coordinate index.")
  stage <- tempfile(".gact-marker-db-", tmpdir = dbdir)
  if (!dir.create(stage, mode = "0700")) stop("Cannot create marker database staging directory.")
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  ga <- createDB(version, dbdir = stage)
  marker_dir <- ga$dirs[["marker"]]
  catalog <- file.path(marker_dir, "markers.txt")
  index <- file.path(marker_dir, "cpra.txt")
  data.table::fwrite(markers, catalog)
  cpra <- with(markers, paste(chr, pos, ea, nea, sep = "_"))
  data.table::fwrite(data.frame(cpra), index, col.names = FALSE)
  for (path in c(catalog, index)) {
    R.utils::gzip(path, destname = paste0(path, ".gz"), remove = TRUE)
  }
  ga$rsids <- markers$rsids
  ga$dbdir <- target
  ga$dirs[] <- file.path(target, names(ga$dirs))
  files <- paste0(c(catalog, index), ".gz")
  checksums <- stats::setNames(as.list(unname(tools::md5sum(files))), basename(files))
  ga$metadata$catalog <- list(build = "GRCh38", effect_allele = "ALT",
    frequency = "IGSR global alternate allele count divided by allele number",
    autosomes = autosomes, n_input = n_input, n_markers = nrow(markers),
    n_repeated_position = sum(repeated_position),
    n_repeated_rsid = sum(repeated_rsid),
    n_dropped = n_input - nrow(markers), checksums = checksums)
  jsonlite::write_json(c(ga$metadata$catalog, list(inputs = provenance)),
    file.path(marker_dir, "marker_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
  stage_db <- file.path(stage, version)
  saveRDS(ga, file.path(stage_db, paste0("GAlist_", version, ".rds")))
  if (file.exists(target) || !file.rename(stage_db, target)) {
    stop("Cannot publish marker database at ", target)
  }
  invisible(ga)
}
