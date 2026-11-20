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
  ancestries <- setdiff(intersect(c("EUR", "EAS", "SAS", "AFR", "AMR"), study_ancestry), names(panels))
  suffixes <- c(".markers.tsv", ".meta.txt", ".row_ptr.u64.bin",
                ".col_idx.u32.0based.bin", ".values.f32.bin")
  for (ancestry in ancestries) {
    stem <- paste0("g1000_", tolower(ancestry), ".subset.aligned")
    default <- paste0("g1000_", tolower(ancestry), ".1cm.aligned")
    if (any(file.exists(file.path(marker_dir, paste0(default, suffixes))))) stem <- default
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

# LD score preparation and retrieval --------

.gact_register_prepared_ld <- function(GAlist, ancestry, directory, metadata) {
  GAlist$ldPanels[[ancestry]] <- list(
    markers = normalizePath(file.path(directory, metadata$sparse$markers), winslash = "/", mustWork = TRUE),
    assembly = metadata$build,
    descriptor = file.path(directory, metadata$sparse$descriptor))
  GAlist
}

# Key file stamps by receipt name so a moved file or a missing stamp cannot
# accidentally take a different file's cached checksum.
.gact_ld_stamps <- function(files) {
  info <- file.info(unname(files))
  if (anyNA(info$size) || any(info$isdir)) stop("Missing LD file for checksum.")
  stats::setNames(sprintf("%.0f/%.17g", info$size, as.numeric(info$mtime)), names(files))
}

.gact_ld_hashes <- function(files, stamps, previous = NULL, previous_stamps = NULL) {
  keys <- names(files)
  hashes <- as.character(unlist(previous)[keys])
  cached <- unlist(previous_stamps)[keys]
  if (!length(hashes)) hashes <- rep(NA_character_, length(keys))
  reuse <- rep(FALSE, length(keys))
  if (length(cached)) reuse <- !is.na(cached) & cached == stamps & !is.na(hashes)
  if (any(!reuse)) hashes[!reuse] <- unname(tools::md5sum(files[!reuse]))
  if (anyNA(hashes)) stop("Cannot checksum LD files.")
  as.list(stats::setNames(hashes, keys))
}

#' Prepare signed sparse LD and ordinary LD scores from an aligned BED panel
#'
#' Uses the installed ancestry panel in physical BED order. Each pair within
#' the inclusive BIM genetic-distance window on one chromosome contributes
#' when its raw squared correlation is at least `r2`. The score is one plus
#' the sum of retained off-diagonal squared correlations, without a finite-sample
#' correction. One sparse gprep call computes both outputs per ancestry.
#' @param GAlist Database descriptor pointing to aligned ancestry BED panels.
#'   BIM A1 must be the database effect allele. GRCh38 NYGC panels built or
#'   downloaded by gact already use this convention.
#' @param ancestry One or more of EUR, AFR, EAS, SAS or AMR.
#' @param nthreads Number of computation threads.
#' @param window_cm Positive maximum distance in centimorgans, default 1.
#'   Requires a populated, nondecreasing genetic map within each chromosome.
#' @param r2 Inclusive minimum squared correlation, from 0 through 1.
#' @param rsids Optional unique marker IDs. Both ends of a pair must belong
#'   to the selection. Unknown IDs raise an error. NULL uses the entire BED.
#' @param overwrite Recompute an existing result with these settings.
#' @param progress Report preparation, computation and reuse.
#' @param verify Rehash all files on reuse. Otherwise reuse recorded checksums
#'   for files with unchanged size and modification time.
#' @return Updated GAlist with sparse resources in `ldPanels` and matching
#'   scores in `ldscorePanels`. Save the descriptor to persist that selection.
#' @details Each selection and combination of settings has its own output files
#'   under marker/. Unchanged results reuse their checksum receipt. The indexed
#'   BIM/FAM metadata is shared across selections of the same panel.
#' @seealso [getLDscores()], [writeMarkerStat()]
#' @export
prepareLD <- function(GAlist, ancestry = "EUR", nthreads = 1L, window_cm = 1,
                      r2 = 0, rsids = NULL, overwrite = FALSE,
                      progress = TRUE, verify = FALSE) {
  build <- .gact_catalog_build(GAlist)
  if (length(build) != 1L || !build %in% c("GRCh37", "GRCh38"))
    stop("GAlist must specify GRCh37 or GRCh38.")
  directory <- GAlist$dirs[["marker"]]
  if (length(directory) != 1L || !dir.exists(directory))
    stop("GAlist must point to an existing marker directory.")
  ancestry <- toupper(ancestry)
  if (!length(ancestry) || anyNA(ancestry) || anyDuplicated(ancestry) ||
      any(!ancestry %in% c("EUR", "AFR", "EAS", "SAS", "AMR")))
    stop("Choose unique ancestries from EUR, AFR, EAS, SAS and AMR.")
  if (length(nthreads) != 1L || !is.numeric(nthreads) || is.na(nthreads) ||
      !is.finite(nthreads) || nthreads < 1 || nthreads > .Machine$integer.max || nthreads != floor(nthreads))
    stop("nthreads must be a positive integer.")
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite) ||
      !is.logical(progress) || length(progress) != 1L || is.na(progress))
    stop("overwrite and progress must be TRUE or FALSE.")
  if (!is.logical(verify) || length(verify) != 1L || is.na(verify)) stop("verify must be TRUE or FALSE.")
  if (!is.numeric(window_cm) || length(window_cm) != 1L || !is.finite(window_cm) || window_cm <= 0)
    stop("window_cm must be a positive finite distance in centimorgans.")
  if (!is.numeric(r2) || length(r2) != 1L || !is.finite(r2) || r2 < 0 || r2 > 1)
    stop("r2 must be a finite number between 0 and 1.")
  if (!is.null(rsids) && (!is.character(rsids) || !length(rsids) || anyNA(rsids) ||
      any(!nzchar(rsids)) || anyDuplicated(rsids) || any(grepl("[[:space:]]", rsids))))
    stop("rsids must contain unique, nonempty marker IDs.")
  selection_md5 <- if (is.null(rsids)) "all" else
    as.character(openssl::md5(paste(sort(rsids, method = "radix"), collapse = "\n")))
  settings <- list(window_cm = as.double(window_cm), min_r2 = as.double(r2), corrected = FALSE)
  tag <- "1cm"
  if (window_cm != 1 || r2 != 0 || !is.null(rsids)) {
    key <- paste(sprintf("%.17g", c(window_cm, r2)), selection_md5, collapse = "/")
    tag <- paste0("ld-", as.character(openssl::md5(key)))
  }
  for (population in ancestry) {
    prefix <- file.path(directory, paste0("g1000_", tolower(population)))
    inputs <- c(bed = paste0(prefix, ".bed"),
                bim = paste0(prefix, ".bim"), fam = paste0(prefix, ".fam"))
    if (any(!file.exists(inputs))) stop("Missing reference inputs: ",
      paste(basename(inputs[!file.exists(inputs)]), collapse = ", "))
    output <- file.path(directory, paste0("markers_1000G_", tolower(population), "_filtered.", tag, ".txt.gz"))
    receipt <- paste0(output, ".json")
    ld_prefix <- file.path(directory, paste0("g1000_", tolower(population), ".", tag, ".aligned"))
    descriptor <- file.path(directory, paste0("Glist_1000G_", tolower(population), "_filtered.", tag, ".rds"))
    lock <- paste0(output, ".lock")
    if (!dir.create(lock, showWarnings = FALSE)) stop("LD-score preparation is locked: ", lock)
    tryCatch({
      old <- if (!overwrite && file.exists(receipt)) jsonlite::read_json(receipt)
      panel_directory <- normalizePath(directory, winslash = "/")
      cached <- if (!verify && identical(old$sparse$directory, panel_directory))
        old$file_stamps
      input_stamps <- .gact_ld_stamps(inputs)
      hashes <- .gact_ld_hashes(inputs, input_stamps, old$inputs, cached$inputs)
      if (file.exists(output) || file.exists(receipt)) {
        if (!overwrite) {
          if (!file.exists(receipt)) stop("Existing scores have no provenance; use overwrite = TRUE to rebuild.")
          if (!identical(old$inputs, hashes) || !identical(old$selection_md5, selection_md5) ||
              !isTRUE(all.equal(old$settings, settings, tolerance = 0)) ||
              !identical(old$build, build) || !identical(old$ancestry, population))
            stop("LD-score inputs or settings changed; use overwrite = TRUE to rebuild.")
          score_file <- c(scores = output)
          output_stamps <- .gact_ld_stamps(score_file)
          actual <- .gact_ld_hashes(score_file, output_stamps,
            list(scores = old$md5), cached$output)
          if (!identical(actual$scores, old$md5)) stop("LD-score checksum differs from its provenance manifest.")
          if (is.null(old$sparse)) stop("Missing sparse LD receipt.")
          keys <- names(old$sparse$checksums)
          files <- stats::setNames(file.path(directory, keys), keys)
          if (!length(files) || any(!file.exists(files))) stop("Missing sparse LD files for checksum.")
          sparse_stamps <- .gact_ld_stamps(files)
          if (!identical(.gact_ld_hashes(files, sparse_stamps,
                old$sparse$checksums, cached$sparse), old$sparse$checksums))
            stop("Sparse LD checksum differs or files are missing; use overwrite = TRUE to rebuild.")
          if (is.null(cached$sparse)) {
            saved <- readRDS(descriptor)
            gbase::validate_Glist(saved)
            if (!identical(normalizePath(dirname(saved$ldfiles$markers), winslash = "/"),
                           panel_directory))
              stop("Sparse LD descriptor points to a different marker directory.")
            rm(saved)
          }
          verified_stamps <- list(inputs = as.list(input_stamps),
            output = as.list(output_stamps), sparse = as.list(sparse_stamps))
          if (!identical(old$file_stamps, verified_stamps) ||
              !identical(old$sparse$directory, panel_directory)) {
            old$file_stamps <- verified_stamps
            old$sparse$directory <- panel_directory
            staged_receipt <- file.path(lock, "verified.json")
            jsonlite::write_json(old, staged_receipt, auto_unbox = TRUE, pretty = TRUE, digits = NA)
            if (!file.rename(staged_receipt, receipt)) stop("Cannot update verified LD receipt.")
          }
          GAlist <- .gact_register_prepared_ld(GAlist, population, directory, old)
          if (progress) message(population, ": reusing verified sparse LD and scores")
          GAlist$ldscorePanels[[population]] <- list(file = output, manifest = receipt, build = build)
          next
        }
      }
      if (progress) message(population, ": opening aligned BED panel")
      if (!overwrite && (file.exists(descriptor) ||
          any(file.exists(paste0(ld_prefix, c(".markers.tsv", ".meta.txt", ".row_ptr.u64.bin",
                                              ".col_idx.u32.0based.bin", ".values.f32.bin"))))))
        stop("Unverified sparse LD outputs already exist; use overwrite = TRUE to rebuild.")
      # Reuse the existing gsuite indexed reader. Its metadata columns stay
      # disk-backed. The BED already uses the database effect allele.
      # Content-address the immutable index so changed panels cannot reuse it.
      index_key <- as.character(openssl::md5(paste(unlist(hashes[c("bed", "bim", "fam")]), collapse = "")))
      index_dir <- file.path(directory, paste0("g1000_", tolower(population), ".metadata-", index_key))
      if (dir.exists(index_dir)) {
        if (progress) message(population, ": reopening indexed reference metadata")
        gl <- gbase::open_mapped_Glist(index_dir)
      } else {
        if (progress) message(population, ": building reference metadata index once")
        gl <- gbase::gprep_mapped(bedfiles = inputs[["bed"]],
          bimfiles = inputs[["bim"]], famfiles = inputs[["fam"]],
          out_dir = index_dir, study = paste("1000G", population))
      }
      if (progress) message(population, ": computing sparse LD and scores (", window_cm,
                            " cM, minimum r^2 ", r2, ")")
      if (overwrite) invisible(gc())
      gl <- gbase::gprep(gl, task = "sparseld", rsids = rsids, window_cm = window_cm,
        r2 = r2, out_prefix = ld_prefix, overwrite = overwrite, nthreads = nthreads,
        progress = progress)
      gl$assembly <- build
      ids <- gl$rsidsLD[[1L]]
      scores <- unname(gl$ldscores[[1L]])
      if (length(ids) != length(scores) || any(!is.finite(scores)) || any(scores < 1))
        stop("Invalid LD-score result from gbase.")
      stage <- tempfile(".ldscores-", tmpdir = lock)
      metadata <- list(build = build, ancestry = population, settings = settings,
        inputs = hashes, selection_md5 = selection_md5,
        n_markers = length(ids), n_samples = length(gl$idsLD),
        gbase_version = as.character(utils::packageVersion("gbase")),
        formula = "1 + sum of retained in-window off-diagonal r^2; uncorrected")
      if (progress) message(population, ": saving sparse LD descriptor and receipt checksums")
      saveRDS(gl, paste0(stage, ".rds"), compress = FALSE)
      if (!file.copy(paste0(stage, ".rds"), descriptor, overwrite = overwrite))
        stop("Cannot save the sparse LD descriptor.")
      sidecars <- vapply(gl$sparseLD$metadata_storage, function(x)
        if (is.null(x)) "" else x$path, character(1L))
      index_files <- list.files(index_dir, recursive = TRUE, full.names = TRUE)
      files <- c(unlist(gl$ldfiles, use.names = FALSE), sidecars[nzchar(sidecars)], descriptor, index_files)
      relative <- substring(normalizePath(files, winslash = "/", mustWork = TRUE),
                            nchar(normalizePath(directory, winslash = "/", mustWork = TRUE)) + 2L)
      # gbase already recorded all five LD file hashes. Its mapped writer
      # streams the three CSR hashes. Preserve the existing receipt keys.
      known <- c(unname(gl$sparseLD$checksums[names(gl$ldfiles)]),
        vapply(Filter(Negate(is.null), gl$sparseLD$metadata_storage), "[[", "", "md5"))
      remaining <- c(descriptor, index_files)
      checksums <- as.list(setNames(c(known, unname(tools::md5sum(remaining))), relative))
      sparse_stamps <- .gact_ld_stamps(setNames(files, relative))
      if (anyNA(unlist(checksums))) stop("Cannot checksum sparse LD outputs.")
      metadata$sparse <- list(directory = panel_directory,
        descriptor = basename(descriptor),
        markers = basename(gl$ldfiles$markers), checksums = checksums)
      marker_file <- gl$ldfiles$markers
      rm(gl)
      invisible(gc())
      final_stamps <- .gact_ld_stamps(inputs)
      if (!identical(input_stamps, final_stamps) &&
          !identical(.gact_ld_hashes(inputs, final_stamps, hashes, input_stamps), hashes))
        stop("Reference inputs changed during LD preparation.")
      if (progress) message(population, ": writing compressed LD-score table")
      fields <- c("marker_id", "chromosome", "position_bp", "allele1", "allele2")
      markers <- data.table::fread(marker_file, select = fields)
      data.table::setnames(markers, fields, c("rsids", "chr", "pos", "ea", "nea"))
      if (!identical(markers$rsids, ids)) stop("Sparse LD sidecar and score marker order differ.")
      data.table::set(markers, j = "chr", value = sub("^chr", "", markers$chr, ignore.case = TRUE))
      data.table::set(markers, j = "ldscores", value = scores)
      data.table::fwrite(markers, stage, sep = "\t", compress = "none")
      R.utils::gzip(stage, destname = paste0(stage, ".gz"), remove = TRUE)
      metadata$md5 <- unname(tools::md5sum(paste0(stage, ".gz")))
      if (!file.copy(paste0(stage, ".gz"), output, overwrite = overwrite))
        stop("Cannot publish LD scores.")
      metadata$file_stamps <- list(inputs = as.list(final_stamps),
        output = as.list(.gact_ld_stamps(c(scores = output))), sparse = as.list(sparse_stamps))
      jsonlite::write_json(metadata, paste0(stage, ".json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
      if (!file.copy(paste0(stage, ".json"), receipt, overwrite = overwrite))
        stop("Cannot publish LD scores and provenance.")
      GAlist$ldscorePanels[[population]] <- list(file = output, manifest = receipt, build = build)
      GAlist <- .gact_register_prepared_ld(GAlist, population, directory, metadata)
      if (progress) message(population, ": saved ", length(ids), " LD scores")
      rm(ids, scores, markers)
      invisible(gc())
    }, finally = unlink(lock, recursive = TRUE))
  }
  GAlist
}

#' Prepare the local HapMap3 reference for ordinary LDSC
#'
#' @param GAlist Human database with an installed aligned 1000G ancestry panel and genetic map.
#' @param ancestry One of EUR, AFR, EAS, SAS or AMR.
#' @param nthreads Number of genotype computation threads.
#' @param hm3_rsids Optional HapMap3 marker IDs. By default use the rsIDs from
#'   the checksum-pinned Pan-UKB EUR HapMap3 resource matched to this release.
#' @param overwrite Recompute an existing reference.
#' @return GAlist with the reference path under `ldscPanels`.
#' @details For each HapMap3 regression SNP, sum finite-sample corrected squared
#'   correlations within 1 cM over variable panel markers with MAF above 1%.
#'   A second sum over regression SNPs supplies the LD weights. Exclude the
#'   extended MHC (chr6, 25-34 Mb inclusive) from regression and weight SNPs.
#'   The common-marker count is the number of reference SNPs with MAF above 5%.
#' @export
prepareLDSC <- function(GAlist, ancestry = "EUR", nthreads = 1L, hm3_rsids = NULL,
                        overwrite = FALSE) {
  ancestry <- match.arg(ancestry, c("EUR", "AFR", "EAS", "SAS", "AMR"))
  build <- .gact_catalog_build(GAlist)
  if (length(build) != 1L || !build %in% c("GRCh37", "GRCh38"))
    stop("prepareLDSC requires a human panel with a known build.")
  GAlist$metadata$build <- build
  if (length(nthreads) != 1L || !is.numeric(nthreads) || !is.finite(nthreads) ||
      nthreads < 1 || nthreads != floor(nthreads)) stop("Expected a positive thread count.")
  directory <- GAlist$dirs[["marker"]]
  prefix <- file.path(directory, paste0("g1000_", tolower(ancestry)))
  files <- stats::setNames(paste0(prefix, c(".bed", ".bim", ".fam")), c("bed", "bim", "fam"))
  stamps <- .gact_ld_stamps(files)
  source <- "caller"
  if (is.null(hm3_rsids)) {
    GAlist <- .gact_install_panukb_scores(GAlist, "EUR")
    hm3_rsids <- data.table::fread(.gact_panukb_files(directory, "EUR"), select = "rsids")$rsids
    source <- "Pan-UKB EUR HapMap3 rsIDs"
  }
  if (!is.character(hm3_rsids) || !length(hm3_rsids) || anyNA(hm3_rsids) ||
      any(!nzchar(hm3_rsids)) || anyDuplicated(hm3_rsids)) stop("Expected unique HapMap3 rsIDs.")
  selection <- intersect(hm3_rsids, GAlist$rsids)
  if (length(selection) < 3L) stop("Too few HapMap3 SNPs match the catalog.")
  settings <- list(version = 1L, build = build, ancestry = ancestry, source = source,
    hm3_md5 = as.character(openssl::md5(paste(sort(selection, method = "radix"), collapse = "\n"))),
    window_cm = 1, min_maf = 0.01, count_maf = 0.05, corrected = TRUE)
  output <- file.path(directory, paste0("ldsc_1000g_", tolower(ancestry), "_hm3.rds"))
  lock <- paste0(output, ".lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Cannot acquire LDSC preparation lock: ", lock)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  GAlist$ldscPanels[[ancestry]] <- list(file = output, build = build)
  previous <- if (file.exists(output)) readRDS(output)
  hashes <- .gact_ld_hashes(files, stamps, previous$inputs, previous$stamps)
  if (!is.null(previous) && !overwrite) {
    if (!identical(previous$settings, settings) || !identical(previous$inputs, hashes))
      stop("LDSC reference inputs changed. Use overwrite = TRUE to recompute.")
    return(GAlist)
  }
  # Bound genotype metadata and scratch disk by the largest chromosome.
  chr <- sub("^chr", "", data.table::fread(files["bim"], header = FALSE, select = 1L)[[1L]])
  runs <- rle(chr)
  if (anyDuplicated(runs$values) || any(!runs$values %in% as.character(1:22)))
    stop("Expected contiguous autosomes in the aligned BIM.")
  n <- nrow(data.table::fread(files["fam"], header = FALSE, select = 1L))
  stride <- ceiling(n / 4)
  if (n < 3 || file.size(files["bed"]) != 3 + length(chr) * stride)
    stop("Invalid reference BED extent or sample count.")
  bed <- file(files["bed"], "rb")
  on.exit(close(bed), add = TRUE)
  magic <- readBin(bed, "raw", 3L)
  if (!identical(magic, as.raw(c(108, 27, 1)))) stop("Expected SNP-major PLINK BED.")
  bim <- file(files["bim"], "rt")
  on.exit(close(bim), add = TRUE)
  tmp <- file.path(lock, "chr")
  rows <- vector("list", length(runs$values))
  M <- 0
  for (i in seq_along(rows)) {
    chromosome <- runs$values[i]
    message(ancestry, ": preparing HapMap3 LDSC scores for chr", chromosome)
    len <- runs$lengths[i]
    writeLines(readLines(bim, n = len), paste0(tmp, ".bim"))
    out <- file(paste0(tmp, ".bed"), "wb")
    remaining <- len * stride
    tryCatch({
      writeBin(magic, out)
      while (remaining > 0) {
        bytes <- readBin(bed, "raw", min(4194304, remaining))
        if (!length(bytes)) stop("Truncated BED.")
        writeBin(bytes, out)
        remaining <- remaining - length(bytes)
      }
    }, finally = close(out))
    args <- list(bedfiles = paste0(tmp, ".bed"), famfiles = files[["fam"]], nthreads = nthreads)
    if ("companion" %in% names(formals(gbase::gprep))) args["companion"] <- list(NULL)
    gl <- do.call(gbase::gprep, args)
    cm <- gl$map[[1L]]
    if (any(!is.finite(cm)) || is.unsorted(cm) || (len > 1 && cm[1L] == cm[len]))
      stop("Expected a populated, ordered genetic map on chromosome ", chromosome, ".")
    # Integer allele counts make strict 1% and 5% boundaries strand-independent.
    copies <- 2 * (gl$n0[[1L]] + gl$n1[[1L]] + gl$n2[[1L]])
    minor <- 2 * pmin(gl$n0[[1L]], gl$n2[[1L]]) + gl$n1[[1L]]
    variable <- gl$het[[1L]] < 1 & 100 * minor > copies
    M <- M + sum(variable & 20 * minor > copies)
    ids <- gl$rsids[[1L]][variable]
    pos <- gl$pos[[1L]][variable]
    hm <- ids %in% selection & !(chromosome == "6" & pos >= 25e6 & pos <= 34e6)
    if (!any(hm)) next
    annotation <- cbind(reference = rep(1, length(ids)), hapmap3 = as.numeric(hm))
    rownames(annotation) <- ids
    gl <- gbase::gprep(gl, task = "ldscores", rsids = ids, annotation = annotation,
      powers = c(U = -1), corrected = TRUE, r2 = 0, window_cm = 1, nthreads = nthreads)
    scores <- gl$ldscore_models$scores
    rows[[i]] <- data.frame(rsids = ids[hm], chr = as.integer(chromosome),
      ld_score = unname(scores[hm, "reference", "U"]),
      weight_ld = unname(scores[hm, "hapmap3", "U"]))
  }
  rows <- do.call(rbind, rows)
  if (is.null(rows) || nrow(rows) < 3 || M < 1 || anyDuplicated(rows$rsids) ||
      any(!is.finite(rows$ld_score) | !is.finite(rows$weight_ld))) stop("Invalid LDSC reference.")
  if (!identical(stamps, .gact_ld_stamps(files))) stop("Panel changed during LDSC preparation.")
  staged <- file.path(lock, "scores.rds")
  saveRDS(list(settings = settings, inputs = hashes, stamps = stamps, n_samples = n,
    marker_count = M, rows = rows), staged)
  if (!file.rename(staged, output)) stop("Cannot publish LDSC reference.")
  message(ancestry, ": saved ", nrow(rows), " HapMap3 scores; M_5_50 = ", M)
  GAlist
}

# Pan-UK Biobank's published, ancestry-specific HapMap3 LD scores.
.gact_panukb_panel <- function(ancestry) {
  ancestry <- sub("^PanUKBB-", "", ancestry)
  manifest <- system.file("extdata", "panukb-ldscores.tsv", package = "gact", mustWork = TRUE)
  panels <- utils::read.delim(manifest, stringsAsFactors = FALSE)
  if (length(ancestry) != 1L || !ancestry %in% panels$ancestry)
    stop("Choose one Pan-UKB ancestry: ", paste(panels$ancestry, collapse = ", "), ".")
  panels[panels$ancestry == ancestry, , drop = FALSE]
}

.gact_panukb_files <- function(directory, ancestry) {
  file.path(directory, paste0("markers_panukb_", tolower(ancestry), "_hm3.txt.gz"))
}

.gact_panukb_table <- function(path, catalog) {
  scores <- data.table::fread(path, select = c("SNP", "L2"))
  if (!identical(names(scores), c("SNP", "L2")) || !is.numeric(scores$L2))
    stop("Invalid Pan-UKB score columns.")
  repeated <- duplicated(scores$SNP) | duplicated(scores$SNP, fromLast = TRUE)
  valid <- grepl("^rs[0-9]+$", scores$SNP) & is.finite(scores$L2)
  counts <- list(input = nrow(scores), invalid = sum(!valid),
                 ambiguous_source = sum(valid & repeated))
  scores <- scores[valid & !repeated]
  markers <- data.table::fread(catalog, select = c("rsids", "chr", "pos"))
  if (!identical(names(markers), c("rsids", "chr", "pos"))) stop("Invalid marker catalog.")
  markers <- markers[markers$rsids %in% scores$SNP]
  at <- match(scores$SNP, markers$rsids)
  ambiguous <- scores$SNP %in% markers$rsids[duplicated(markers$rsids)]
  counts$outside_backbone <- sum(is.na(at))
  counts$ambiguous_backbone <- sum(ambiguous)
  keep <- !is.na(at) & !ambiguous
  if (!any(keep)) stop("No unambiguous Pan-UKB rsIDs match the marker catalog.")
  rows <- markers[at[keep]]
  chr <- suppressWarnings(as.integer(sub("^chr", "", rows$chr, ignore.case = TRUE)))
  if (anyNA(chr) || any(!chr %in% 1:22) || anyNA(rows$pos) || any(rows$pos <= 0))
    stop("Invalid autosomal coordinates in matched marker catalog.")
  data.table::set(rows, j = "chr", value = chr)
  data.table::set(rows, j = "ldscores", value = scores$L2[keep])
  data.table::setorder(rows, chr, pos, rsids)
  counts$matched <- nrow(rows)
  list(rows = rows, counts = counts)
}

.gact_install_panukb_scores <- function(GAlist, ancestry) {
  spec <- .gact_panukb_panel(ancestry)
  ancestry <- spec$ancestry
  build <- .gact_catalog_build(GAlist)
  if (is.null(build) || !build %in% c("GRCh37", "GRCh38")) stop("Expected a human marker catalog.")
  directory <- GAlist$dirs[["marker"]]
  catalog <- file.path(directory, "markers.txt.gz")
  if (!file.exists(catalog)) stop("Install the release marker catalog first.")
  output <- .gact_panukb_files(directory, ancestry)
  receipt <- paste0(output, ".json")
  names <- paste0("UKBB.", ancestry, c(".rsid.l2.ldscore.gz", ".l2.M_5_50"))
  urls <- paste0("https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/", names)
  expected <- c(spec$scores_sha256, spec$count_sha256)
  catalog_md5 <- unname(tools::md5sum(catalog))
  old <- if (file.exists(receipt)) jsonlite::read_json(receipt, simplifyVector = TRUE)
  reusable <- file.exists(output) && identical(old$source, "PanUKB") &&
    identical(old$build, build) && identical(old$ancestry, ancestry) &&
    identical(old$catalog_md5, catalog_md5) && identical(old$input_sha256, expected) &&
    identical(old$marker_count, spec$marker_count) && identical(old$md5, unname(tools::md5sum(output)))
  if (!reusable) {
    stage <- tempfile(".panukb-", tmpdir = directory)
    if (!dir.create(stage)) stop("Cannot stage Pan-UKB scores.")
    on.exit(unlink(stage, recursive = TRUE))
    inputs <- file.path(stage, names)
    for (i in seq_along(urls)) {
      message("Downloading ", names[i])
      curl::curl_download(urls[i], inputs[i], mode = "wb")
      if (.gact_sha256_file(inputs[i]) != expected[i]) stop("SHA-256 mismatch for ", names[i])
    }
    count <- scan(inputs[2L], what = integer(), quiet = TRUE)
    if (!identical(count, spec$marker_count)) stop("Invalid Pan-UKB M_5_50.")
    result <- .gact_panukb_table(inputs[1L], catalog)
    plain <- file.path(stage, "scores.tsv")
    data.table::fwrite(result$rows, plain, sep = "\t")
    compressed <- paste0(plain, ".gz")
    R.utils::gzip(plain, destname = compressed)
    metadata <- list(source = "PanUKB", ancestry = ancestry, build = build,
      urls = urls, input_sha256 = expected, catalog_md5 = catalog_md5,
      marker_count = count, counts = result$counts, md5 = unname(tools::md5sum(compressed)))
    staged_receipt <- file.path(stage, "receipt.json")
    jsonlite::write_json(metadata, staged_receipt, auto_unbox = TRUE, pretty = TRUE)
    if (!file.rename(compressed, output) || !file.rename(staged_receipt, receipt))
      stop("Could not publish Pan-UKB scores and receipt. Retry downloadDB().")
  }
  message(if (reusable) "Reusing " else "Installed ", basename(output))
  GAlist <- log_download(GAlist, "panukb-ldscores", "PanUKB", urls, directory, "success")
  entries <- tail(seq_len(nrow(GAlist$downloads)), 2L)
  GAlist$downloads$dest_file[entries] <- output
  GAlist$downloads$file_size_MB[entries] <- file.size(output) / 1e6
  GAlist
}

#' Retrieve Linkage Disequilibrium (LD) Scores for Specific Genetic Markers
#'
#' This function retrieves LD scores for specified genetic markers based on the population ancestry and version of the reference panel. The function reads in marker files and extracts LD scores, allowing for optional filtering by specific rsIDs.
#'
#' @param GAlist A list containing directory paths where marker data is stored. It is assumed that `GAlist$dirs["marker"]` points to the correct directory.
#' @param chr Optional chromosome filter.
#' @param ancestry Ancestry population: EUR, EAS, SAS, AFR or AMR. AFR and AMR
#'   support 1000G. Use `ancestry = "PanUKBB-EUR"` to select the Pan-UKB EUR scores;
#'   AFR, AMR, CSA, EAS and MID are also supported with this prefix.
#' @param version Choose "HapMap3" for the regression score column saved by
#'   prepareLDSC(), "1000G" (default) for the ordinary scores selected by
#'   prepareLD(), "Original" or "PanUKB". Legacy HapMap3 tables are also supported.
#' @param rsids An optional character vector of rsIDs to filter the LD scores. If not specified, all LD scores will be returned.
#'
#' @return A named vector of LD scores, where the names are rsIDs.
#'
#' @details
#' Use [prepareLD()] to prepare the `1000G` tables from matching BED
#' panels. Prepared tables are checked against their build, ancestry and
#' checksum manifest. GRCh38 tables require this provenance; existing GRCh37
#' tables without a manifest remain supported. No computation or download is
#' triggered by this reader.
#'
#' The function reads in precomputed LD score files for the selected ancestry and version, extracts the LD scores, and returns them as a named vector. If the `rsids` argument is provided, the returned vector will only include LD scores for the specified rsIDs.
#'
#' @examples
#' # Example usage:
#' GAlist <- list(dirs = c(marker = "/path/to/marker/files"))
#' ldscores <- getLDscores(GAlist, ancestry = "EUR", version = "HapMap3", rsids = c("rs123", "rs456"))
#'
#' @export
#'
getLDscores <- function(GAlist=NULL, chr=NULL, ancestry="EUR", version="1000G", rsids=NULL) {
 # Check if GAlist is provided
 if (is.null(GAlist) || !("marker" %in% names(GAlist$dirs))) {
  stop("Error: GAlist or GAlist$dirs['marker'] is missing.")
 }

 if (isTRUE(startsWith(as.character(ancestry), "PanUKBB-"))) {
   if (!missing(version) && !identical(version, "PanUKB"))
     stop("The PanUKBB- ancestry prefix conflicts with the requested LD-score version.")
   version <- "PanUKB"
 }
 if (identical(version, "HapMap3") && ancestry %in% c("EUR", "AFR", "EAS", "SAS", "AMR")) {
   path <- file.path(GAlist$dirs[["marker"]], paste0("ldsc_1000g_", tolower(ancestry), "_hm3.rds"))
   if (file.exists(path)) {
     ref <- readRDS(path)
     if (!identical(ref$settings$build, .gact_catalog_build(GAlist)) ||
         !identical(ref$settings$ancestry, ancestry)) stop("LDSC reference build or ancestry differs.")
     rows <- ref$rows
     if (!is.null(chr)) rows <- rows[rows$chr %in% chr, , drop = FALSE]
     if (!is.null(rsids)) rows <- rows[rows$rsids %in% rsids, , drop = FALSE]
     return(stats::setNames(rows$ld_score, rows$rsids))
   }
 }
 if (identical(version, "PanUKB")) {
   spec <- .gact_panukb_panel(ancestry)
   ancestry <- spec$ancestry
   marker_file_path <- .gact_panukb_files(GAlist$dirs[["marker"]], ancestry)
   if (!file.exists(paste0(marker_file_path, ".json")) || !file.exists(marker_file_path))
     stop("Install the panel with downloadDB(GAlist, what = 'panukb-ldscores', ancestry = '", ancestry, "').")
 } else {
   if (!ancestry %in% c("EUR", "EAS", "SAS", "AFR", "AMR"))
     stop("Error: Invalid ancestry. Choose from EUR, EAS, SAS, AFR or AMR.")
   if (!version %in% c("HapMap3", "1000G", "Original"))
     stop("Error: Invalid version. Choose HapMap3, 1000G, Original or PanUKB.")
   file_name <- switch(
     paste(ancestry, version, sep = "_"),
     "EUR_HapMap3" = "markers_1000G_eur_hm3.txt.gz",
     "EAS_HapMap3" = "markers_1000G_eas_hm3.txt.gz",
     "SAS_HapMap3" = "markers_1000G_sas_hm3.txt.gz",
     "EUR_1000G" = "markers_1000G_eur_filtered.txt.gz",
     "EAS_1000G" = "markers_1000G_eas_filtered.txt.gz",
     "SAS_1000G" = "markers_1000G_sas_filtered.txt.gz",
     "AFR_1000G" = "markers_1000G_afr_filtered.txt.gz",
     "AMR_1000G" = "markers_1000G_amr_filtered.txt.gz",
     "EUR_Original" = "markers_1000G_eur_w_ld.txt.txt.gz",
     stop("Error: No marker file found for the selected ancestry and version.")
    )
   marker_file_path <- file.path(GAlist$dirs["marker"], file_name)
 }
 if (version == "1000G") {
   registered <- GAlist$ldscorePanels[[ancestry]]$file
   default <- sub("[.]txt[.]gz$", ".1cm.txt.gz", marker_file_path)
   if (!is.null(registered)) marker_file_path <- registered
   else if (file.exists(default)) marker_file_path <- default
 }
 if (!file.exists(marker_file_path)) {
  stop("LD-score table is missing: ", marker_file_path,
       ". Run prepareLD(GAlist, ancestry = '", ancestry, "').")
 }

 receipt <- paste0(marker_file_path, ".json")
 build <- .gact_catalog_build(GAlist)
 if (file.exists(receipt)) {
  metadata <- jsonlite::read_json(receipt, simplifyVector = TRUE)
  if (identical(version, "PanUKB") &&
      (!identical(metadata$source, version) || !identical(metadata$marker_count, spec$marker_count) ||
       !identical(metadata$input_sha256, c(spec$scores_sha256, spec$count_sha256))))
    stop("Invalid Pan-UKB LD-score receipt.")
  if ((!is.null(build) && !identical(metadata$build, build)) ||
      !identical(metadata$ancestry, ancestry))
   stop("LD-score build or ancestry differs from the database.")
  if (!identical(metadata$md5, unname(tools::md5sum(marker_file_path))))
   stop("LD-score checksum differs from its provenance manifest.")
 } else if (identical(build, "GRCh38")) {
  stop("GRCh38 LD scores require provenance; run prepareLD().")
 }

 marker <- fread(marker_file_path, data.table = FALSE)
 # Check if 'ldscores' and 'rsids' columns exist
 if (!all(c("ldscores", "rsids") %in% colnames(marker))) {
  stop("Error: The marker file does not contain the required columns 'ldscores' and 'rsids'.")
 }
 if (anyNA(marker$rsids) || any(!nzchar(marker$rsids)) || anyDuplicated(marker$rsids) ||
     !is.numeric(marker$ldscores) || any(!is.finite(marker$ldscores)))
  stop("LD-score table must contain unique rsIDs and finite numeric scores.")
 if (!is.null(chr)) {
  if (!"chr" %in% names(marker)) stop("LD-score table has no chr column.")
  marker <- marker[marker$chr %in% chr, , drop = FALSE]
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
#' For Z-only imports, the writer reads the saved default-matrix column
#' and selects markers in the requested backbone order.
#' All studies use the gact catalog effect allele. Reference panels must use
#' that same dosage allele. Each panel's binary map indexes Z rows in LD order.
#' For EUR, EAS, SAS, AFR and AMR studies, the writer looks for sparse LD with prefix
#' `g1000_<ancestry>.1cm.aligned` (or the older `.subset.aligned` prefix) in the marker directory, using lowercase
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
  z_only <- which(studies$storage == "z")
  if (length(z_only)) {
    dir <- file.path(GAlist$dirs[["gstat"]], "matrix")
    stored_ids <- data.table::fread(file.path(dir, "studies.tsv"), select = "study_id")$study_id
    stored_markers <- data.table::fread(file.path(dir, "markers.tsv"), select = "rsids",
                                        colClasses = "character")$rsids
    rows <- match(rsids, stored_markers)
    cols <- match(ids[z_only], stored_ids)
    binary <- file.path(dir, "z.f64.bin")
    if (anyNA(rows) || anyNA(cols) || !file.exists(binary) ||
        file.size(binary) != 8 * length(stored_markers) * length(stored_ids))
      stop("Default Z matrix does not cover the selected markers and studies.")
    files[z_only] <- binary
    offsets <- rep(-1, length(files))
    offsets[z_only] <- (cols - 1) * 8 * length(stored_markers)
  }
  if (anyNA(files) || any(!file.exists(files) | dir.exists(files)))
    stop("Missing aligned study file for ", ids[which(
      is.na(files) | !file.exists(files) | dir.exists(files))[1L]], ".", call. = FALSE)
  files <- enc2utf8(normalizePath(files, winslash = "/", mustWork = TRUE))
  if (length(z_only)) {
    attr(files, "matrix_offsets") <- offsets
    attr(files, "matrix_rows") <- rows
  }
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
        files,
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
