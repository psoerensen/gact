#' Build whole-autosome NYGC GRCh38 genotype panels
#'
#' Build checksummed chromosome parts with [prepareNYGCGenotypes()] and merge
#' them into one BED/BIM/FAM triplet per ancestry. The staged parts are removed
#' when the build finishes. To retain parts for parallel processing or reuse,
#' prepare them separately and call [mergeNYGCGenotypes()].
#'
#' @inheritParams prepareNYGCGenotypes
#' @param out_dir New standalone output directory. Its parent must exist.
#' @param regions Optional data frame with integer-valued `chromosome`, `start`,
#'   and `end` columns and one inclusive interval per chromosome. Intended for
#'   bounded tests. The default builds all 22 autosomes.
#' @return Invisibly, a data frame with the file prefix and sample and SNP
#'   counts for each ancestry.
#' @export
buildNYGCGenotypes <- function(input_dir, marker_catalog, out_dir, regions = NULL) {
  out_dir <- .gact_nygc_new_directory(out_dir)
  intervals <- .gact_nygc_intervals(regions)
  parts <- tempfile(".nygc-chromosomes-", tmpdir = dirname(out_dir))
  if (!dir.create(parts)) stop("Could not stage NYGC chromosome files.")
  on.exit(unlink(parts, recursive = TRUE), add = TRUE)
  for (chr in names(intervals)) {
    prepareNYGCGenotypes(input_dir, marker_catalog, as.integer(chr), parts, intervals[[chr]])
  }
  mergeNYGCGenotypes(parts, out_dir, regions)
}

#' Merge NYGC chromosome panels
#'
#' Verify chromosome receipts, PLINK BED layout, marker order, and identical
#' FAM files within each ancestry. Merge in ascending chromosome order and
#' publish one `g1000_<ancestry>` triplet per ancestry and a JSON receipt with
#' checksums, provenance, and counts. Existing directories are never reused.
#'
#' @param chromosome_dir Directory containing outputs of [prepareNYGCGenotypes()].
#' @inheritParams buildNYGCGenotypes
#' @return Invisibly, a data frame with the file prefix and sample and SNP
#'   counts for each ancestry.
#' @export
mergeNYGCGenotypes <- function(chromosome_dir, out_dir, regions = NULL) {
  chromosome_dir <- normalizePath(chromosome_dir, winslash = "/", mustWork = TRUE)
  out_dir <- .gact_nygc_new_directory(out_dir)
  intervals <- .gact_nygc_intervals(regions)
  lock <- paste0(out_dir, ".lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Could not lock NYGC output ", out_dir)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  stage <- file.path(lock, "panel")
  if (!dir.create(stage)) stop("Could not stage merged NYGC files.")

  paths <- lapply(names(intervals), function(chr) {
    .gact_nygc_chromosome_files(chr, intervals[[chr]])
  })
  receipts <- lapply(seq_along(paths), function(i) {
    files <- paths[[i]]
    chr <- names(intervals)[i]
    receipt_path <- file.path(chromosome_dir, files$receipt)
    if (!file.exists(receipt_path)) stop("Missing NYGC chromosome receipt: ", files$receipt)
    receipt <- jsonlite::read_json(receipt_path, simplifyVector = TRUE)
    if (!identical(receipt$build, "GRCh38") ||
        !identical(receipt$release, "20220422_3202_phased_SNV_INDEL_SV") ||
        !identical(as.integer(receipt$chromosome), as.integer(chr)) ||
        !identical(as.numeric(receipt$region), as.numeric(intervals[[chr]])) ||
        length(receipt$inputs) != 4L || anyNA(unlist(receipt$inputs))) {
      stop("Unexpected NYGC provenance in ", files$receipt)
    }
    observed <- unname(tools::md5sum(file.path(chromosome_dir, files$files)))
    expected <- unname(unlist(receipt$outputs)[files$files])
    if (anyNA(observed) || !identical(observed, expected)) {
      stop("NYGC chromosome files are missing or differ from ", files$receipt)
    }
    receipt
  })
  shared_inputs <- lapply(receipts, function(receipt) unname(tail(unlist(receipt$inputs), 2L)))
  if (!all(vapply(shared_inputs, identical, logical(1), shared_inputs[[1L]]))) {
    stop("NYGC chromosome receipts use different sample panels or marker catalogs.")
  }

  populations <- c(AFR = 661L, AMR = 347L, EAS = 504L, EUR = 503L, SAS = 489L)
  panels <- lapply(seq_along(populations), function(i) {
    pop <- names(populations)[i]
    prefix <- paste0("g1000_", tolower(pop))
    parts <- vapply(paths, function(x) x$prefixes[i], character(1))
    counts <- .gact_merge_nygc_panel(chromosome_dir, parts, file.path(stage, prefix),
                                     intervals, populations[[i]])
    expected <- vapply(receipts, function(receipt) {
      row <- match(pop, receipt$panels$ancestry)
      as.integer(receipt$panels$markers[row])
    }, integer(1))
    if (!identical(counts, expected)) stop("NYGC marker counts disagree for ", pop)
    data.frame(ancestry = pop, prefix = prefix, samples = populations[[i]],
               markers = sum(counts))
  })
  panels <- do.call(rbind, panels)
  files <- as.vector(outer(panels$prefix, c(".bed", ".bim", ".fam"), paste0))
  hashes <- unname(tools::md5sum(file.path(stage, files)))
  if (anyNA(hashes)) stop("Could not checksum merged NYGC files.")
  jsonlite::write_json(list(
    build = receipts[[1L]]$build, release = receipts[[1L]]$release,
    chromosomes = as.integer(names(intervals)), regions = regions, panels = panels,
    inputs = lapply(receipts, `[[`, "inputs"),
    outputs = as.list(stats::setNames(hashes, files))
  ), file.path(stage, "nygc-panels.json"), pretty = TRUE, auto_unbox = TRUE)
  if (file.exists(out_dir) || !file.rename(stage, out_dir)) {
    stop("Could not publish merged NYGC files to ", out_dir)
  }
  invisible(panels)
}

.gact_nygc_new_directory <- function(path) {
  parent <- normalizePath(dirname(path), winslash = "/", mustWork = TRUE)
  path <- file.path(parent, basename(path))
  link <- Sys.readlink(path)
  if (!dir.exists(parent) || file.exists(path) || (!is.na(link) && nzchar(link))) {
    stop("NYGC out_dir must be a new directory with an existing parent: ", path)
  }
  path
}

.gact_nygc_intervals <- function(regions) {
  if (is.null(regions)) return(stats::setNames(vector("list", 22L), seq_len(22L)))
  columns <- c("chromosome", "start", "end")
  if (!is.data.frame(regions) || !nrow(regions) || !all(columns %in% names(regions)) ||
      !all(vapply(columns, function(column) {
        x <- regions[[column]]
        is.numeric(x) && all(is.finite(x) & x >= 1 & x <= .Machine$integer.max & x == floor(x))
      }, logical(1))) ||
      any(!regions$chromosome %in% seq_len(22L)) || anyDuplicated(regions$chromosome) ||
      any(regions$start > regions$end)) stop("Specify one valid interval per autosome in regions.")
  rows <- order(regions$chromosome)
  stats::setNames(lapply(rows, function(i) as.integer(c(regions$start[i], regions$end[i]))),
                  as.integer(regions$chromosome[rows]))
}

.gact_merge_nygc_panel <- function(chromosome_dir, parts, prefix, intervals, samples) {
  parts <- file.path(chromosome_dir, parts)
  fam <- paste0(parts, ".fam")
  fam_hash <- unname(tools::md5sum(fam))
  if (anyNA(fam_hash) || any(fam_hash != fam_hash[1L]) ||
      nrow(data.table::fread(fam[1L], header = FALSE)) != samples ||
      !file.copy(fam[1L], paste0(prefix, ".fam"))) {
    stop("Inconsistent NYGC FAM files for ", basename(prefix))
  }
  bed <- file(paste0(prefix, ".bed"), "wb")
  on.exit(close(bed), add = TRUE)
  bim <- file(paste0(prefix, ".bim"), "wb")
  on.exit(close(bim), add = TRUE)
  magic <- as.raw(c(0x6c, 0x1b, 0x01))
  writeBin(magic, bed)
  snps <- lapply(seq_along(parts), function(i) {
    part <- parts[i]
    markers <- data.table::fread(paste0(part, ".bim"), header = FALSE)
    chr <- names(intervals)[i]
    region <- intervals[[i]]
    if (ncol(markers) != 6L || !nrow(markers) || anyNA(markers[[2L]]) ||
        !identical(sub("^chr", "", markers[[1L]]), rep(chr, nrow(markers))) ||
        anyNA(markers[[4L]]) ||
        is.unsorted(markers[[4L]]) || (!is.null(region) &&
        any(markers[[4L]] < region[1L] | markers[[4L]] > region[2L]))) {
      stop("Invalid chromosome marker metadata in ", basename(part), ".bim")
    }
    source_bed <- paste0(part, ".bed")
    size <- 3 + ceiling(samples / 4) * nrow(markers)
    if (file.size(source_bed) != size) stop("Invalid BED size for ", basename(part))
    .gact_append_nygc_file(source_bed, bed, magic)
    .gact_append_nygc_file(paste0(part, ".bim"), bim)
    markers[[2L]]
  })
  if (anyDuplicated(unlist(snps, use.names = FALSE))) {
    stop("Duplicate SNP IDs in merged NYGC panel ", basename(prefix))
  }
  lengths(snps)
}

.gact_append_nygc_file <- function(path, output, header = raw()) {
  input <- file(path, "rb")
  on.exit(close(input))
  if (!identical(readBin(input, "raw", n = length(header)), header)) {
    stop("Invalid variant-major BED header in ", basename(path))
  }
  repeat {
    chunk <- readBin(input, "raw", n = 4L * 1024L * 1024L)
    if (!length(chunk)) break
    writeBin(chunk, output)
  }
}
