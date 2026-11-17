# Raw reference inputs and chromosome preparation --------

#' Download raw GRCh38 reference inputs
#'
#' Downloads the source files for the IGSR marker backbone or NYGC genotype
#' panels. VCFs, indexes and FASTAs are checked against their release manifests.
#' The NYGC sample panel is checked for the expected sample count, unique
#' IDs and ancestry groups. Existing verified files are reused and conflicting files cause an
#' error. IGSR also prepares rsID lookup tables from gnomAD v4.1 genome VCFs.
#'
#' @param panel One of "IGSR" or "NYGC".
#' @param reference_dir Parent directory for the raw downloads. Defaults to
#'   `reference` in the working directory. For a database created elsewhere,
#'   use `GAlist$dirs[["reference"]]`.
#' @param chromosomes Autosomes to download, from 1 to 22.
#' @param include_rsids For IGSR, prepare gnomAD rsID lookup tables. Set to
#'   FALSE to use separately prepared tables. gnomAD VCFs are processed
#'   one chromosome at a time and discarded after preparation.
#' @return Invisibly, the input paths for the corresponding preparation functions.
#' @seealso [prepareIGSRMarkers()], [buildMarkerDB()], [buildNYGCGenotypes()]
#' @export
downloadReferenceInputs <- function(panel = c("IGSR", "NYGC"), reference_dir = "reference",
                                    chromosomes = 1:22, include_rsids = TRUE) {
  panel <- match.arg(panel)
  if (!length(chromosomes) || !all(chromosomes %in% 1:22) || anyDuplicated(chromosomes))
    stop("Supply unique autosomes from 1 to 22.")
  if (!dir.exists(reference_dir) && !dir.create(reference_dir, recursive = TRUE))
    stop("Cannot create reference directory: ", reference_dir)
  reference_dir <- normalizePath(reference_dir, winslash = "/", mustWork = TRUE)
  lock <- file.path(reference_dir, paste0(".gact-", tolower(panel), "-download"))
  if (!dir.create(lock, mode = "0700", showWarnings = FALSE))
    stop("Cannot lock the ", panel, " download directory: ", lock)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)

  fetch <- function(spec, directory) {
    if (!dir.exists(directory) && !dir.create(directory))
      stop("Cannot create input directory: ", directory)
    manifest <- .gact_fetch_reference_file(spec$url, spec$manifest, directory, lock,
      function(path) .gact_reference_checksums(path, spec))
    checksums <- .gact_reference_checksums(manifest, spec)
    for (i in seq_along(spec$files)) {
      verify <- function(path) {
        if (!identical(unname(tools::md5sum(path)), checksums[i]))
          stop("Release checksum mismatch for ", spec$files[i], " at ", path)
      }
      .gact_fetch_reference_file(spec$url, spec$files[i], directory, lock, verify)
    }
  }
  if (panel == "IGSR") {
    dirs <- list(panel_dir = file.path(reference_dir, "igsr-20190312"),
                 fasta_dir = file.path(reference_dir, "ucsc-hg38"))
    fetch(.gact_reference_spec("IGSR", chromosomes), dirs$panel_dir)
    fetch(.gact_reference_spec("UCSC", chromosomes), dirs$fasta_dir)
    if (isTRUE(include_rsids)) {
      dirs$rsid_dir <- file.path(reference_dir, "gnomad-rsids")
      if (!dir.exists(dirs$rsid_dir) && !dir.create(dirs$rsid_dir))
        stop("Cannot create rsID lookup directory: ", dirs$rsid_dir)
      for (chr in chromosomes) .gact_prepare_gnomad_lookup(chr, dirs$rsid_dir, lock)
    }
  } else {
    dirs <- list(input_dir = file.path(reference_dir, "nygc-30x"))
    fetch(.gact_reference_spec("NYGC", chromosomes), dirs$input_dir)
    .gact_fetch_reference_file(.gact_nygc_sample_url, .gact_nygc_sample_file,
                              dirs$input_dir, lock, .gact_read_nygc_samples)
  }
  invisible(dirs)
}

.gact_nygc_sample_file <- "integrated_call_samples_v3.20130502.ALL.panel"
.gact_nygc_sample_url <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
.gact_nygc_populations <- c("AFR", "AMR", "EAS", "EUR", "SAS")

.gact_read_nygc_samples <- function(path) {
  samples <- data.table::fread(path, header = TRUE, sep = "\t", fill = TRUE,
                               select = c("sample", "super_pop", "gender"))
  if (nrow(samples) != 2504L || anyNA(samples) || any(!nzchar(samples$sample)) ||
      anyDuplicated(samples$sample) || !setequal(samples$super_pop, .gact_nygc_populations))
    stop("Expected 2,504 unrelated 1000G samples across five ancestries.")
  samples
}

.gact_reference_spec <- function(source, chromosomes) {
  igsr <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
  switch(source,
    IGSR = list(
      url = paste0(igsr, "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/"),
      manifest = "20190312_biallelic_SNV_and_INDEL_MANIFEST.txt",
      files = sprintf("ALL.chr%d.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz",
                      as.integer(chromosomes)),
      file_column = 1L, hash_column = 3L),
    UCSC = list(url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/",
      manifest = "md5sum.txt", files = paste0("chr", chromosomes, ".fa.gz"),
      file_column = 2L, hash_column = 1L),
    NYGC = {
      vcf <- paste0("1kGP_high_coverage_Illumina.chr", chromosomes,
                    ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
      list(url = paste0(igsr, "1000G_2504_high_coverage/working/",
                         "20220422_3202_phased_SNV_INDEL_SV/"),
           manifest = "20220804_manifest.txt", files = c(vcf, paste0(vcf, ".tbi")),
           file_column = 1L, hash_column = 2L)
    })
}

.gact_reference_checksums <- function(path, spec) {
  rows <- strsplit(trimws(readLines(path)), "[[:space:]]+")
  rows <- rows[lengths(rows) >= max(spec$file_column, spec$hash_column)]
  names <- vapply(rows, function(row) sub("^\\./", "", row[spec$file_column]), character(1L))
  vapply(spec$files, function(file) {
    match <- which(names == file)
    if (length(match) != 1L) stop("Expected one release checksum for ", file)
    hash <- tolower(rows[[match]][spec$hash_column])
    if (!grepl("^[0-9a-f]{32}$", hash)) stop("Invalid release checksum for ", file)
    hash
  }, character(1L), USE.NAMES = FALSE)
}

.gact_fetch_reference_file <- function(base, name, directory, stage, verify) {
  target <- file.path(directory, name)
  if (!file.exists(target)) {
    message("Downloading ", name)
    incoming <- file.path(stage, name)
    curl::curl_download(paste0(base, name), incoming, quiet = FALSE, mode = "wb")
    verify(incoming)
    if (file.exists(target) || !file.rename(incoming, target))
      stop("Cannot publish reference input: ", target)
  } else {
    verify(target)
  }
  invisible(target)
}

.gact_prepare_gnomad_lookup <- function(chr, directory, stage) {
  label <- paste0("chr", chr)
  target <- file.path(directory, paste0("gnomad_rsid_", label, ".tsv.bgz"))
  receipt <- paste0(target, ".json")
  if (file.exists(target) || file.exists(receipt)) {
    if (!all(file.exists(c(target, receipt)))) stop("Incomplete gnomAD lookup: ", target)
    saved <- jsonlite::read_json(receipt)
    if (!identical(saved$version, "4.1") || !isTRUE(saved$chromosome == chr) ||
        !identical(saved$output_md5, unname(tools::md5sum(target))))
      stop("Existing gnomAD lookup differs from its receipt: ", target)
    return(invisible(target))
  }

  .gact_require_vcf_scan()
  bucket <- "gcp-public-data--gnomad"
  object <- paste0("release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.", label, ".vcf.bgz")
  storage <- "https://storage.googleapis.com/"
  url <- paste0(storage, bucket, "/", object)
  metadata_url <- paste0(storage, "storage/v1/b/", bucket, "/o/",
                         utils::URLencode(object, reserved = TRUE), "?fields=md5Hash")
  response <- curl::curl_fetch_memory(metadata_url)
  if (response$status_code != 200L) stop("Cannot read gnomAD checksum for ", label)
  hash <- jsonlite::fromJSON(rawToChar(response$content))$md5Hash
  if (length(hash) != 1L || is.na(hash) || length(openssl::base64_decode(hash)) != 16L)
    stop("Missing gnomAD VCF MD5 for ", label)
  expected <- paste(sprintf("%02x", as.integer(openssl::base64_decode(hash))), collapse = "")
  vcf <- file.path(stage, basename(object))
  message("Downloading gnomAD v4.1 genome sites for ", label)
  curl::curl_download(url, vcf, quiet = FALSE, mode = "wb")
  if (!identical(unname(tools::md5sum(vcf)), expected))
    stop("gnomAD VCF checksum mismatch for ", label)

  message("Preparing gnomAD rsID lookup for ", label)
  lookup <- file.path(stage, basename(target))
  con <- gzfile(lookup, "wt", compression = 1)
  on.exit(if (!is.null(con)) close(con), add = TRUE)
  writeLines("#CHROM\tPOS\tID\tREF\tALT", con)
  last_position <- 0
  bases <- c("A", "C", "G", "T")
  .gact_vcf_scan(vcf, function(batch) {
    if (any(batch$chrom != label) || is.unsorted(c(last_position, batch$pos)))
      stop("Unexpected chromosome or position order in gnomAD ", label)
    last_position <<- tail(batch$pos, 1L)
    keep <- batch$ref %in% bases & batch$alt %in% bases & batch$ref != batch$alt &
      grepl("(^|;)rs[0-9]+(;|$)", batch$id)
    if (any(keep)) {
      writeLines(paste(batch$chrom[keep], as.integer(batch$pos[keep]), batch$id[keep],
                       batch$ref[keep], batch$alt[keep], sep = "\t"), con)
    }
  })
  close(con)
  con <- NULL
  staged_receipt <- file.path(stage, basename(receipt))
  jsonlite::write_json(list(version = "4.1", chromosome = as.integer(chr), url = url,
    source_md5 = expected, output_md5 = unname(tools::md5sum(lookup))),
    staged_receipt, auto_unbox = TRUE, pretty = TRUE)
  suspendInterrupts({
    if (!file.rename(lookup, target)) stop("Cannot publish gnomAD lookup: ", target)
    if (!file.rename(staged_receipt, receipt)) {
      unlink(target)
      stop("Cannot publish gnomAD lookup receipt: ", receipt)
    }
  })
  unlink(vcf)
  invisible(target)
}

#' Prepare one IGSR GRCh38 chromosome for the marker catalog
#'
#' Stream the IGSR VCF using gsim, retain PASS biallelic SNPs, and check
#' the reference allele and the global minor allele frequency. A variant
#' must resolve to exactly one distinct rsID across the source VCF and the
#' exact chromosome, position, REF and ALT matches in the rsID lookup.
#'
#' @param chromosome An autosome number from 1 to 22.
#' @param panel_dir Directory containing the IGSR 2019-03-12 VCFs and
#'   their release checksum manifest.
#' @param fasta_dir Directory containing the UCSC GRCh38 chromosome FASTAs
#'   and `md5sum.txt`.
#' @param rsid_dir Directory containing position-sorted
#'   `gnomad_rsid_chrN.tsv.bgz` lookup tables with columns
#'   `#CHROM`, `POS`, `ID`, `REF`, and `ALT`.
#' @param out_dir Directory for the chromosome table and provenance receipt.
#'   The output files must be new.
#' @param min_maf Minimum global minor allele frequency.
#' @return Invisibly, the chromosome preparation receipt.
#' @seealso [buildMarkerDB()]
#' @export
prepareIGSRMarkers <- function(chromosome, panel_dir, fasta_dir, rsid_dir,
                               out_dir, min_maf = 0.01) {
  if (length(chromosome) != 1L || !chromosome %in% seq_len(22L)) {
    stop("Specify an autosome from 1 to 22.")
  }
  if (!is.numeric(min_maf) || length(min_maf) != 1L ||
      !is.finite(min_maf) || min_maf < 0 || min_maf > 0.5) {
    stop("min_maf must be between 0 and 0.5.")
  }
  label <- paste0("chr", chromosome)
  panel <- .gact_reference_spec("IGSR", chromosome)
  reference_spec <- .gact_reference_spec("UCSC", chromosome)
  vcf <- file.path(panel_dir, panel$files)
  fasta <- file.path(fasta_dir, reference_spec$files)
  rsids <- file.path(rsid_dir, paste0("gnomad_rsid_", label, ".tsv.bgz"))
  panel_manifest_path <- file.path(panel_dir, panel$manifest)
  fasta_manifest_path <- file.path(fasta_dir, reference_spec$manifest)
  inputs <- c(vcf, fasta, rsids)
  required <- c(inputs, panel_manifest_path, fasta_manifest_path)
  missing <- required[!file.exists(required)]
  if (length(missing)) stop("Missing IGSR preparation inputs for ", label, ":\n",
                            paste(missing, collapse = "\n"))
  expected <- c(.gact_reference_checksums(panel_manifest_path, panel),
                .gact_reference_checksums(fasta_manifest_path, reference_spec))
  checksums <- unname(tools::md5sum(inputs))
  if (anyNA(checksums) || !identical(checksums[1:2], expected)) {
    stop("Missing input or release checksum mismatch for ", label)
  }
  urls <- c(paste0(panel$url, panel$files),
            paste0(reference_spec$url, reference_spec$files),
            "https://gnomad.broadinstitute.org/downloads#v4")
  sources <- unname(Map(function(path, url, md5)
    list(file = basename(path), url = url, md5 = md5), inputs, urls, checksums))
  lookup_receipt <- paste0(rsids, ".json")
  if (file.exists(lookup_receipt)) {
    source <- jsonlite::read_json(lookup_receipt)
    if (!identical(source$output_md5, checksums[3L]))
      stop("gnomAD lookup checksum differs from its receipt for ", label)
    sources[[3L]]$upstream <- source
  }

  if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE)) {
    stop("Cannot create chromosome output directory.")
  }
  targets <- file.path(out_dir, paste0(label, c(".tsv.gz", ".json")))
  stage <- file.path(out_dir, paste0(".", label, ".lock"))
  if (!dir.create(stage, showWarnings = FALSE)) stop("Chromosome preparation is already staged.")
  on.exit(unlink(stage, recursive = TRUE, expand = FALSE), add = TRUE)
  if (any(file.exists(targets))) stop("Chromosome output already exists.")

  read_fasta <- function() {
    input <- gzfile(fasta, "rt")
    on.exit(close(input))
    if (!identical(trimws(readLines(input, n = 1L)), paste0(">", label))) {
      stop("Unexpected FASTA chromosome.")
    }
    charToRaw(toupper(paste0(readLines(input), collapse = "")))
  }
  reference <- read_fasta()
  bases <- c("A", "C", "G", "T")
  base_bytes <- as.integer(charToRaw("ACGT"))
  lookup <- gzfile(rsids, "rt")
  on.exit(close(lookup), add = TRUE)
  if (!identical(trimws(readLines(lookup, n = 1L)), "#CHROM\tPOS\tID\tREF\tALT")) {
    stop("Unexpected rsID lookup header.")
  }
  annotations <- data.table::data.table(pos = integer(), id = character(),
                                        ref = character(), alt = character())
  lookup_last <- 0L
  lookup_eof <- FALSE
  vcf_last <- 0
  counts <- c(input = 0L, not_pass_biallelic_snv = 0L, below_min_maf = 0L,
              reference_mismatch = 0L, unresolved_rsid = 0L, retained = 0L)
  plain <- file.path(stage, paste0(label, ".tsv"))
  writeLines("rsids\tchr\tpos\tref\talt\taf", plain)

  consume <- function(batch) {
    if (any(batch$chrom != as.character(chromosome) | batch$pos > length(reference)) ||
        is.unsorted(c(vcf_last, batch$pos))) {
      stop("Unexpected chromosome or unsorted/out-of-range VCF position.")
    }
    vcf_last <<- tail(batch$pos, 1L)
    counts["input"] <<- counts["input"] + nrow(batch)
    pass <- batch$filter == "PASS" & batch$ref %in% bases &
      batch$alt %in% bases & batch$ref != batch$alt
    counts["not_pass_biallelic_snv"] <<- counts["not_pass_biallelic_snv"] + sum(!pass)
    batch <- batch[pass, , drop = FALSE]
    ac <- suppressWarnings(as.numeric(batch$info_AC))
    an <- suppressWarnings(as.numeric(batch$info_AN))
    if (any(!grepl("^[0-9]+$", batch$info_AC) | !grepl("^[0-9]+$", batch$info_AN) |
            !is.finite(ac) | !is.finite(an) | an <= 0 | ac > an)) {
      stop("Invalid IGSR AC or AN.")
    }
    af <- ac / an
    frequent <- pmin(af, 1 - af) >= min_maf
    counts["below_min_maf"] <<- counts["below_min_maf"] + sum(!frequent)
    batch <- batch[frequent, , drop = FALSE]
    af <- af[frequent]
    same_ref <- as.integer(reference[batch$pos]) == base_bytes[match(batch$ref, bases)]
    counts["reference_mismatch"] <<- counts["reference_mismatch"] + sum(!same_ref)
    batch <- batch[same_ref, , drop = FALSE]
    af <- af[same_ref]
    if (!nrow(batch)) return(invisible(NULL))

    first <- batch$pos[1L]
    last <- tail(batch$pos, 1L)
    if (nrow(annotations)) annotations <<- annotations[annotations$pos >= first]
    while (!lookup_eof && lookup_last <= last) {
      lines <- readLines(lookup, n = 65536L)
      if (!length(lines)) {
        lookup_eof <<- TRUE
        break
      }
      chunk <- data.table::fread(text = paste0(paste(lines, collapse = "\n"), "\n"),
        sep = "\t", header = FALSE, quote = "", na.strings = NULL,
        colClasses = c("character", "integer",
                                                               rep("character", 3L)),
        col.names = c("chrom", "pos", "id", "ref", "alt"))
      if (nrow(chunk) != length(lines) || anyNA(chunk) || any(chunk$chrom != label) ||
          any(chunk$pos < 1L) || is.unsorted(c(lookup_last, chunk$pos))) {
        stop("Invalid or unsorted rsID lookup records.")
      }
      lookup_last <<- tail(chunk$pos, 1L)
      annotations <<- data.table::rbindlist(list(annotations,
        chunk[chunk$pos >= first, c("pos", "id", "ref", "alt")]))
    }

    keys <- data.table::data.table(pos = batch$pos, ref = batch$ref, alt = batch$alt)
    hits <- annotations[keys, on = c("pos", "ref", "alt"), which = TRUE,
                        nomatch = NULL, allow.cartesian = TRUE]
    # x[i] may return several annotations for a VCF row. Group by allele key
    # before matching back so ties and repeated positions remain unambiguous.
    lookup_matches <- unique(annotations[hits])
    tokens <- strsplit(lookup_matches$id, ";", fixed = TRUE)
    called <- data.table::data.table(
      pos = rep(lookup_matches$pos, lengths(tokens)),
      ref = rep(lookup_matches$ref, lengths(tokens)),
      alt = rep(lookup_matches$alt, lengths(tokens)),
      id = as.character(unlist(tokens)))
    called <- called[grepl("^rs[0-9]+$", called$id)]
    distinct <- called[, list(ids = list(unique(id))), by = c("pos", "ref", "alt")]
    matched <- distinct[keys, on = c("pos", "ref", "alt")]$ids
    ids <- vapply(seq_len(nrow(batch)), function(i) {
      from_vcf <- batch$id[i]
      candidates <- unique(c(matched[[i]], from_vcf[grepl("^rs[0-9]+$", from_vcf)]))
      candidates <- candidates[!is.na(candidates)]
      if (length(candidates) == 1L) candidates else NA_character_
    }, character(1L))
    resolved <- !is.na(ids)
    counts["unresolved_rsid"] <<- counts["unresolved_rsid"] + sum(!resolved)
    counts["retained"] <<- counts["retained"] + sum(resolved)
    if (any(resolved)) {
      data.table::fwrite(data.frame(rsids = ids[resolved], chr = chromosome,
        pos = batch$pos[resolved], ref = batch$ref[resolved], alt = batch$alt[resolved],
        af = sprintf("%.17g", af[resolved])), plain, sep = "\t", append = TRUE,
        col.names = FALSE)
    }
    invisible(NULL)
  }
  header <- .gact_vcf_scan(vcf, consume, info = c("AC", "AN"))
  reference_header <- header[startsWith(header, "##reference=")]
  if (length(reference_header) != 1L || !grepl("GRCh38", reference_header, fixed = TRUE)) {
    stop("Expected a GRCh38 IGSR VCF.")
  }
  if (counts["retained"] == 0L) stop("No eligible variants on ", label)
  staged <- file.path(stage, basename(targets))
  R.utils::gzip(plain, destname = staged[1L], remove = TRUE)
  receipt <- list(chromosome = as.integer(chromosome), build = "GRCh38",
    min_maf = min_maf, counts = as.list(counts), sources = sources,
    output = list(file = basename(targets[1L]), md5 = unname(tools::md5sum(staged[1L]))))
  jsonlite::write_json(receipt, staged[2L], pretty = TRUE, auto_unbox = TRUE)
  published <- FALSE
  committed <- FALSE
  on.exit(if (published && !committed) unlink(targets[1L]), add = TRUE)
  suspendInterrupts({
    .publish_gwas_import_file(staged[1L], targets[1L])
    published <- TRUE
    .publish_gwas_import_file(staged[2L], targets[2L])
    committed <- TRUE
  })
  invisible(receipt)
}

.gact_require_vcf_scan <- function() {
  if (!requireNamespace("gsim", quietly = TRUE) ||
      !"gsim_scan_vcf" %in% getNamespaceExports("gsim")) {
    stop("Install gsim with gsim_scan_vcf() support.")
  }
}

.gact_vcf_scan <- function(...) {
  .gact_require_vcf_scan()
  gsim::gsim_scan_vcf(...)
}

# Marker catalog build and release packaging --------

.gact_grch38_marker_tables <- function(dbdir) {
 roots <- unique(c(dbdir, dirname(dbdir)))
 reference <- file.path(roots, "reference")
 runs <- as.character(unlist(lapply(reference, list.dirs, full.names = TRUE, recursive = FALSE)))
 runs <- runs[startsWith(basename(runs), "grch38-prep.")]
 candidates <- c(file.path(reference, "prepared_grch38"),
                 file.path(roots, "prepared_grch38"), file.path(runs, "tables"))
 files <- paste0("chr", rep(1:22, 2L), rep(c(".tsv.gz", ".json"), each = 22L))
 ready <- vapply(candidates, function(path) all(file.exists(file.path(path, files))),
                logical(1L))
 if (!any(ready)) return(NULL)
 tables <- unique(normalizePath(candidates[ready], winslash = "/", mustWork = TRUE))
 if (length(tables) > 1L)
  stop("Found multiple complete GRCh38 preparations: ", paste(tables, collapse = ", "))
 tables
}

.gact_prepare_grch38_marker_tables <- function(reference_dir) {
  out_dir <- file.path(reference_dir, "grch38-prep.auto", "tables")
  if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE, showWarnings = FALSE) &&
      !dir.exists(out_dir)) {
    stop("Cannot create GRCh38 preparation directory: ", out_dir)
  }
  lock <- file.path(out_dir, ".prepare.lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Cannot lock GRCh38 preparation: ", lock)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)

  chromosomes <- seq_len(22L)
  labels <- paste0("chr", chromosomes)
  has_table <- file.exists(file.path(out_dir, paste0(labels, ".tsv.gz")))
  has_receipt <- file.exists(file.path(out_dir, paste0(labels, ".json")))
  incomplete <- labels[xor(has_table, has_receipt)]
  if (length(incomplete)) {
    stop("Incomplete GRCh38 chromosome pairs in ", out_dir, ": ",
         paste(incomplete, collapse = ", "),
         ". Move those partial outputs aside and retry.")
  }
  missing <- chromosomes[!has_table]
  if (!length(missing)) return(out_dir)
  .gact_require_vcf_scan()
  message("Preparing GRCh38 chromosome tables in ", out_dir, " (",
          length(missing), " of ", length(chromosomes), " to prepare).")
  for (chr in missing) {
    message("Chromosome ", chr, "/22: downloading or verifying reference inputs")
    inputs <- downloadReferenceInputs("IGSR", reference_dir, chromosomes = chr)
    message("Chromosome ", chr, "/22: preparing marker table")
    prepareIGSRMarkers(chr, panel_dir = inputs$panel_dir, fasta_dir = inputs$fasta_dir,
                       rsid_dir = inputs$rsid_dir, out_dir = out_dir)
    message("Chromosome ", chr, "/22: marker table ready")
  }
  out_dir
}

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

#' Package a GRCh38 marker database for release
#'
#' Create a Zenodo-ready gzip tar archive containing the marker catalog,
#' coordinate-and-allele index, and provenance. Verify the catalog and index
#' against the build metadata, and include SHA-256 checksums for the archived
#' files and the archive itself.
#'
#' @param dbdir Existing GRCh38 database directory containing the marker directory.
#' @param out_dir Output directory outside the source database.
#' @return Invisibly, a list containing the archive path, checksum file path,
#'   and archive SHA-256.
#' @seealso [buildMarkerDB()]
#' @export
packageMarkerDB <- function(dbdir, out_dir = getwd()) {
  dbdir <- normalizePath(dbdir, winslash = "/", mustWork = TRUE)
  out_dir <- .gact_gwas_output_path(path.expand(out_dir))
  source_key <- .gact_gwas_path_key(dbdir)
  output_key <- .gact_gwas_path_key(out_dir)
  if (output_key == source_key || startsWith(output_key, paste0(source_key, "/"))) {
    stop("out_dir must be outside the source database")
  }
  filenames <- c("markers.txt.gz", "cpra.txt.gz", "marker_manifest.json")
  source <- file.path(dbdir, "marker", filenames)
  if (!all(file.exists(source) & !dir.exists(source))) {
    stop("The source database must contain markers.txt.gz, cpra.txt.gz, and marker_manifest.json")
  }
  if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create output directory: ", out_dir)
  }
  name <- "gact-hsa.GRCh38.0.0.1-marker.tar.gz"
  archive <- file.path(out_dir, name)
  checksum_file <- paste0(archive, ".sha256")
  stage <- file.path(out_dir, paste0(".", name, ".lock"))
  if (!dir.create(stage, mode = "0700", showWarnings = FALSE)) stop("Cannot stage marker archive: ", stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  if (any(file.exists(c(archive, checksum_file)))) stop("Marker archive or checksum already exists in ", out_dir)
  payload <- file.path(stage, "payload")
  if (!dir.create(file.path(payload, "marker"), recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create archive payload directory")
  }
  members <- paste0("marker/", filenames)
  copies <- file.path(payload, members)
  if (!all(file.copy(source, copies))) stop("Cannot copy marker release files")
  metadata <- jsonlite::read_json(copies[3L])
  expected <- unlist(metadata$checksums)
  if (!identical(metadata$build, "GRCh38") ||
      !setequal(names(expected), filenames[1:2]) || length(expected) != 2L ||
      !identical(unname(tools::md5sum(file.path(payload, "marker", names(expected)))), unname(expected))) {
    stop("GRCh38 marker metadata or build checksums do not match the source files")
  }
  hashes <- vapply(copies, .gact_sha256_file, character(1), USE.NAMES = FALSE)
  writeLines(paste(hashes, members, sep = "  "), file.path(payload, "SHA256SUMS"))
  staged_archive <- file.path(stage, name)
  old <- setwd(payload)
  status <- tryCatch(
    utils::tar(staged_archive, c(members, "SHA256SUMS"), compression = "gzip", tar = "internal"),
    finally = setwd(old)
  )
  if (status != 0L) stop("Cannot create marker archive")
  unpack <- file.path(stage, "verify")
  if (!dir.create(unpack) || utils::untar(staged_archive, exdir = unpack, tar = "internal") != 0L ||
      !identical(vapply(file.path(unpack, members), .gact_sha256_file, character(1), USE.NAMES = FALSE), hashes) ||
      !identical(readLines(file.path(unpack, "SHA256SUMS")), paste(hashes, members, sep = "  "))) {
    stop("Marker archive did not pass readback verification")
  }
  checksum <- .gact_sha256_file(staged_archive)
  staged_checksum <- paste0(staged_archive, ".sha256")
  writeLines(paste(checksum, name, sep = "  "), staged_checksum)
  suspendInterrupts({
    if (!file.rename(staged_archive, archive)) stop("Cannot publish marker archive")
    if (!file.rename(staged_checksum, checksum_file)) {
      unlink(archive)
      stop("Cannot publish marker archive checksum")
    }
  })
  invisible(list(archive = archive, checksum_file = checksum_file, sha256 = checksum))
}

# NYGC genotype preparation and merging --------

#' Prepare ancestry-specific NYGC GRCh38 genotype files
#'
#' Converts one chromosome from the NYGC 30x phased release to PLINK
#' BED/BIM/FAM using gsim. Retains biallelic SNPs that match the catalog by
#' chromosome, position, REF and ALT. ALT is the BED dosage allele.
#' Genotypes come from the 2,504 unrelated Phase 3 samples, split by
#' super-population. The catalog and GAlist are read only.
#'
#' @param input_dir Directory with the official NYGC VCF, tabix index,
#'   `20220804_manifest.txt` and `integrated_call_samples_v3.20130502.ALL.panel`.
#' @param marker_catalog Path to the GRCh38 `markers.txt.gz` file with ALT
#'   as the effect allele.
#' @param chromosome An autosome number.
#' @param out_dir Existing staging directory for the chromosome files.
#'   Existing files are never overwritten.
#' @param region Optional inclusive `c(start, end)` for a bounded test.
#' @return Invisibly, a per-ancestry summary of sample and marker counts.
#' @details Requires gsim with site scanning and genotype conversion.
#'   Region requests scan the chromosome VCF sequentially. Incomplete or
#'   unphased genotypes are excluded per ancestry by the gsim importer.
#' @export
prepareNYGCGenotypes <- function(input_dir, marker_catalog, chromosome, out_dir,
                                 region = NULL) {
  if (!requireNamespace("gsim", quietly = TRUE) ||
      !all(c("gsim_import_vcf", "gsim_simulate_pedigree", "gsim_pedigree_from_table") %in%
           getNamespaceExports("gsim"))) stop("Install gsim with phased VCF and BED support.")
  if (length(chromosome) != 1L || !chromosome %in% seq_len(22L)) {
    stop("Specify an autosome from 1 to 22.")
  }
  if (!is.null(region) && (!is.numeric(region) || length(region) != 2L ||
      any(!is.finite(region) | region < 1L | region > .Machine$integer.max |
          region != floor(region)) ||
      region[1L] > region[2L])) stop("Specify region as c(start, end).")
  if (!is.null(region)) region <- as.integer(region)
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  if (!dir.exists(out_dir)) stop("out_dir must be an existing directory.")
  label <- paste0("chr", chromosome)
  spec <- .gact_reference_spec("NYGC", chromosome)
  inputs <- file.path(input_dir, spec$files)
  vcf <- inputs[1L]
  manifest_file <- file.path(input_dir, spec$manifest)
  sample_file <- file.path(input_dir, .gact_nygc_sample_file)
  required <- c(inputs, manifest_file, sample_file)
  missing <- required[!file.exists(required)]
  if (length(missing)) stop("Missing NYGC preparation inputs for ", label, ":\n",
                            paste(missing, collapse = "\n"))
  checksum <- unname(tools::md5sum(inputs))
  expected <- .gact_reference_checksums(manifest_file, spec)
  if (anyNA(checksum) || !identical(checksum, expected)) {
    stop("NYGC VCF or index checksum differs from the release manifest.")
  }
  samples <- .gact_read_nygc_samples(sample_file)
  populations <- .gact_nygc_populations

  paths <- .gact_nygc_chromosome_files(chromosome, region)
  prefixes <- paths$prefixes
  files <- c(paths$files, paths$receipt)
  lock <- file.path(out_dir, paste0(".", tail(files, 1L), ".lock"))
  if (!dir.create(lock, showWarnings = FALSE)) stop("Another genotype preparation holds ", lock)
  on.exit(unlink(lock, recursive = TRUE, expand = FALSE), add = TRUE)
  targets <- file.path(out_dir, files)
  if (any(file.exists(targets))) stop("Genotype output files already exist in ", out_dir)

  catalog <- data.table::fread(marker_catalog, select = c("rsids", "chr", "pos", "ea", "nea"))
  catalog <- catalog[catalog$chr == chromosome]
  if (!is.null(region)) catalog <- catalog[catalog$pos >= region[1L] & catalog$pos <= region[2L]]
  if (!nrow(catalog) || anyDuplicated(catalog$rsids) || anyNA(catalog)) {
    stop("Expected unique catalog markers in the requested region.")
  }
  if (anyDuplicated(catalog, by = c("pos", "nea", "ea"))) {
    stop("Catalog allele keys must be unique.")
  }
  filtered <- file.path(lock, "matched.vcf.gz")
  variants <- .gact_select_nygc_variants(vcf, catalog, label, region, filtered)
  if (!nrow(variants) || anyDuplicated(variants$id) || anyDuplicated(variants$rsid) ||
      any(!is.finite(variants$genetic_position_cm))) {
    stop("Expected unique matching NYGC alleles and a complete genetic map.")
  }
  map_columns <- c("chromosome", "base_pair_position", "genetic_position_cm")
  map <- as.data.frame(variants[, map_columns, with = FALSE])
  same_position <- duplicated(map$base_pair_position)
  if (any(map$genetic_position_cm != map$genetic_position_cm[
      match(map$base_pair_position, map$base_pair_position)])) {
    stop("NYGC provides conflicting map distances at the same position.")
  }
  map <- map[!same_position, , drop = FALSE]
  stats <- lapply(seq_along(populations), function(i) {
    pop <- populations[i]
    panel <- samples[samples$super_pop == pop]
    hap <- gsim::gsim_import_vcf(filtered, map, file.path(lock, paste0("hap-", pop)),
                                samples = panel$sample, unsupported = "skip")
    ped <- gsim::gsim_pedigree_from_table(data.frame(animal = panel$sample, sire = NA_character_,
      dam = NA_character_, sex = ifelse(panel$gender == "male", "M", "F")))
    prefix <- file.path(lock, prefixes[i])
    gsim::gsim_simulate_pedigree(hap, ped, seed = 1L, output = prefix, format = "bed")
    bim <- data.table::fread(paste0(prefix, ".bim"), header = FALSE)
    fam <- data.table::fread(paste0(prefix, ".fam"), header = FALSE)
    matches <- match(bim[[2L]], variants$id)
    if (anyNA(matches) || !identical(as.character(bim[[1L]]), variants$chromosome[matches]) ||
        !identical(as.integer(bim[[4L]]), as.integer(variants$base_pair_position[matches])) ||
        !identical(bim[[5L]], variants$alt[matches]) ||
        !identical(bim[[6L]], variants$ref[matches]) ||
        !setequal(fam[[2L]], panel$sample) || nrow(fam) != nrow(panel)) {
      stop("BED metadata does not match NYGC samples and catalog alleles for ", pop)
    }
    bim[[2L]] <- variants$rsid[matches]
    data.table::fwrite(bim, paste0(prefix, ".bim"), sep = "\t", col.names = FALSE)
    data.frame(ancestry = pop, samples = nrow(fam), markers = nrow(bim))
  })
  result <- do.call(rbind, stats)
  jsonlite::write_json(list(build = "GRCh38", release = "20220422_3202_phased_SNV_INDEL_SV",
    chromosome = chromosome, region = region, matched = nrow(variants), panels = result,
    inputs = as.list(tools::md5sum(c(inputs, sample_file, marker_catalog))),
    outputs = as.list(stats::setNames(unname(tools::md5sum(file.path(lock, head(files, -1L)))),
                                      head(files, -1L)))),
    file.path(lock, tail(files, 1L)), pretty = TRUE, auto_unbox = TRUE)
  published <- character()
  committed <- FALSE
  on.exit(if (!committed) unlink(published, expand = FALSE), add = TRUE)
  for (i in seq_along(files)) {
    .publish_gwas_import_file(file.path(lock, files[i]), targets[i])
    published <- c(published, targets[i])
  }
  committed <- TRUE
  invisible(result)
}

.gact_select_nygc_variants <- function(vcf, catalog, chromosome, region, output) {
  data.table::setkeyv(catalog, c("pos", "nea", "ea"))
  selected <- list()
  bases <- c("A", "C", "G", "T")
  .gact_vcf_scan(vcf, function(batch) {
    keys <- data.table::as.data.table(batch[c("pos", "ref", "alt")])
    at <- catalog[keys, on = c(pos = "pos", nea = "ref", ea = "alt"), which = TRUE]
    keep <- !is.na(at) & batch$filter %in% c("PASS", ".") &
      batch$ref %in% bases & batch$alt %in% bases & batch$ref != batch$alt
    if (any(keep)) {
      selected[[length(selected) + 1L]] <<- data.table::data.table(
        chromosome = batch$chrom[keep], base_pair_position = batch$pos[keep],
        id = batch$id[keep], ref = batch$ref[keep], alt = batch$alt[keep],
        genetic_position_cm = suppressWarnings(as.numeric(batch$info_CM[keep])),
        rsid = catalog$rsids[at[keep]])
    }
    keep
  }, info = "CM", chromosome = chromosome, region = region, batch_size = 1024L,
  output = output)
  data.table::rbindlist(selected)
}

.gact_nygc_chromosome_files <- function(chromosome, region) {
  label <- paste0("chr", chromosome)
  suffix <- if (is.null(region)) "" else paste(region, collapse = "-")
  prefixes <- paste0("g1000_", c("afr", "amr", "eas", "eur", "sas"), ".", label,
                     if (nzchar(suffix)) paste0(".", suffix))
  list(prefixes = prefixes,
       files = as.vector(outer(prefixes, c(".bed", ".bim", ".fam"), paste0)),
       receipt = paste0("nygc-", label, if (nzchar(suffix)) paste0("-", suffix), ".json"))
}

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

# Published NYGC panel installation --------

.gact_nygc_release_hashes <- function(record) {
  manifest <- system.file("extdata", paste0("nygc-", record, "-sha256.txt"),
                          package = "gact", mustWork = TRUE)
  checksums <- utils::read.table(manifest, colClasses = "character")
  stats::setNames(checksums[[1L]], checksums[[2L]])
}

.gact_sha256_file <- function(path) {
  connection <- file(path, "rb")
  on.exit(close(connection))
  as.character(openssl::sha256(connection))
}

.gact_download_nygc_archive <- function(record, archive, dest) {
  url <- paste0("https://zenodo.org/api/records/", record, "/files/", archive, "/content")
  curl::curl_download(url, dest, mode = "wb",
                      handle = curl::new_handle(connecttimeout = 30, timeout = 3000))
}

.gact_install_nygc_panels <- function(dest, doi) {
  if (length(dest) != 1L || !dir.exists(dest)) {
    stop("The GRCh38 GAlist must have an existing marker directory.")
  }
  record <- sub("^10\\.5281/zenodo\\.", "", doi)
  sha256 <- .gact_nygc_release_hashes(record)
  ancestries <- c("afr", "amr", "eas", "eur", "sas")
  triplets <- outer(c(".bed", ".bim", ".fam"), paste0("g1000_", ancestries),
                    function(extension, prefix) paste0(prefix, extension))
  filenames <- as.vector(triplets)
  targets <- file.path(dest, filenames)

  stage <- file.path(dest, ".nygc-install.lock")
  if (!dir.create(stage, showWarnings = FALSE)) stop("Cannot lock NYGC installation: ", stage)
  on.exit(unlink(stage, recursive = TRUE, expand = FALSE), add = TRUE)

  existing <- file.exists(targets)
  links <- Sys.readlink(targets)
  if (any(!is.na(links) & nzchar(links))) stop("NYGC destination contains a symbolic link.")
  for (i in which(existing)) {
    if (dir.exists(targets[i]) ||
        .gact_sha256_file(targets[i]) != sha256[[filenames[i]]]) {
      stop("Existing NYGC destination differs from the published release: ", targets[i])
    }
  }
  if (all(existing)) return(invisible(targets))

  for (i in which(colSums(matrix(existing, nrow = 3L)) != 3L)) {
    archive <- paste0("gact-nygc30x-grch38-", ancestries[i], ".tar.gz")
    archive_path <- file.path(stage, archive)
    message("Downloading GRCh38 NYGC panel: ", toupper(ancestries[i]))
    .gact_download_nygc_archive(record, archive, archive_path)
    if (.gact_sha256_file(archive_path) != sha256[[archive]]) {
      stop("SHA-256 mismatch for ", archive)
    }
    members <- triplets[, i]
    status <- utils::untar(archive_path, files = members, exdir = stage, tar = "internal")
    if (status != 0L) stop("Cannot unpack ", archive)
    for (member in members) {
      path <- file.path(stage, member)
      if (!file.exists(path) || .gact_sha256_file(path) != sha256[[member]]) {
        stop("SHA-256 mismatch for extracted file ", member)
      }
    }
    unlink(archive_path)
  }

  published <- character()
  committed <- FALSE
  on.exit(if (!committed) unlink(published, expand = FALSE), add = TRUE, after = FALSE)
  for (i in which(!existing)) {
    if (file.exists(targets[i]) || !file.rename(file.path(stage, filenames[i]), targets[i])) {
      stop("Cannot install NYGC file: ", targets[i])
    }
    published <- c(published, targets[i])
  }
  committed <- TRUE
  invisible(targets)
}
