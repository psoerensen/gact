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
  panel_name <- sprintf("ALL.chr%d.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz",
                        as.integer(chromosome))
  vcf <- file.path(panel_dir, panel_name)
  fasta <- file.path(fasta_dir, paste0(label, ".fa.gz"))
  rsids <- file.path(rsid_dir, paste0("gnomad_rsid_", label, ".tsv.bgz"))
  panel_manifest <- strsplit(trimws(readLines(file.path(panel_dir,
    "20190312_biallelic_SNV_and_INDEL_MANIFEST.txt"))), "[[:space:]]+")
  fasta_manifest <- strsplit(trimws(readLines(file.path(fasta_dir,
    "md5sum.txt"))), "[[:space:]]+")
  find_checksum <- function(rows, name, file_column, checksum_column, width) {
    matches <- Filter(function(row) length(row) == width &&
      sub("^\\./", "", row[file_column]) == name, rows)
    if (length(matches) != 1L) stop("Expected one checksum for ", name)
    matches[[1L]][checksum_column]
  }
  expected <- c(find_checksum(panel_manifest, panel_name, 1L, 3L, 3L),
                find_checksum(fasta_manifest, basename(fasta), 2L, 1L, 2L))
  inputs <- c(vcf, fasta, rsids)
  checksums <- unname(tools::md5sum(inputs))
  if (anyNA(checksums) || !identical(checksums[1:2], expected)) {
    stop("Missing input or release checksum mismatch for ", label)
  }
  igsr <- paste0("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/",
                 "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/")
  urls <- c(paste0(igsr, panel_name),
            paste0("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/", basename(fasta)),
            "https://gnomad.broadinstitute.org/downloads#v4")
  sources <- unname(Map(function(path, url, md5)
    list(file = basename(path), url = url, md5 = md5), inputs, urls, checksums))

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

.gact_vcf_scan <- function(...) {
  if (!requireNamespace("gsim", quietly = TRUE) ||
      !"gsim_scan_vcf" %in% getNamespaceExports("gsim")) {
    stop("Install gsim with gsim_scan_vcf() support.")
  }
  gsim::gsim_scan_vcf(...)
}
