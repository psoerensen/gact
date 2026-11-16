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
#' @details Requires bcftools on PATH for indexed VCF selection and genetic
#'   map extraction, and gsim for conversion. Incomplete or unphased
#'   genotypes are excluded per ancestry by the gsim importer.
#' @export
prepareNYGCGenotypes <- function(input_dir, marker_catalog, chromosome, out_dir,
                                 region = NULL) {
  if (!requireNamespace("gsim", quietly = TRUE) ||
      !all(c("gsim_import_vcf", "gsim_simulate_pedigree", "gsim_pedigree_from_table") %in%
           getNamespaceExports("gsim"))) stop("Install gsim with phased VCF and BED support.")
  bcftools <- Sys.which("bcftools")
  if (!nzchar(bcftools)) stop("Install bcftools to extract the indexed VCF and map.")
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
  name <- paste0("1kGP_high_coverage_Illumina.", label,
                 ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
  vcf <- file.path(input_dir, name)
  inputs <- c(vcf, paste0(vcf, ".tbi"))
  checksum <- unname(tools::md5sum(inputs))
  manifest <- read.table(file.path(input_dir, "20220804_manifest.txt"),
                         colClasses = "character")
  expected <- manifest[[2L]][match(basename(inputs), manifest[[1L]])]
  if (anyNA(checksum) || !identical(checksum, expected)) {
    stop("NYGC VCF or index checksum differs from the release manifest.")
  }
  sample_file <- file.path(input_dir, "integrated_call_samples_v3.20130502.ALL.panel")
  samples <- data.table::fread(sample_file, header = TRUE, sep = "\t", fill = TRUE,
                              select = c("sample", "super_pop", "gender"))
  populations <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  if (nrow(samples) != 2504L || anyNA(samples) || anyDuplicated(samples$sample) ||
      !setequal(samples$super_pop, populations)) stop("Expected 2,504 unrelated 1000G samples.")

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
  catalog_keys <- with(catalog, paste(pos, nea, ea))
  if (anyDuplicated(catalog_keys)) stop("Catalog allele keys must be unique.")
  targets_file <- file.path(lock, "positions.tsv")
  data.table::fwrite(data.frame(chr = label, pos = sort(unique(catalog$pos))), targets_file,
                     sep = "\t", col.names = FALSE)
  all_samples <- file.path(lock, "samples.txt")
  writeLines(samples$sample, all_samples)
  selected <- file.path(lock, "selected.bcf")
  where <- if (is.null(region)) label else paste0(label, ":", paste(region, collapse = "-"))
  .gact_run_bcftools(bcftools, c("view", "-r", where, "-T", targets_file, "-S", all_samples,
    "-f", "PASS,.", "-m2", "-M2", "-v", "snps", "-Ob", "-o", selected, vcf))
  query <- file.path(lock, "variants.tsv")
  .gact_run_bcftools(bcftools, c("query", "-f", "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/CM\\n",
                                "-o", query, selected))
  variants <- data.table::fread(query, header = FALSE,
    col.names = c("chromosome", "base_pair_position", "id", "ref", "alt", "genetic_position_cm"))
  at <- match(with(variants, paste(base_pair_position, ref, alt)), catalog_keys)
  variants$rsid <- catalog$rsids[at]
  variants <- variants[!is.na(variants$rsid)]
  if (!nrow(variants) || anyDuplicated(variants$id) || anyDuplicated(variants$rsid) ||
      anyNA(variants$genetic_position_cm) || any(!is.finite(variants$genetic_position_cm))) {
    stop("Expected unique matching NYGC alleles and a complete genetic map.")
  }
  ids <- file.path(lock, "variants.txt")
  writeLines(variants$id, ids)
  filtered <- file.path(lock, "matched.vcf.gz")
  .gact_run_bcftools(bcftools, c("view", "-i", paste0("ID=@", ids), "-Oz", "-o", filtered, selected))
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

.gact_run_bcftools <- function(binary, args) {
  status <- system2(binary, shQuote(args))
  if (!identical(status, 0L)) stop("bcftools ", args[1L], " failed with status ", status)
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
