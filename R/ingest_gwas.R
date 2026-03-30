#' Detect GWAS summary-statistics schema using column names and values
#'
#' Scores candidate columns for canonical GWAS fields and returns a mapping
#' with confidence scores, alternatives, and warnings.
#'
#' Canonical fields are: `marker`, `chr`, `pos`, `ea`, `nea`, `eaf`, `b`,
#' `seb`, `p`, `n`, `ncase`, `ncontrol`, and `info`.
#'
#' @param stat A `data.frame` containing summary statistics.
#' @param sample_n Optional integer; number of rows to sample for detection.
#' @return A list with elements `mapping`, `confidence`, `candidates`,
#'   `warnings`, and `details`.
#' @export

detectStatSchema <- function(stat, sample_n = 100000) {
  if (is.null(stat)) stop("`stat` is NULL")
  if (!is.data.frame(stat)) stat <- as.data.frame(stat)
  if (!nrow(stat)) stop("`stat` has zero rows")

  if (!is.null(sample_n) && is.finite(sample_n) && nrow(stat) > sample_n) {
    idx <- sort(sample.int(nrow(stat), sample_n))
    stat <- stat[idx, , drop = FALSE]
  }

  col_names <- colnames(stat)
  clean_names <- tolower(gsub("[^a-z0-9]", "", col_names))

  name_patterns <- list(
    marker = c("^rsid$", "^rsids$", "^snp$", "^marker$", "^variantid$", "^id$"),
    chr = c("^chr$", "^chrom$", "^chromosome$"),
    pos = c("^pos$", "^bp$", "^position$", "^basepairposition$", "^bppos$"),
    ea = c("^a1$", "^ea$", "^effectallele$", "^alt$", "^allele1$"),
    nea = c("^a2$", "^nea$", "^otherallele$", "^ref$", "^non\.?effectallele$", "^allele2$"),
    eaf = c("^eaf$", "^effectallelefreq", "^effectallelefrequency$", "^af$", "^maf$", "^freq$"),
    b = c("^b$", "^beta$", "^effect$", "^estimate$", "^logor$", "^or$"),
    seb = c("^se$", "^seb$", "^stderr$", "^standarderror$", "^sebeta$"),
    p = c("^p$", "^pval$", "^pvalue$", "^pvalue", "^pvalue\.?$"),
    n = c("^n$", "^totalsamplesize$", "^samplesize$", "^neff$"),
    ncase = c("^ncase$", "^cases$", "^ncases$"),
    ncontrol = c("^ncontrol$", "^controls$", "^ncontrols$"),
    info = c("^info$", "^imputationinfo$", "^infoscore$")
  )

  score_name <- function(field) {
    pats <- name_patterns[[field]]
    sapply(clean_names, function(nm) {
      if (any(grepl(paste(pats, collapse = "|"), nm, perl = TRUE))) 1 else 0
    })
  }

  get_col <- function(j) stat[[j]]

  is_numeric_col <- function(x) is.numeric(x) || is.integer(x)

  value_scores <- lapply(seq_along(col_names), function(j) {
    x <- get_col(j)
    x_char <- toupper(trimws(as.character(x)))
    x_num <- suppressWarnings(as.numeric(x))

    non_na_num <- x_num[!is.na(x_num)]
    non_na_char <- x_char[!is.na(x_char) & nzchar(x_char)]

    frac_numeric <- if (!length(x)) 0 else mean(!is.na(x_num))
    uniq <- length(unique(non_na_char))

    chr_score <- 0
    if (length(non_na_num)) {
      chr_score <- mean(non_na_num %in% c(1:26, 23, 24, 25, 26))
    }

    pos_score <- 0
    if (length(non_na_num)) {
      pos_score <- mean(non_na_num > 1000 & non_na_num < 4e8)
    }

    allele_score <- 0
    if (length(non_na_char)) {
      allele_score <- mean(non_na_char %in% c("A", "C", "G", "T"))
    }

    p_score <- 0
    if (length(non_na_num)) {
      p_score <- mean(non_na_num >= 0 & non_na_num <= 1)
      tiny <- mean(non_na_num < 1e-6)
      p_score <- p_score + 0.2 * tiny
    }

    eaf_score <- 0
    if (length(non_na_num)) {
      in01 <- mean(non_na_num >= 0 & non_na_num <= 1)
      med <- stats::median(non_na_num)
      shape <- as.numeric(med > 0.01 && med < 0.99)
      eaf_score <- in01 + 0.2 * shape
    }

    seb_score <- 0
    if (length(non_na_num)) {
      seb_score <- mean(non_na_num > 0 & non_na_num < 1)
    }

    b_score <- 0
    if (length(non_na_num)) {
      b_score <- as.numeric(any(non_na_num < 0)) + 0.5 * mean(abs(non_na_num) < 5)
    }

    n_score <- 0
    if (length(non_na_num)) {
      n_score <- mean(non_na_num >= 100)
    }

    info_score <- 0
    if (length(non_na_num)) {
      info_score <- mean(non_na_num >= 0 & non_na_num <= 1) +
        0.3 * as.numeric(stats::median(non_na_num) > 0.7)
    }

    marker_score <- as.numeric(!is_numeric_col(x)) * as.numeric(uniq > 1000)

    list(
      chr = chr_score,
      pos = pos_score,
      ea = allele_score,
      nea = allele_score,
      p = p_score,
      eaf = eaf_score,
      seb = seb_score,
      b = b_score,
      n = n_score,
      ncase = n_score,
      ncontrol = n_score,
      info = info_score,
      marker = marker_score,
      frac_numeric = frac_numeric,
      uniq = uniq
    )
  })

  fields <- c("marker", "chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n", "ncase", "ncontrol", "info")
  score_tbl <- matrix(0, nrow = length(fields), ncol = length(col_names),
                      dimnames = list(fields, col_names))

  for (f in fields) {
    nscore <- score_name(f)
    vscore <- vapply(value_scores, function(v) v[[f]], numeric(1))
    score_tbl[f, ] <- 0.65 * nscore + 0.35 * pmin(vscore, 1.25)
  }

  pick_best <- function(field, exclude = character()) {
    scores <- score_tbl[field, ]
    if (length(exclude)) scores[names(scores) %in% exclude] <- -Inf
    ord <- order(scores, decreasing = TRUE)
    best <- names(scores)[ord[1]]
    top <- scores[ord[1]]
    second <- if (length(ord) > 1) scores[ord[2]] else NA_real_
    confidence <- if (is.na(second) || is.infinite(second)) top else top - second
    list(best = if (is.finite(top) && top > 0.25) best else NA_character_,
         confidence = max(0, confidence),
         candidates = names(scores)[ord[seq_len(min(3, length(ord)))]] )
  }

  mapping <- setNames(rep(NA_character_, length(fields)), fields)
  confidence <- setNames(rep(0, length(fields)), fields)
  candidates <- vector("list", length(fields)); names(candidates) <- fields

  used <- character(0)
  for (f in c("marker", "chr", "pos", "p", "b", "seb", "eaf", "n", "ncase", "ncontrol", "info")) {
    p <- pick_best(f, exclude = used)
    mapping[f] <- p$best
    confidence[f] <- p$confidence
    candidates[[f]] <- p$candidates
    if (!is.na(p$best)) used <- c(used, p$best)
  }

  ea_pick <- pick_best("ea")
  nea_scores <- score_tbl["nea", ]
  if (!is.na(ea_pick$best)) nea_scores[ea_pick$best] <- -Inf
  nea_ord <- order(nea_scores, decreasing = TRUE)
  nea_best <- names(nea_scores)[nea_ord[1]]

  mapping["ea"] <- if (!is.na(ea_pick$best) && score_tbl["ea", ea_pick$best] > 0.25) ea_pick$best else NA_character_
  mapping["nea"] <- if (is.finite(nea_scores[nea_best]) && nea_scores[nea_best] > 0.25) nea_best else NA_character_
  confidence["ea"] <- ea_pick$confidence
  confidence["nea"] <- if (length(nea_ord) > 1) nea_scores[nea_ord[1]] - nea_scores[nea_ord[2]] else nea_scores[nea_ord[1]]
  candidates[["ea"]] <- ea_pick$candidates
  candidates[["nea"]] <- names(nea_scores)[nea_ord[seq_len(min(3, length(nea_ord)))]]

  warnings <- character(0)
  required <- c("chr", "pos", "ea", "nea", "p")
  missing_required <- required[is.na(mapping[required])]
  if (length(missing_required)) warnings <- c(warnings, paste("Missing required fields:", paste(missing_required, collapse = ", ")))
  low_conf <- names(confidence)[confidence < 0.1 & !is.na(mapping)]
  if (length(low_conf)) warnings <- c(warnings, paste("Low-confidence mappings:", paste(low_conf, collapse = ", ")))

  list(
    mapping = mapping,
    confidence = confidence,
    candidates = candidates,
    warnings = unique(warnings),
    details = as.data.frame(t(score_tbl), stringsAsFactors = FALSE)
  )
}

#' Normalize GWAS summary statistics to gact canonical column names
#'
#' @param stat A `data.frame` with raw summary statistics.
#' @param schema Output from [detectStatSchema()] or a named character vector.
#' @param drop_unmapped Logical; if `TRUE`, return only canonical + passthrough columns.
#' @param require_fields Canonical fields required to be present after normalization.
#' @return A normalized `data.frame` with canonical columns where possible.
#' @export
normalizeStatSchema <- function(stat,
                                schema,
                                drop_unmapped = FALSE,
                                require_fields = c("chr", "pos", "ea", "nea", "p")) {
  if (is.null(stat)) stop("`stat` is NULL")
  if (!is.data.frame(stat)) stat <- as.data.frame(stat)

  mapping <- if (is.list(schema) && !is.null(schema$mapping)) schema$mapping else schema
  if (is.null(names(mapping))) stop("`schema` mapping must be a named vector")

  out <- stat
  canon <- c("marker", "chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n", "ncase", "ncontrol", "info")

  for (nm in canon) {
    src <- mapping[[nm]]
    if (!is.null(src) && !is.na(src) && src %in% colnames(out)) {
      out[[nm]] <- out[[src]]
    }
  }

  if (!is.null(out$chr)) {
    chr <- as.character(out$chr)
    chr <- gsub("^chr", "", chr, ignore.case = TRUE)
    chr[chr %in% c("X", "x")] <- "23"
    chr[chr %in% c("Y", "y")] <- "24"
    chr[chr %in% c("MT", "M", "mt", "m")] <- "25"
    out$chr <- suppressWarnings(as.integer(chr))
  }

  if (!is.null(out$pos)) out$pos <- suppressWarnings(as.integer(out$pos))
  if (!is.null(out$ea)) out$ea <- toupper(as.character(out$ea))
  if (!is.null(out$nea)) out$nea <- toupper(as.character(out$nea))

  for (nm in c("eaf", "b", "seb", "p", "n", "ncase", "ncontrol", "info")) {
    if (!is.null(out[[nm]])) out[[nm]] <- suppressWarnings(as.numeric(out[[nm]]))
  }

  if (!is.null(out$p)) {
    out$p[out$p == 0] <- .Machine$double.xmin
  }

  if (is.null(out$marker) && "rsids" %in% colnames(out)) out$marker <- out$rsids
  if (is.null(out$marker) && "SNP" %in% colnames(out)) out$marker <- out$SNP

  missing <- require_fields[!require_fields %in% colnames(out) | sapply(require_fields, function(x) all(is.na(out[[x]])))]
  if (length(missing)) {
    stop("Missing required canonical fields after normalization: ", paste(missing, collapse = ", "))
  }

  if (drop_unmapped) {
    keep <- unique(c(canon, setdiff(colnames(out), names(mapping))))
    out <- out[, keep[keep %in% colnames(out)], drop = FALSE]
  }

  out
}

#' Validate normalized GWAS summary statistics
#'
#' @param stat A normalized `data.frame` (typically from [normalizeStatSchema()]).
#' @return A list with `ok`, `errors`, `warnings`, and `metrics`.
#' @export
validateStatSchema <- function(stat) {
  if (is.null(stat)) stop("`stat` is NULL")
  if (!is.data.frame(stat)) stat <- as.data.frame(stat)

  errors <- character(0)
  warnings <- character(0)

  required <- c("chr", "pos", "ea", "nea", "p")
  miss <- required[!required %in% colnames(stat)]
  if (length(miss)) errors <- c(errors, paste("Missing required columns:", paste(miss, collapse = ", ")))

  if ("p" %in% colnames(stat)) {
    bad_p <- sum(!(stat$p >= 0 & stat$p <= 1), na.rm = TRUE)
    if (bad_p > 0) errors <- c(errors, sprintf("Invalid p-values outside [0,1]: %s", bad_p))
  }

  if (all(c("ea", "nea") %in% colnames(stat))) {
    valid <- c("A", "C", "G", "T")
    bad_ea <- sum(!toupper(stat$ea) %in% valid, na.rm = TRUE)
    bad_nea <- sum(!toupper(stat$nea) %in% valid, na.rm = TRUE)
    if (bad_ea > 0) warnings <- c(warnings, sprintf("EA non-ACGT rows: %s", bad_ea))
    if (bad_nea > 0) warnings <- c(warnings, sprintf("NEA non-ACGT rows: %s", bad_nea))
  }

  if (all(c("chr", "pos", "ea", "nea") %in% colnames(stat))) {
    cpra <- paste(stat$chr, stat$pos, stat$ea, stat$nea, sep = "_")
    dup <- sum(duplicated(cpra))
    if (dup > 0) warnings <- c(warnings, sprintf("Duplicate chr-pos-allele rows: %s", dup))
  }

  metrics <- list(
    n_rows = nrow(stat),
    n_cols = ncol(stat),
    p_missing = if ("p" %in% colnames(stat)) mean(is.na(stat$p)) else NA_real_
  )

  list(ok = !length(errors), errors = unique(errors), warnings = unique(warnings), metrics = metrics)
}

#' Ingest GWAS summary statistics with AI-ready schema detection pipeline
#'
#' Runs schema detection, normalization, validation, and then delegates
#' storage/QC to [updateStatDB()].
#'
#' @param GAlist A gact database object.
#' @param stat A raw summary-statistics `data.frame`.
#' @param schema Optional schema output from [detectStatSchema()] or named mapping.
#' @param min_confidence Minimum per-field confidence used for auto mode.
#' @param strict If `TRUE`, stop when low confidence exists for required fields.
#' @param write_to_db Passed to [updateStatDB()] as `writeStatDB`.
#' @param ... Additional metadata arguments passed to [updateStatDB()] (e.g. `source`, `trait`, `type`).
#' @return A list with updated `GAlist`, `schema`, and `validation`.
#' @export
ingestStatDB <- function(GAlist,
                         stat,
                         schema = NULL,
                         min_confidence = 0.1,
                         strict = TRUE,
                         write_to_db = TRUE,
                         ...) {
  if (is.null(schema)) schema <- detectStatSchema(stat)

  mapping <- if (is.list(schema) && !is.null(schema$mapping)) schema$mapping else schema
  confidence <- if (is.list(schema) && !is.null(schema$confidence)) schema$confidence else setNames(rep(1, length(mapping)), names(mapping))

  required <- c("chr", "pos", "ea", "nea", "p")
  low_conf <- required[!is.na(mapping[required]) & confidence[required] < min_confidence]
  if (strict && length(low_conf)) {
    stop("Low-confidence required mappings: ", paste(low_conf, collapse = ", "),
         ". Provide `schema=` explicitly or set strict = FALSE.")
  }

  stat_norm <- normalizeStatSchema(stat, schema = mapping)
  validation <- validateStatSchema(stat_norm)
  if (!validation$ok) stop(paste(validation$errors, collapse = " | "))

  updated <- updateStatDB(
    GAlist = GAlist,
    stat = stat_norm,
    writeStatDB = write_to_db,
    ...
  )

  list(
    GAlist = updated,
    schema = schema,
    validation = validation,
    normalized_stat = stat_norm
  )
}
