scoreTemplateMatch <- function(colnames, template_mapping) {

 clean <- function(x) tolower(gsub("[^a-z0-9]", "", x))
 cols <- clean(colnames)

 score <- 0
 total <- 0

 for (f in names(template_mapping)) {

  aliases <- template_mapping[[f]]

  if (is.null(aliases) || all(is.na(aliases)) || length(aliases) == 0) next

  aliases <- clean(aliases)

  total <- total + 1

  if (any(aliases %in% cols)) {
   score <- score + 1
  }
 }

 if (total == 0) return(0)

 score / total
}

detectTemplate <- function(stat,
                           templates = c("UKB","FinnGen","SAIGE"),
                           threshold = 0.6) {

 coln <- colnames(stat)

 scores <- sapply(templates, function(tpl) {
  mapping <- getTemplateMapping(tpl)
  scoreTemplateMatch(coln, mapping)
 })

 if (all(scores == 0)) {
  return(list(template = NA, score = 0, all_scores = scores))
 }

 best <- names(which.max(scores))
 best_score <- max(scores)

 list(
  template = if (best_score >= threshold) best else NA,
  score = best_score,
  all_scores = scores
 )
}

resolveTemplateMapping <- function(stat, template) {

 tpl <- getTemplateMapping(template)
 cols <- colnames(stat)

 clean <- function(x) tolower(gsub("[^a-z0-9]", "", x))
 clean_cols <- clean(cols)

 mapping <- setNames(vector("list", length(tpl)), names(tpl))

 for (f in names(tpl)) {

  aliases <- tpl[[f]]

  if (is.null(aliases) || all(is.na(aliases))) {
   mapping[[f]] <- NA_character_
   next
  }

  aliases_clean <- clean(aliases)

  match_idx <- which(clean_cols %in% aliases_clean)

  mapping[[f]] <- if (length(match_idx)) cols[match_idx[1]] else NA_character_
 }

 mapping
}

autoDetectSchema <- function(stat,
                             template_threshold = 0.6,
                             verbose = TRUE) {

 tpl <- detectTemplate(stat, threshold = template_threshold)

 # -------------------------
 # TEMPLATE MODE
 # -------------------------
 if (!is.na(tpl$template)) {

  if (verbose) {
   message("Template detected: ", tpl$template,
           " (score=", round(tpl$score,2), ")")
  }

  mapping <- resolveTemplateMapping(stat, tpl$template)

  return(list(
   mode = "template",
   template = tpl$template,
   mapping = mapping,
   confidence = setNames(rep(tpl$score, length(mapping)), names(mapping)),
   schema = NULL
  ))
 }

 # -------------------------
 # AUTO MODE
 # -------------------------
 if (verbose) message("Falling back to schema detection")

 schema <- detectStatSchema(stat)

 return(list(
  mode = "auto",
  template = NA,
  mapping = schema$mapping,
  confidence = schema$confidence,
  schema = schema
 ))
}


ingestWithAuto <- function(GAlist, stat, ..., confidence_threshold = 0.15) {

 det <- autoDetectSchema(stat)

 mapping <- det$mapping
 confidence <- det$confidence

 # -------------------------
 # Confidence check
 # -------------------------
 low_conf <- names(mapping)[
  !is.na(mapping) & confidence < confidence_threshold
 ]

 if (length(low_conf)) {
  stop("Low confidence fields: ",
       paste(low_conf, collapse = ", "),
       ". Please review mapping.")
 }

 # -------------------------
 # Ingest
 # -------------------------
 ingestStatDB(
  GAlist = GAlist,
  stat = stat,
  schema = mapping,
  ...
 )
}


cleanRecipe <- function(recipe) {

 `%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || is.na(a) || identical(a, "")) b else a
 }

 recipe$source    <- recipe$source    %||% "unknown"
 recipe$ancestry  <- recipe$ancestry  %||% "unknown"
 recipe$build     <- recipe$build     %||% "GRCh37"
 recipe$gender    <- recipe$gender    %||% "both"
 recipe$comments  <- recipe$comments  %||% ""

 # Ensure mapping is named list
 if (!is.null(recipe$mapping)) {
  recipe$mapping <- as.list(recipe$mapping)
 }

 recipe
}



validateRecipe <- function(recipe, level = c("minimal","standard","strict")) {

 level <- match.arg(level)

 errors <- character()
 warnings <- character()
 info <- character()

 is_empty <- function(x) {
  is.null(x) || length(x) == 0 || all(is.na(x)) || identical(x, "")
 }

 # -----------------------------
 # REQUIRED
 # -----------------------------
 if (is_empty(recipe$trait)) {
  errors <- c(errors, "Missing trait")
 }

 if (is_empty(recipe$type)) {
  errors <- c(errors, "Missing type (quantitative/binary)")
 }

 # -----------------------------
 # Mapping checks
 # -----------------------------
 if (is.null(recipe$mapping)) {

  errors <- c(errors, "Missing mapping")

 } else {

  required_fields <- c("marker","chr","pos","ea","nea","b","p")

  for (f in required_fields) {

   val <- recipe$mapping[[f]]

   if (is_empty(val)) {
    errors <- c(errors, paste("Mapping missing for:", f))
   } else if (!is.character(val)) {
    errors <- c(errors, paste("Mapping not character for:", f))
   }
  }

  # duplicates
  mapped_cols <- unlist(recipe$mapping)
  mapped_cols <- mapped_cols[!is.na(mapped_cols)]

  dup <- mapped_cols[duplicated(mapped_cols)]

  if (length(dup)) {
   warnings <- c(warnings,
                 paste("Duplicate column mapping:", paste(unique(dup), collapse=", ")))
  }

  # recommended fields
  if (is_empty(recipe$mapping$seb)) {
   warnings <- c(warnings, "Missing SE (seb) → limits many models")
  }

  if (is_empty(recipe$mapping$eaf)) {
   warnings <- c(warnings, "Missing EAF → limits QC / LD-aware models")
  }
 }

 # -----------------------------
 # STANDARD
 # -----------------------------
 if (level %in% c("standard","strict")) {

  if (is_empty(recipe$source)) {
   warnings <- c(warnings, "Missing source → fallback used")
  }

  if (is_empty(recipe$ancestry)) {
   warnings <- c(warnings, "Missing ancestry → set to 'unknown'")
  }

  if (is_empty(recipe$build)) {
   warnings <- c(warnings, "Missing genome build")
  }

  if (is_empty(recipe$n) && is_empty(recipe$ncase)) {
   warnings <- c(warnings, "Missing sample size (n or ncase/ncontrol)")
  }
 }

 # -----------------------------
 # STRICT
 # -----------------------------
 if (level == "strict") {

  if (is_empty(recipe$category)) {
   warnings <- c(warnings, "Missing category")
  }

  if (is_empty(recipe$efo)) {
   warnings <- c(warnings, "Missing EFO term")
  }

  if (is_empty(recipe$reference)) {
   warnings <- c(warnings, "Missing reference (PMID)")
  }
 }

 # -----------------------------
 # Validation block
 # -----------------------------
 if (!is.null(recipe$validation)) {

  if (!isTRUE(recipe$validation$ok)) {
   errors <- c(errors, "Underlying data validation failed")
  }

  if (!is.null(recipe$validation$metrics$n_rows)) {
   if (recipe$validation$metrics$n_rows < 10) {
    warnings <- c(warnings, "Very small dataset preview (<10 rows)")
   }
  }
 }

 list(
  ok = length(errors) == 0,
  level = level,
  errors = unique(errors),
  warnings = unique(warnings),
  info = info
 )
}


# recipes <- list.files("recipes/", full.names = TRUE)
#
# for (f in recipes) {
#  r <- jsonlite::read_json(f)
#
#  v <- validateRecipe(r, "standard")
#
#  if (!v$ok) {
#   message("Skipping:", f)
#   next
#  }
#
#  runRecipe(r)
# }

runRecipe <- function(recipe, GAlist, db_path,
                      validate_level = "standard",
                      verbose = TRUE) {

 # -----------------------------
 # 1. Validate recipe
 # -----------------------------
 v <- validateRecipe(recipe, level = validate_level)

 if (!v$ok) {
  stop("Recipe failed validation:\n", paste(v$errors, collapse = "\n"))
 }

 if (length(v$warnings) > 0 && verbose) {
  message("Warnings:\n", paste(v$warnings, collapse = "\n"))
 }

 # -----------------------------
 # 2. Clean recipe
 # -----------------------------
 recipe <- cleanRecipe(recipe)

 # -----------------------------
 # 3. Load GWAS
 # -----------------------------
 has_file <- !is.null(recipe$file_name) &&
  nzchar(recipe$file_name) &&
  file.exists(recipe$file_name)

 has_url  <- !is.null(recipe$url) && nzchar(recipe$url)

 if (verbose) {
  message("Reading GWAS: ",
          if (has_file) recipe$file_name else recipe$url)
 }

 if (!has_file && has_url) {

  tmp <- tempfile(fileext = ".gz")

  tryCatch({
   download.file(recipe$url, tmp, mode = "wb")
  }, error = function(e) {
   stop("Failed to download GWAS from URL")
  })

  stat <- data.table::fread(tmp)

  # checksum
  recipe$file_md5 <- unname(tools::md5sum(tmp))

 } else if (has_file) {

  stat <- data.table::fread(recipe$file_name)
  recipe$file_md5 <- unname(tools::md5sum(recipe$file_name))

 } else {
  stop("No valid file or URL provided in recipe")
 }

 # -----------------------------
 # 4. Normalize
 # -----------------------------
 if (verbose) message("Normalizing...")

 stat_norm <- normalizeStatSchema(
  stat = stat,
  schema = recipe$mapping
 )

 # -----------------------------
 # 5. Validate stats
 # -----------------------------
 val <- validateStatSchema(stat_norm)

 if (!val$ok) {
  stop("Stat validation failed:\n", paste(val$errors, collapse = "\n"))
 }

 if (verbose) {
  message("Rows:", nrow(stat_norm))
 }

 # -----------------------------
 # 6. Ingest
 # -----------------------------
 if (verbose) message("Updating database...")

 G_updated <- updateStatDB(
  GAlist = GAlist,
  stat = stat_norm,
  source = recipe$source,
  trait = recipe$trait,
  type = recipe$type,
  category = recipe$category,
  efo = recipe$efo,
  gender = recipe$gender,
  ancestry = recipe$ancestry,
  build = recipe$build,
  n = recipe$n,
  ncase = recipe$ncase,
  ncontrol = recipe$ncontrol,
  reference = recipe$reference,
  comments = recipe$comments,
  recipe = recipe,
  writeStatDB = TRUE
 )

 # -----------------------------
 # 7. Save DB
 # -----------------------------
 saveRDS(G_updated, db_path, compress = FALSE)

 if (verbose) message("Done:", recipe$trait)

 return(G_updated)
}

createRecipesFromManifest <- function(meta, recipes_dir) {

 dir.create(recipes_dir, showWarnings = FALSE, recursive = TRUE)

 for (i in seq_len(nrow(meta))) {

  row <- meta[i, ]

  recipe <- list(
   file_name = row$filename,        # or NA if remote only
   url = row$aws_link,             # ⭐ important

   trait = row$Description,
   type = ifelse(row$trait_type == "biomarkers", "quantitative", "binary"),
   category = row$Category,
   efo = row$trait_efo_terms,

   gender = row$pheno_sex,
   ancestry = "EUR",               # or from column if available
   build = "GRCh37",

   n = row$n_cases_full_cohort_both_sexes,
   ncase = row$n_cases_EUR,
   ncontrol = row$n_controls_EUR,

   source = "UKB",
   reference = row$trait_efo_terms,
   comments = paste("phenocode:", row$phenocode),

   # ⚠️ mapping can be default or template-based
   mapping = list(
    marker = "rsid",
    chr = "#CHROM",
    pos = "POS",
    ea = "ALT",
    nea = "REF",
    eaf = "AF",
    b = "BETA",
    seb = "SE",
    p = "P",
    n = NA,
    ncase = NA,
    ncontrol = NA,
    info = NA
   ),

   created_at = as.character(Sys.time())
  )

  file_name <- file.path(
   recipes_dir,
   paste0(row$phenocode, ".json")
  )

  jsonlite::write_json(recipe, file_name, pretty = TRUE, auto_unbox = TRUE)
 }
}

validateAllRecipes <- function(recipes_dir) {

 files <- list.files(recipes_dir, "*.json", full.names = TRUE)

 for (f in files) {

  r <- jsonlite::read_json(f)
  v <- validateRecipe(r, "standard")

  if (!v$ok) {
   message("❌ Invalid:", f)
   print(v$errors)
  } else {
   message("✅ OK:", f)
  }
 }
}

runRecipes <- function(recipes_dir, GAlist, db_path) {

 files <- list.files(recipes_dir, "*.json", full.names = TRUE)

 for (f in files) {

  message("Processing:", f)

  recipe <- jsonlite::read_json(f)

  tryCatch({
   GAlist <- runRecipe(recipe, GAlist, db_path)
  }, error = function(e) {
   message("FAILED:", f)
  })
 }

 return(GAlist)
}

getTemplateMapping <- function(template = c("UKB","FinnGen","SAIGE")) {

 template <- match.arg(template)

 # canonical fields
 canon <- c(
  "marker","chr","pos","ea","nea","eaf",
  "b","seb","p","n","ncase","ncontrol","info"
 )

 mapping <- switch(template,

                   # -------------------------
                   # UK Biobank (Pan-UKB / BOLT / Neale)
                   # -------------------------
                   UKB = list(
                    marker   = c("rsid","variant","markername"),
                    chr      = c("#chrom","chr"),
                    pos      = c("pos","position"),
                    ea       = c("alt","allele1","effect_allele"),
                    nea      = c("ref","allele0","other_allele"),
                    eaf      = c("af","alt_af","eaf"),
                    b        = c("beta","b","effect"),
                    seb      = c("se","standard_error"),
                    p        = c("pval","p","p_bolt_lmm","p.value"),
                    n        = c("n","n_complete_samples"),
                    ncase    = c("ncase"),
                    ncontrol = c("ncontrol"),
                    info     = c("info")
                   ),

                   # -------------------------
                   # FinnGen
                   # -------------------------
                   FinnGen = list(
                    marker   = c("rsid"),
                    chr      = c("#chrom","chr"),
                    pos      = c("pos"),
                    ea       = c("alt"),
                    nea      = c("ref"),
                    eaf      = c("af_alt","eaf"),
                    b        = c("beta"),
                    seb      = c("sebeta","se"),
                    p        = c("pval","p"),
                    n        = c("n_total"),
                    ncase    = c("n_cases"),
                    ncontrol = c("n_controls"),
                    info     = c("info")
                   ),

                   # -------------------------
                   # SAIGE
                   # -------------------------
                   SAIGE = list(
                    marker   = c("markerid","snpid"),
                    chr      = c("chr"),
                    pos      = c("pos"),
                    ea       = c("allele2","alt"),
                    nea      = c("allele1","ref"),
                    eaf      = c("af","eaf"),
                    b        = c("beta"),
                    seb      = c("se"),
                    p        = c("p.value","pval"),
                    n        = c("n"),
                    ncase    = c("ncase"),
                    ncontrol = c("ncontrol"),
                    info     = c("info")
                   )
 )

 # ensure all canonical fields exist
 missing <- setdiff(canon, names(mapping))
 if (length(missing) > 0) {
  for (m in missing) mapping[[m]] <- NA
 }

 return(mapping)
}


#' Read GWAS alias dictionary from YAML
#'
#' @param path Optional path to YAML file. If `NULL`, uses
#'   `inst/extdata/gwas_aliases.yml` from the package.
#' @param must_exist Logical; if `TRUE`, throws an error when file is missing.
#' @return A normalized alias dictionary list.
#' @keywords internal
read_aliases_yaml <- function(path = NULL, must_exist = TRUE) {
  if (is.null(path)) {
    path <- system.file("extdata", "gwas_aliases.yml", package = "gact")
  }

  if (!nzchar(path) || !file.exists(path)) {
    if (must_exist) stop("Alias YAML not found: ", path)
    return(NULL)
  }

  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read alias dictionaries.")
  }

  x <- yaml::read_yaml(path)
  if (is.null(x$fields) || !is.list(x$fields)) {
    stop("Alias YAML must contain a top-level 'fields:' list.")
  }

  clean_name <- function(z) tolower(gsub("[^a-z0-9]", "", as.character(z)))
  `%||%` <- function(a, b) if (is.null(a)) b else a

  out <- list(version = x$version %||% NA, fields = list())
  for (fld in names(x$fields)) {
    node <- x$fields[[fld]]
    aliases <- unique(clean_name(node$aliases %||% character(0)))
    aliases <- aliases[nzchar(aliases)]
    regex <- unique(as.character(node$regex %||% character(0)))
    out$fields[[fld]] <- list(aliases = aliases, regex = regex)
  }
  out
}

#' Apply alias dictionary to column names and return alias score matrix
#'
#' @param col_names Character vector of column names from incoming data.
#' @param alias_dict Alias dictionary from [read_aliases_yaml()].
#' @param fields Canonical fields to score.
#' @return A list with `score` matrix (`field x column`) and `hits`.
#' @keywords internal
apply_alias_dictionary <- function(col_names, alias_dict, fields) {
  clean_name <- function(z) tolower(gsub("[^a-z0-9]", "", as.character(z)))
  clean_cols <- clean_name(col_names)

  score <- matrix(0, nrow = length(fields), ncol = length(col_names),
                  dimnames = list(fields, col_names))
  hits <- setNames(vector("list", length(fields)), fields)

  if (is.null(alias_dict) || is.null(alias_dict$fields)) {
    return(list(score = score, hits = hits))
  }

  for (f in fields) {
    fnode <- alias_dict$fields[[f]]
    if (is.null(fnode)) next

    alias_hits <- rep(FALSE, length(clean_cols))
    regex_hits <- rep(FALSE, length(clean_cols))

    if (length(fnode$aliases)) alias_hits <- clean_cols %in% fnode$aliases

    if (length(fnode$regex)) {
      regex_hits <- vapply(col_names, function(nm) {
        any(vapply(fnode$regex, function(rx) {
          #grepl(rx, nm, ignore.case = TRUE, perl = TRUE)
         clean_nm <- tolower(gsub("[^a-z0-9]", "", nm))
         grepl(rx, clean_nm, perl = TRUE)
        }, logical(1)))
      }, logical(1))
    }

    sc <- numeric(length(col_names))
    sc[alias_hits] <- pmax(sc[alias_hits], 1.0)
    sc[regex_hits] <- pmax(sc[regex_hits], 0.7)

    score[f, ] <- sc
    hits[[f]] <- col_names[alias_hits | regex_hits]
  }

  list(score = score, hits = hits)
}

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
#' @param alias_path Optional path to alias YAML. Defaults to package extdata file.
#' @param alias_dict Optional pre-loaded alias dictionary from [read_aliases_yaml()].
#' @param w_name Weight for name-pattern score.
#' @param w_value Weight for value-distribution score.
#' @param w_alias Weight for alias-dictionary score.
#' @return A list with elements `mapping`, `confidence`, `candidates`,
#'   `warnings`, `details`, and `alias_hits`.
#' @export

detectStatSchema <- function(stat, sample_n = 100000,
                             alias_path = NULL, alias_dict = NULL,
                             w_name = 0.45, w_value = 0.30, w_alias = 0.25) {
 if (is.null(stat)) stop("`stat` is NULL")
 if (!is.data.frame(stat)) stat <- as.data.frame(stat)
 if (!nrow(stat)) stop("`stat` has zero rows")

 if (!is.null(sample_n) && is.finite(sample_n) && nrow(stat) > sample_n) {
  idx <- sort(sample.int(nrow(stat), sample_n))
  stat <- stat[idx, , drop = FALSE]
 }

 col_names <- colnames(stat)
 clean_names <- tolower(gsub("[^a-z0-9]", "", col_names))

 # ---- FIX 1: load alias_dict BEFORE use ----
 if (is.null(alias_dict)) {
  alias_dict <- tryCatch(
   read_aliases_yaml(path = alias_path, must_exist = FALSE),
   error = function(e) NULL
  )
 }

 fields <- c("marker", "chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n", "ncase", "ncontrol", "info")

 # ---- FIX 2: safe name_patterns construction ----
 name_patterns <- setNames(vector("list", length(fields)), fields)

 if (!is.null(alias_dict) && !is.null(alias_dict$fields)) {
  for (f in intersect(names(alias_dict$fields), fields)) {
   fnode <- alias_dict$fields[[f]]
   if (length(fnode$aliases)) {
    name_patterns[[f]] <- paste0("^", fnode$aliases, "$")
   }
  }
 }

 # score_name <- function(field) {
 #  pats <- name_patterns[[field]]
 #  if (length(pats) == 0) return(rep(0, length(clean_names)))
 #  sapply(clean_names, function(nm) {
 #   if (any(grepl(paste(pats, collapse = "|"), nm, perl = TRUE))) 1 else 0
 #  })
 # }
 score_name <- function(field) {
  pats <- name_patterns[[field]]
  fnode <- alias_dict$fields[[field]]

  sapply(seq_along(col_names), function(i) {
   nm <- col_names[i]
   clean_nm <- clean_names[i]

   alias_hit <- FALSE
   if (length(pats)) {
    alias_hit <- any(grepl(paste(pats, collapse = "|"), clean_nm, perl = TRUE))
   }

   regex_hit <- FALSE
   if (!is.null(fnode$regex) && length(fnode$regex)) {
    regex_hit <- any(vapply(fnode$regex, function(rx) {
     grepl(rx, nm, ignore.case = TRUE, perl = TRUE) ||
      grepl(rx, clean_nm, perl = TRUE)
    }, logical(1)))
   }

   if (alias_hit) return(1.0)
   if (regex_hit) return(0.8)

   return(0)
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

  chr_score <- if (length(non_na_num)) mean(non_na_num %in% c(1:26, 23, 24, 25, 26)) else 0
  pos_score <- if (length(non_na_num)) mean(non_na_num > 1000 & non_na_num < 4e8) else 0

  allele_score <- if (length(non_na_char)) mean(grepl("^[ACGT]+$", non_na_char)) else 0

  # ---- FIX 3: better p detection ----
  p_score <- 0
  if (length(non_na_num)) {
   if (max(non_na_num, na.rm = TRUE) > 20 && median(non_na_num, na.rm = TRUE) > 5) {
    p_score <- 0.3  # likely -log10(p)
   } else {
    p_score <- mean(non_na_num >= 0 & non_na_num <= 1)
   }
  }

  eaf_score <- 0
  if (length(non_na_num)) {
   in01 <- mean(non_na_num >= 0 & non_na_num <= 1)
   med <- stats::median(non_na_num)
   shape <- as.numeric(med > 0.01 && med < 0.99)
   eaf_score <- in01 + 0.2 * shape
  }

  seb_score <- if (length(non_na_num)) mean(non_na_num > 0 & non_na_num < 1) else 0

  # ---- FIX 4: OR vs beta heuristic ----
  b_score <- 0
  if (length(non_na_num)) {
   if (all(non_na_num > 0, na.rm = TRUE) && median(non_na_num, na.rm = TRUE) > 0.5) {
    b_score <- 0.3  # OR
   } else if (any(non_na_num < 0, na.rm = TRUE)) {
    b_score <- 0.5  # beta
   }
  }

  n_score <- if (length(non_na_num)) mean(non_na_num >= 100) else 0

  info_score <- 0
  if (length(non_na_num)) {
   info_score <- mean(non_na_num >= 0 & non_na_num <= 1) +
    0.3 * as.numeric(stats::median(non_na_num) > 0.7)
  }

  # ---- FIX 5: marker_score restored ----
  marker_score <- as.numeric(!is_numeric_col(x)) * as.numeric(uniq > 1000)
  rs_frac <- if (length(non_na_char)) mean(grepl("^rs[0-9]+$", non_na_char)) else 0
  marker_score <- max(marker_score, rs_frac)

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

 score_tbl <- matrix(0, nrow = length(fields), ncol = length(col_names),
                     dimnames = list(fields, col_names))

 alias_res <- apply_alias_dictionary(col_names = col_names,
                                     alias_dict = alias_dict,
                                     fields = fields)
 alias_score_tbl <- alias_res$score

 for (f in fields) {
  nscore <- score_name(f)
  vscore <- vapply(value_scores, function(v) v[[f]], numeric(1))
  score_tbl[f, ] <- w_name * nscore +
   w_value * pmin(vscore, 1.25) +
   w_alias * alias_score_tbl[f, ]
 }

 warnings <- character(0)

 pick_best <- function(field, exclude = character()) {
  scores <- score_tbl[field, ]
  if (length(exclude)) scores[names(scores) %in% exclude] <- -Inf

  ord <- order(scores, decreasing = TRUE)
  top <- scores[ord[1]]
  second <- if (length(ord) > 1) scores[ord[2]] else NA_real_

  # ---- FIX 6: tie-breaking ----
  if (!is.na(second) && top == second) {
   alias_bonus <- alias_score_tbl[field, ord]
   ord <- ord[order(alias_bonus[ord], decreasing = TRUE)]
  }

  best <- names(scores)[ord[1]]
  confidence <- if (is.na(second) || is.infinite(second)) top else top - second
  ambiguous <- !is.na(second) && abs(top - second) < 0.05

  list(best = if (is.finite(top) && top > 0.25) best else NA_character_,
       confidence = max(0, confidence),
       candidates = names(scores)[ord[seq_len(min(3, length(ord)))]],
       ambiguous = ambiguous)
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
  if (p$ambiguous) warnings <- c(warnings, paste("Ambiguous mapping for", f))
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

 required <- c("chr", "pos", "ea", "nea", "p")
 missing_required <- required[is.na(mapping[required])]
 if (length(missing_required)) {
  warnings <- c(warnings, paste("Missing required fields:", paste(missing_required, collapse = ", ")))
 }

 low_conf <- names(confidence)[confidence < 0.1 & !is.na(mapping)]
 if (length(low_conf)) {
  warnings <- c(warnings, paste("Low-confidence mappings:", paste(low_conf, collapse = ", ")))
 }

 list(
  mapping = mapping,
  confidence = confidence,
  candidates = candidates,
  alias_hits = alias_res$hits,
  warnings = unique(warnings),
  details = as.data.frame(t(score_tbl), stringsAsFactors = FALSE)
 )
}

# detectStatSchema <- function(stat, sample_n = 100000,
#                              alias_path = NULL, alias_dict = NULL,
#                              w_name = 0.45, w_value = 0.30, w_alias = 0.25) {
#   if (is.null(stat)) stop("`stat` is NULL")
#   if (!is.data.frame(stat)) stat <- as.data.frame(stat)
#   if (!nrow(stat)) stop("`stat` has zero rows")
#
#   if (!is.null(sample_n) && is.finite(sample_n) && nrow(stat) > sample_n) {
#     idx <- sort(sample.int(nrow(stat), sample_n))
#     stat <- stat[idx, , drop = FALSE]
#   }
#
#   col_names <- colnames(stat)
#   clean_names <- tolower(gsub("[^a-z0-9]", "", col_names))
#
#   name_patterns <- lapply(alias_dict$fields, function(fnode) {
#    if (length(fnode$aliases)) {
#     paste0("^", fnode$aliases, "$")
#    } else character(0)
#   })
#
#   score_name <- function(field) {
#     pats <- name_patterns[[field]]
#     sapply(clean_names, function(nm) {
#       if (any(grepl(paste(pats, collapse = "|"), nm, perl = TRUE))) 1 else 0
#     })
#   }
#
#   get_col <- function(j) stat[[j]]
#
#   is_numeric_col <- function(x) is.numeric(x) || is.integer(x)
#
#   value_scores <- lapply(seq_along(col_names), function(j) {
#     x <- get_col(j)
#     x_char <- toupper(trimws(as.character(x)))
#     x_num <- suppressWarnings(as.numeric(x))
#
#     non_na_num <- x_num[!is.na(x_num)]
#     non_na_char <- x_char[!is.na(x_char) & nzchar(x_char)]
#
#     frac_numeric <- if (!length(x)) 0 else mean(!is.na(x_num))
#     uniq <- length(unique(non_na_char))
#
#     chr_score <- 0
#     if (length(non_na_num)) {
#       chr_score <- mean(non_na_num %in% c(1:26, 23, 24, 25, 26))
#     }
#
#     pos_score <- 0
#     if (length(non_na_num)) {
#       pos_score <- mean(non_na_num > 1000 & non_na_num < 4e8)
#     }
#
#     allele_score <- 0
#     if (length(non_na_char)) {
#       #allele_score <- mean(non_na_char %in% c("A", "C", "G", "T"))
#       allele_score <- mean(grepl("^[ACGT]+$", non_na_char))
#     }
#
#     p_score <- 0
#     #if (length(non_na_num)) {
#     #  p_score <- mean(non_na_num >= 0 & non_na_num <= 1)
#     #  tiny <- mean(non_na_num < 1e-6)
#     #  p_score <- p_score + 0.2 * tiny
#     #}
#     if (length(non_na_num)) {
#      if (max(non_na_num, na.rm = TRUE) > 50) {
#       # likely -log10(p)
#       p_score <- 0.3
#      } else {
#       p_score <- mean(non_na_num >= 0 & non_na_num <= 1)
#      }
#     }
#
#     eaf_score <- 0
#     if (length(non_na_num)) {
#       in01 <- mean(non_na_num >= 0 & non_na_num <= 1)
#       med <- stats::median(non_na_num)
#       shape <- as.numeric(med > 0.01 && med < 0.99)
#       eaf_score <- in01 + 0.2 * shape
#     }
#
#     seb_score <- 0
#     if (length(non_na_num)) {
#       seb_score <- mean(non_na_num > 0 & non_na_num < 1)
#     }
#
#     b_score <- 0
#     #if (length(non_na_num)) {
#     #  b_score <- as.numeric(any(non_na_num < 0)) + 0.5 * mean(abs(non_na_num) < 5)
#     #}
#     if (length(non_na_num)) {
#      if (all(non_na_num > 0, na.rm = TRUE) && median(non_na_num, na.rm = TRUE) > 0.5) {
#       # likely OR
#       b_score <- b_score + 0.3
#      } else if (any(non_na_num < 0, na.rm = TRUE)) {
#       # likely beta
#       b_score <- b_score + 0.5
#      }
#     }
#
#     n_score <- 0
#     if (length(non_na_num)) {
#       n_score <- mean(non_na_num >= 100)
#     }
#
#     info_score <- 0
#     if (length(non_na_num)) {
#       info_score <- mean(non_na_num >= 0 & non_na_num <= 1) +
#         0.3 * as.numeric(stats::median(non_na_num) > 0.7)
#     }
#
#     #marker_score <- as.numeric(!is_numeric_col(x)) * as.numeric(uniq > 1000)
#     rs_frac <- mean(grepl("^rs[0-9]+$", non_na_char))
#     marker_score <- max(marker_score, rs_frac)
#
#     list(
#       chr = chr_score,
#       pos = pos_score,
#       ea = allele_score,
#       nea = allele_score,
#       p = p_score,
#       eaf = eaf_score,
#       seb = seb_score,
#       b = b_score,
#       n = n_score,
#       ncase = n_score,
#       ncontrol = n_score,
#       info = info_score,
#       marker = marker_score,
#       frac_numeric = frac_numeric,
#       uniq = uniq
#     )
#   })
#
#   fields <- c("marker", "chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n", "ncase", "ncontrol", "info")
#   score_tbl <- matrix(0, nrow = length(fields), ncol = length(col_names),
#                       dimnames = list(fields, col_names))
#
#   if (is.null(alias_dict)) {
#     alias_dict <- tryCatch(
#       read_aliases_yaml(path = alias_path, must_exist = FALSE),
#       error = function(e) NULL
#     )
#   }
#   alias_res <- apply_alias_dictionary(col_names = col_names,
#                                       alias_dict = alias_dict,
#                                       fields = fields)
#   alias_score_tbl <- alias_res$score
#
#   for (f in fields) {
#     nscore <- score_name(f)
#     vscore <- vapply(value_scores, function(v) v[[f]], numeric(1))
#     score_tbl[f, ] <- w_name * nscore +
#       w_value * pmin(vscore, 1.25) +
#       w_alias * alias_score_tbl[f, ]
#   }
#
#   pick_best <- function(field, exclude = character()) {
#     scores <- score_tbl[field, ]
#     if (length(exclude)) scores[names(scores) %in% exclude] <- -Inf
#     ord <- order(scores, decreasing = TRUE)
#     best <- names(scores)[ord[1]]
#     top <- scores[ord[1]]
#     second <- if (length(ord) > 1) scores[ord[2]] else NA_real_
#     if (top == second) {
#      # prefer alias-supported column
#      alias_bonus <- alias_score_tbl[field, ord]
#      ord <- ord[order(alias_bonus[ord], decreasing = TRUE)]
#     }
#     confidence <- if (is.na(second) || is.infinite(second)) top else top - second
#     if (!is.na(second) && abs(top - second) < 0.05) {
#      warnings <- c(warnings, paste("Ambiguous mapping for", field))
#     }
#     list(best = if (is.finite(top) && top > 0.25) best else NA_character_,
#          confidence = max(0, confidence),
#          candidates = names(scores)[ord[seq_len(min(3, length(ord)))]] )
#   }
#
#   mapping <- setNames(rep(NA_character_, length(fields)), fields)
#   confidence <- setNames(rep(0, length(fields)), fields)
#   candidates <- vector("list", length(fields)); names(candidates) <- fields
#
#   used <- character(0)
#   for (f in c("marker", "chr", "pos", "p", "b", "seb", "eaf", "n", "ncase", "ncontrol", "info")) {
#     p <- pick_best(f, exclude = used)
#     mapping[f] <- p$best
#     confidence[f] <- p$confidence
#     candidates[[f]] <- p$candidates
#     if (!is.na(p$best)) used <- c(used, p$best)
#   }
#
#   ea_pick <- pick_best("ea")
#   nea_scores <- score_tbl["nea", ]
#   #if (!is.na(ea_pick$best)) nea_scores[ea_pick$best] <- -Inf
#   if (!is.na(ea_pick$best) && ea_pick$best == nea_best) {
#    nea_scores[nea_best] <- -Inf
#    nea_best <- names(sort(nea_scores, decreasing = TRUE))[1]
#   }
#   nea_ord <- order(nea_scores, decreasing = TRUE)
#   nea_best <- names(nea_scores)[nea_ord[1]]
#
#   mapping["ea"] <- if (!is.na(ea_pick$best) && score_tbl["ea", ea_pick$best] > 0.25) ea_pick$best else NA_character_
#   mapping["nea"] <- if (is.finite(nea_scores[nea_best]) && nea_scores[nea_best] > 0.25) nea_best else NA_character_
#   confidence["ea"] <- ea_pick$confidence
#   confidence["nea"] <- if (length(nea_ord) > 1) nea_scores[nea_ord[1]] - nea_scores[nea_ord[2]] else nea_scores[nea_ord[1]]
#   candidates[["ea"]] <- ea_pick$candidates
#   candidates[["nea"]] <- names(nea_scores)[nea_ord[seq_len(min(3, length(nea_ord)))]]
#
#   warnings <- character(0)
#   required <- c("chr", "pos", "ea", "nea", "p")
#   missing_required <- required[is.na(mapping[required])]
#   if (length(missing_required)) warnings <- c(warnings, paste("Missing required fields:", paste(missing_required, collapse = ", ")))
#   low_conf <- names(confidence)[confidence < 0.1 & !is.na(mapping)]
#   if (length(low_conf)) warnings <- c(warnings, paste("Low-confidence mappings:", paste(low_conf, collapse = ", ")))
#
#   list(
#     mapping = mapping,
#     confidence = confidence,
#     candidates = candidates,
#     alias_hits = alias_res$hits,
#     warnings = unique(warnings),
#     details = as.data.frame(t(score_tbl), stringsAsFactors = FALSE)
#   )
# }

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

 # ---- Standardize basic fields ----
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

 # ---- AUTO TRANSFORMATIONS ----
 transformations <- character(0)

 # Detect and convert -log10(p) → p
 if (!is.null(out$p)) {
  p_vals <- out$p[!is.na(out$p)]
  if (length(p_vals)) {
   if (max(p_vals, na.rm = TRUE) > 20 && mean(p_vals > 5, na.rm = TRUE) > 0.5) {
    out$p <- 10^(-out$p)
    transformations <- c(transformations, "-log10(p) -> p")
   }
  }
 }

 # Detect and convert OR → log(OR)
 if (!is.null(out$b)) {
  b_vals <- out$b[!is.na(out$b)]
  if (length(b_vals)) {
   if (all(b_vals > 0, na.rm = TRUE) && stats::median(b_vals, na.rm = TRUE) > 0.5) {
    out$b <- log(out$b)
    transformations <- c(transformations, "OR -> log(OR)")
   }
  }
 }

 # ---- Clamp extreme p-values ----
 if (!is.null(out$p)) {
  out$p[out$p == 0] <- .Machine$double.xmin
  out$p[out$p < .Machine$double.xmin] <- .Machine$double.xmin
 }

 # ---- Optional warning: SE mismatch ----
 if ("OR -> log(OR)" %in% transformations && !is.null(out$seb)) {
  warning("SE may not correspond to log(OR); verify input scale.")
 }

 # ---- Store transformation metadata ----
 if (length(transformations)) {
  attr(out, "transformations") <- transformations
 }

 # ---- Fallback marker handling ----
 if (is.null(out$marker) && "rsids" %in% colnames(out)) out$marker <- out$rsids
 if (is.null(out$marker) && "SNP" %in% colnames(out)) out$marker <- out$SNP

 # ---- Validate required fields ----
 missing <- require_fields[
  !require_fields %in% colnames(out) |
   sapply(require_fields, function(x) all(is.na(out[[x]])))
 ]
 if (length(missing)) {
  stop("Missing required canonical fields after normalization: ",
       paste(missing, collapse = ", "))
 }

 # ---- Optional column filtering ----
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
   normalized_stat = stat_norm,
   transformations = attr(stat_norm, "transformations")
  )
}
