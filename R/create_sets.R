#' Build all or selected GACT set modules
#'
#' @description
#' Creates and updates Ensembl, gene–marker, functional, drug, disease,
#' regulatory, and marker-level mapping sets.
#'
#' You can selectively build specific modules by providing `what`, e.g.:
#' ```
#' createSetsDB(GAlist, what = c("ensembl", "gene", "disease"))
#' ```
#'
#' @param GAlist GACT database list (from [createDB()])
#' @param what Character vector of modules to build. Options:
#'   `"ensembl"`, `"gene"`, `"functional"`, `"drug"`, `"disease"`,
#'   `"regulatory"`, `"marker"`. Defaults to all modules.
#' @param upstream Upstream distance (kb) for gene boundaries (default: 35)
#' @param downstream Downstream distance (kb) for gene boundaries (default: 10)
#' @param min_combined_score Minimum combined score for STRING/STITCH (default: 900)
#' @param min_interactions Minimum number of members per interaction set (default: 5)
#' @param overwrite Logical; if TRUE, rebuilds existing outputs
#' @return Updated GAlist with all relevant mappings
#' @export
createSetsDB <- function(GAlist,
                         what = c("ensembl", "gene", "functional",
                                  "drug", "disease", "regulatory", "marker"),
                         upstream = 35, downstream = 10,
                         min_combined_score = 900, min_interactions = 5,
                         overwrite = FALSE) {

 message("🧱 Building GACT sets for: ", paste(what, collapse = ", "))

 #------------------------------------------------------------
 # Normalize directories
 #------------------------------------------------------------
 for (nm in names(GAlist$dirs))
  GAlist$dirs[[nm]] <- normalizePath(GAlist$dirs[[nm]], winslash = "/", mustWork = FALSE)

 # Helper for checking required files
 check_files <- function(files, label) {
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
   message("⚠️ Skipping ", label, " — missing files:\n  ", paste(basename(missing), collapse = "\n  "))
   return(FALSE)
  }
  TRUE
 }

 #------------------------------------------------------------
 # Selective modular building
 #------------------------------------------------------------

 # ENSEMBL
 if ("ensembl" %in% what) {
  if (check_files(file.path(GAlist$dirs["ensembl"], "GRCh38.109.entrez.tsv.gz"), "Ensembl"))
   GAlist <- safe_call(createEnsemblMaps, GAlist=GAlist, overwrite=overwrite)
 }

 # GENE
 if ("gene" %in% what) {
  gtf_files <- list.files(GAlist$dirs["ensembl"], pattern = "\\.gtf\\.gz$", full.names = TRUE)
  marker_file <- file.path(GAlist$dirs["marker"], "markers.txt.gz")
  if (check_files(c(gtf_files, marker_file), "Gene–marker sets"))
   GAlist <- safe_call(createGeneMarkerSets, GAlist=GAlist, upstream=upstream,
                       downstream=downstream, overwrite=overwrite)
 }

 if ("string" %in% what)
  GAlist <- safe_call(createStringSets, GAlist=GAlist, min_combined_score=min_combined_score,
                      min_interactions=min_interactions, overwrite=overwrite)

 if ("stitch" %in% what)
  GAlist <- safe_call(createStitchSets, GAlist=GAlist, min_combined_score=min_combined_score,
                      min_interactions=min_interactions, overwrite=overwrite)

 # DRUG
 if ("drug" %in% what) {
  dgi_dir <- GAlist$dirs["drugdb"]
  if (check_files(file.path(dgi_dir, "interactions.tsv"), "Drug sets"))
   GAlist <- safe_call(createDrugSets, GAlist=GAlist, overwrite=overwrite)
 }

 # DISEASE
 if ("disease" %in% what) {
  disease_dir <- GAlist$dirs["diseases"]
  files <- list.files(disease_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
  if (length(files) > 0)
   GAlist <- safe_call(createDiseaseSets, GAlist=GAlist, overwrite=overwrite)
  else
   message("⚠️ Skipping Disease sets — no DISEASES files found.")
 }

 # REGULATORY
 if ("regulatory" %in% what) {
  reg_files <- list.files(GAlist$dirs["ensembl"], pattern = "Regulatory_Build.*\\.gff\\.gz$", full.names = TRUE)
  if (check_files(reg_files, "Regulatory sets")) {
   GAlist <- safe_call(createRegulatorySets, GAlist=GAlist, overwrite=overwrite)
   GAlist <- safe_call(createRegulatoryMarkerSets, GAlist, overwrite)
  }
 }

 # MARKER
 if ("marker" %in% what) {
  gsets_dir <- GAlist$dirs["gsets"]
  if (check_files(file.path(gsets_dir, "ensg2rsids.rds"), "Marker sets"))
   GAlist <- safe_call(createMarkerSetsDB, GAlist=GAlist, overwrite=overwrite)
 }

 if ("gwas" %in% what) {
  if (check_files(file.path(GAlist$dirs["gwascatalog"], "gwas-catalog-associations_ontology-annotated.tsv"), "GWAS Catalog"))
   GAlist <- safe_call(createGWASSets, GAlist=GAlist, overwrite=overwrite)
 }

 if ("pubmed" %in% what) {
  if (check_files(file.path(GAlist$dirs["gsets"], "gene2pubmed.gz"), "PubMed"))
   GAlist <- safe_call(createPubMedSets, GAlist=GAlist, overwrite=overwrite)
 }

 if ("gene_annotation" %in% what) {
  if (check_files(file.path(GAlist$dirs["gsets"], "ensg2rsids.rds"), "Gene annotation"))
   GAlist <- safe_call(createGeneAnnotationSets, GAlist=GAlist, overwrite=overwrite)
 }

 #------------------------------------------------------------
 # Metadata update (non-destructive)
 #------------------------------------------------------------
 if (is.null(GAlist$sets_metadata)) GAlist$sets_metadata <- list()
 GAlist$sets_metadata$created <- Sys.time()
 GAlist$sets_metadata$upstream <- upstream
 GAlist$sets_metadata$downstream <- downstream
 GAlist$sets_metadata$min_combined_score <- min_combined_score
 GAlist$sets_metadata$modules <- unique(GAlist$sets_metadata$modules)

 message("✅ GACT set construction completed successfully.")
 return(GAlist)
}



safe_call <- function(fun, ...) {
 args <- list(...)
 formals_names <- names(formals(fun))
 fun_name <- deparse(substitute(fun))

 # Always keep first positional argument (usually GAlist)
 first_arg <- if (length(args) > 0 && is.null(names(args)[1])) list(args[[1]]) else list()

 # Keep named arguments that match the function’s formals
 named_args <- args[names(args) %in% formals_names]

 # Combine them
 all_args <- c(first_arg, named_args)

 start_time <- Sys.time()
 result <- tryCatch(
  do.call(fun, all_args),
  error = function(e) {
   warning("❌ ", fun_name, " failed: ", e$message)
   if (length(args) > 0) return(args[[1]]) else return(NULL)
  }
 )
 end_time <- Sys.time()

 message("⏱ ", fun_name, " completed in ",
         round(difftime(end_time, start_time, units = "secs"), 2), " sec")
 return(result)
}



createEnsemblMaps <- function(GAlist, overwrite = FALSE) {
 file <- file.path(GAlist$dirs["ensembl"], "GRCh38.109.entrez.tsv.gz")
 if (!file.exists(file)) stop("Missing Ensembl annotation file: ", file)

 outfiles <- file.path(GAlist$dirs["gsets"], c(
  "eg2ensg.rds", "eg2ensp.rds", "eg2enst.rds",
  "ensg2eg.rds", "ensg2ensp.rds", "ensg2enst.rds",
  "ensp2eg.rds", "ensp2ensg.rds", "ensp2enst.rds"
 ))

 if (!overwrite && all(file.exists(outfiles))) {
  message("✅ Ensembl maps already exist — skipping.")
  return(GAlist)
 }

 message("📖 Reading Ensembl annotation")
 ensembl <- data.table::fread(file, data.table = FALSE)
 ensembl <- subset(ensembl, protein_stable_id != "-")
 ensembl$xref <- tolower(ensembl$xref)

 # Create mappings
 GAlist$map$eg2ensg  <- split(ensembl$gene_stable_id,  ensembl$xref)
 GAlist$map$eg2ensp  <- split(ensembl$protein_stable_id, ensembl$xref)
 GAlist$map$eg2enst  <- split(ensembl$transcript_stable_id, ensembl$xref)

 GAlist$map$ensp2eg  <- split(ensembl$xref, ensembl$protein_stable_id)
 GAlist$map$ensp2ensg <- split(ensembl$gene_stable_id, ensembl$protein_stable_id)
 GAlist$map$ensp2enst <- split(ensembl$transcript_stable_id, ensembl$protein_stable_id)

 GAlist$map$ensg2eg  <- split(ensembl$xref, ensembl$gene_stable_id)
 GAlist$map$ensg2ensp <- split(ensembl$protein_stable_id, ensembl$gene_stable_id)
 GAlist$map$ensg2enst <- split(ensembl$transcript_stable_id, ensembl$gene_stable_id)

 # Save all in parallel
 for (nm in names(GAlist$map)) {
  saveRDS(GAlist$map[[nm]], file.path(GAlist$dirs["gsets"], paste0(nm, ".rds")))
 }

 message("✅ Ensembl mappings created (", nrow(ensembl), " entries)")
 return(GAlist)
}

createGeneMarkerSets <- function(GAlist, upstream = 35, downstream = 10, overwrite = FALSE) {

 outpath <- file.path(GAlist$dirs["gsets"], "ensg2rsids.rds")
 if (!overwrite && file.exists(outpath)) {
  message("✅ Gene–marker sets already exist — skipping.")
  return(GAlist)
 }

 upstream <- upstream * 1000
 downstream <- downstream * 1000

 gtf <- file.path(GAlist$dirs["ensembl"], "GRCh38.109.gtf.gz")
 if (!file.exists(gtf)) stop("Missing GTF file: ", gtf)

 message("📖 Reading gene annotation from: ", gtf)
 df <- data.table::fread(gtf, skip = 1, header = FALSE)
 colnames(df) <- c("chr","source","type","start","end","score","strand","phase","attributes")

 df <- df[df$type == "gene" & df$source == "ensembl_havana", ]
 df$gene_id <- sub('.*gene_id "(.*?)";.*', '\\1', df$attributes)
 df$chr <- suppressWarnings(as.integer(df$chr))
 df <- subset(df, chr %in% 1:22)

 message("📖 Reading markers")
 markers <- data.table::fread(file.path(GAlist$dirs["marker"], "markers.txt.gz"))
 markers_dt <- data.table(chr = markers$chr, start = markers$pos, end = markers$pos, rsid = markers$rsids)

 genes_dt <- data.table(chr = df$chr,
                        start = pmax(1, df$start - upstream),
                        end = df$end + downstream,
                        gene_id = df$gene_id)

 data.table::setkey(markers_dt, chr, start, end)
 data.table::setkey(genes_dt, chr, start, end)

 message("⚙️ Mapping markers to genes using foverlaps()...")
 overlaps <- data.table::foverlaps(markers_dt, genes_dt, nomatch = 0L)
 ensg2rsids <- split(overlaps$rsid, overlaps$gene_id)
 saveRDS(ensg2rsids, outpath)

 message("✅ Gene–marker mapping done: ", length(ensg2rsids), " genes.")
 return(GAlist)
}

#' Create functional association sets (STRING, STITCH, Reactome)
#'
#' @param GAlist GACT database list
#' @param min_combined_score Minimum combined score to filter STRING/STITCH interactions
#' @param min_interactions Minimum number of members to keep each set
#' @param overwrite Logical; if TRUE, regenerate existing outputs
#'
#' @return Updated GAlist with functional sets saved in gsets directory
#' @export
createFunctionalSets <- function(GAlist,
                                 min_combined_score = 900,
                                 min_interactions = 5,
                                 overwrite = FALSE) {

 dirs <- GAlist$dirs
 map <- GAlist$map

 message("🧬 Creating functional sets (STRING, STITCH, Reactome)")

 # Ensure directories are normalized
 dirs <- lapply(dirs, normalizePath, winslash = "/", mustWork = FALSE)

 # ------------------------------------------------------------
 # 1️⃣ STRING: protein–protein interactions → ENSG sets
 # ------------------------------------------------------------
 string_file <- file.path(dirs[["string"]], "9606.protein.links.v12.0.txt.gz")
 out_string  <- file.path(dirs[["gsets"]], "string2ensg.rds")

 if (!file.exists(string_file)) {
  warning("STRING data missing: ", string_file)
 } else if (!overwrite && file.exists(out_string)) {
  message("✅ STRING sets already exist — skipping.")
 } else {
  message("📖 Reading STRING interactions: ", basename(string_file))
  dt <- data.table::fread(string_file, data.table = TRUE)
  dt[, protein1 := sub("9606\\.", "", protein1)]
  dt[, protein2 := sub("9606\\.", "", protein2)]
  dt <- dt[combined_score >= min_combined_score]

  dt <- dt[protein2 %in% names(map$ensp2ensg)]
  sets <- split(dt$protein2, f = as.factor(dt$protein1))
  sets <- sets[sapply(sets, length) >= min_interactions]

  # Map proteins to genes
  sets <- lapply(sets, function(x) na.omit(unlist(map$ensp2ensg[x])))
  saveRDS(sets, out_string)
  message("✅ STRING sets saved: ", length(sets), " proteins.")
 }

 # ------------------------------------------------------------
 # 2️⃣ STITCH: protein–chemical interactions → ENSG sets
 # ------------------------------------------------------------
 stitch_file <- file.path(dirs[["stitch"]], "9606.protein_chemical.links.v5.0.tsv.gz")
 out_stitch  <- file.path(dirs[["gsets"]], "stitch2ensg.rds")

 if (!file.exists(stitch_file)) {
  warning("STITCH data missing: ", stitch_file)
 } else if (!overwrite && file.exists(out_stitch)) {
  message("✅ STITCH sets already exist — skipping.")
 } else {
  message("📖 Reading STITCH interactions: ", basename(stitch_file))
  dt <- data.table::fread(stitch_file, data.table = TRUE)
  dt[, protein := sub("9606\\.", "", protein)]
  dt <- dt[protein %in% names(map$ensp2ensg)]
  dt <- dt[combined_score >= min_combined_score]

  sets <- split(dt$protein, f = as.factor(dt$chemical))
  sets <- sets[sapply(sets, length) >= min_interactions]

  sets <- lapply(sets, function(x) na.omit(unlist(map$ensp2ensg[x])))
  saveRDS(sets, out_stitch)
  message("✅ STITCH sets saved: ", length(sets), " chemicals.")
 }

 # ------------------------------------------------------------
 # 3️⃣ Reactome: pathways → ENSG sets
 # ------------------------------------------------------------
 reac_map_file <- file.path(dirs[["reactome"]], "Ensembl2Reactome.txt")
 reac_name_file <- file.path(dirs[["reactome"]], "ReactomePathways.txt")
 out_reac2ensg <- file.path(dirs[["gsets"]], "reactome2ensg.rds")
 out_ensg2reac <- file.path(dirs[["gsets"]], "ensg2reactome.rds")
 out_reac2names <- file.path(dirs[["gsets"]], "reactome2names.rds")

 if (!file.exists(reac_map_file) || !file.exists(reac_name_file)) {
  warning("Reactome data missing: ", reac_map_file)
 } else if (!overwrite && all(file.exists(c(out_reac2ensg, out_ensg2reac, out_reac2names)))) {
  message("✅ Reactome sets already exist — skipping.")
 } else {
  message("📖 Reading Reactome mappings")

  reactome <- data.table::fread(reac_map_file, data.table = FALSE, header = FALSE)
  reactome <- reactome[grep("R-HSA", reactome[[2]]), ]
  reac2ensg <- split(reactome[[1]], f = reactome[[2]])
  ensg2reac <- split(reactome[[2]], f = reactome[[1]])

  # Pathway names
  pathway <- data.table::fread(reac_name_file, data.table = FALSE, header = FALSE)
  pathway <- pathway[grep("R-HSA", pathway[[1]]), ]
  reac2names <- pathway[[2]]
  names(reac2names) <- pathway[[1]]

  saveRDS(reac2ensg,  out_reac2ensg)
  saveRDS(ensg2reac,  out_ensg2reac)
  saveRDS(reac2names, out_reac2names)

  message("✅ Reactome sets saved: ", length(reac2ensg), " pathways.")
 }

 # ------------------------------------------------------------
 # Metadata and return
 # ------------------------------------------------------------
 GAlist$sets_metadata$modules <- unique(c(GAlist$sets_metadata$modules, "functional"))
 GAlist
}


#' Create drug–gene and ATC mapping sets
#'
#' @description
#' Builds mappings between drugs, genes, and ATC codes from DGIdb and WHO ATC data.
#' Integrates HGNC identifiers to Ensembl Gene IDs and produces consistent mapping
#' objects for downstream analyses.
#'
#' @param GAlist GACT database list (output of [createDB()])
#' @param overwrite Logical; if TRUE, recompute even if .rds outputs already exist
#' @return Updated GAlist with drug-related mapping objects added and saved in gsets/
#' @export
createDrugSets <- function(GAlist, overwrite = FALSE) {

 dirs <- GAlist$dirs
 dirs <- lapply(dirs, normalizePath, winslash = "/", mustWork = FALSE)
 gsets_dir <- dirs[["gsets"]]
 drugdb_dir <- dirs[["drugdb"]]

 message("💊 Creating drug–gene–ATC mapping sets")

 #--------------------------------------------------
 # 1️⃣ Load DGIdb data
 #--------------------------------------------------
 dgi_files <- c("interactions.tsv", "genes.tsv", "drugs.tsv", "categories.tsv")
 dgi_paths <- file.path(drugdb_dir, dgi_files)
 if (!all(file.exists(dgi_paths))) {
  warning("DGIdb data missing: ", paste(dgi_files[!file.exists(dgi_paths)], collapse = ", "))
  return(GAlist)
 }

 message("📖 Reading DGIdb data")
 interactions <- data.table::fread(dgi_paths[1], quote = "\"", data.table = FALSE)
 genes_df <- data.table::fread(dgi_paths[2], quote = "\"", data.table = FALSE)
 drugs_df <- data.table::fread(dgi_paths[3], quote = "\"", data.table = FALSE)

 #--------------------------------------------------
 # 2️⃣ Load HGNC mapping (used to translate to ENSG)
 #--------------------------------------------------
 hgnc_url <- "https://ftp.ebi.ac.uk/pub/databases/genenames/out_of_date_hgnc/tsv/hgnc_complete_set.txt"
 message("📖 Reading HGNC annotation: ", basename(hgnc_url))
 #hgnc <- data.table::fread(hgnc_url, data.table = FALSE)
 hgnc <- tryCatch(
  data.table::fread(hgnc_url, data.table = FALSE),
  error = function(e) {
   warning("Failed to download HGNC mapping: ", e$message)
   return(NULL)
  }
 )
 if (is.null(hgnc)) return(GAlist)
 hgnc <- subset(hgnc, !is.na(ensembl_gene_id) & ensembl_gene_id != "")
 hgnc$hgnc_id <- tolower(hgnc$hgnc_id)
 hgnc2ensg <- setNames(hgnc$ensembl_gene_id, hgnc$hgnc_id)



 #--------------------------------------------------
 # 3️⃣ Build drug → gene mappings via HGNC
 #--------------------------------------------------
 message("⚙️ Mapping DGIdb drugs to HGNC/Ensembl")
 interactions$gene_concept_id <- tolower(interactions$gene_concept_id)

 drug2hgnc <- split(interactions$gene_concept_id, f = as.factor(interactions$drug_name))
 drug2ensg <- lapply(drug2hgnc, function(x) unique(na.omit(hgnc2ensg[x])))
 drug2ensg <- drug2ensg[sapply(drug2ensg, length) > 0]

 out_drug2ensg <- file.path(gsets_dir, "drug2ensg.rds")
 if (overwrite || !file.exists(out_drug2ensg)) {
  saveRDS(drug2ensg, out_drug2ensg)
  message("✅ Saved drug2ensg (", length(drug2ensg), " drugs).")
 } else {
  message("✅ drug2ensg already exists — skipping.")
 }

 #--------------------------------------------------
 # 4️⃣ Reverse mapping (ensg → drug)
 #--------------------------------------------------
 ensg <- unlist(drug2ensg, use.names = FALSE)
 drugs <- rep(names(drug2ensg), times = sapply(drug2ensg, length))
 ensg2drug <- split(drugs, f = as.factor(ensg))
 out_ensg2drug <- file.path(gsets_dir, "ensg2drug.rds")
 if (overwrite || !file.exists(out_ensg2drug)) saveRDS(ensg2drug, out_ensg2drug)

 #--------------------------------------------------
 # 5️⃣ Attach ATC codes (WHO classification)
 #--------------------------------------------------
 atc_dir <- dirs[["drugdb"]]
 atc_file <- file.path(atc_dir, "atc.rds")
 if (!file.exists(atc_file)) {
  warning("ATC file missing: ", atc_file)
  return(GAlist)
 }


 message("📖 Reading ATC codes")
 atc <- readRDS(atc_file)
 targets_df <- data.frame(Drug = drugs, Target = ensg, stringsAsFactors = FALSE)
 targets_df$ATC <- "Unknown"
 has_atc <- match(tolower(targets_df$Drug), tolower(atc$name))
 targets_df$ATC[!is.na(has_atc)] <- as.character(atc$code[has_atc[!is.na(has_atc)]])

 out_targets <- file.path(gsets_dir, "targets.rds")
 if (overwrite || !file.exists(out_targets)) saveRDS(targets_df, out_targets)

 #--------------------------------------------------
 # 6️⃣ Derive ATC-level groupings (L1–L4)
 #--------------------------------------------------
 message("⚙️ Creating ATC-level sets (L1–L4)")
 known <- subset(targets_df, ATC != "Unknown")
 drug2atc <- setNames(known$ATC, known$Drug)
 saveRDS(drug2atc, file.path(gsets_dir, "drug2atc.rds"))

 # Helper: create sets for each ATC level
 make_atc_sets <- function(level) {
  atc_trunc <- substr(drug2atc, 1, level)
  split(names(atc_trunc), f = atc_trunc)
 }

 for (lvl in 1:4) {
  out <- file.path(gsets_dir, sprintf("atcSets%d.rds", lvl))
  if (overwrite || !file.exists(out))
   saveRDS(make_atc_sets(c(1,3,4,5)[lvl]), out)
 }

 #--------------------------------------------------
 # 7️⃣ Update metadata
 #--------------------------------------------------
 GAlist$sets_metadata$modules <- unique(c(GAlist$sets_metadata$modules, "drug"))
 GAlist
}

#' Create disease–gene association sets
#'
#' @description
#' Processes JensenLab "DISEASES" files to build gene-level disease association sets.
#' Generates mappings between diseases, Ensembl proteins, and Ensembl genes.
#'
#' @param GAlist GACT database list (output of [createDB()])
#' @param overwrite Logical; if TRUE, rebuilds sets even if existing .rds files are present
#' @param min_score Optional numeric; filter associations by score column (if available)
#' @return Updated GAlist with disease mappings added to gsets/
#' @export
createDiseaseSets <- function(GAlist, overwrite = FALSE, min_score = NULL) {

 dirs <- GAlist$dirs
 dirs <- lapply(dirs, normalizePath, winslash = "/", mustWork = FALSE)
 gsets_dir <- dirs[["gsets"]]
 disease_dir <- dirs[["diseases"]]

 message("🧬 Creating disease–gene sets (JensenLab DISEASES)")

 #--------------------------------------------------
 # 1️⃣ Check available files
 #--------------------------------------------------
 patterns <- c(
  "human_disease_textmining_full.tsv.gz",
  "human_disease_textmining_filtered.tsv.gz",
  "human_disease_knowledge_full.tsv.gz",
  "human_disease_knowledge_filtered.tsv.gz",
  "human_disease_experiments_full.tsv.gz",
  "human_disease_experiments_filtered.tsv.gz",
  "human_disease_integrated_full.tsv.gz"
 )
 files <- file.path(disease_dir, patterns)
 existing <- files[file.exists(files)]

 if (length(existing) == 0L) {
  warning("No JensenLab disease files found in: ", disease_dir)
  return(GAlist)
 }

 message("📂 Found ", length(existing), " disease datasets.")

 #--------------------------------------------------
 # 2️⃣ Iterate over each file and build mappings
 #--------------------------------------------------
 ensp2ensg <- GAlist$map$ensp2ensg
 if (is.null(ensp2ensg)) stop("GAlist$map$ensp2ensg missing — run createEnsemblMaps() first.")

 for (file in existing) {
  base <- gsub(".tsv.gz", "", basename(file))
  out_ensp  <- file.path(gsets_dir, paste0("disease2ensp_", base, ".rds"))
  out_ensg  <- file.path(gsets_dir, paste0("disease2ensg_", base, ".rds"))
  out_ensg2 <- file.path(gsets_dir, paste0("ensg2disease_", base, ".rds"))

  if (!overwrite && all(file.exists(c(out_ensp, out_ensg, out_ensg2)))) {
   message("✅ Skipping existing disease set: ", base)
   next
  }

  message("📖 Reading: ", basename(file))
  df <- data.table::fread(file, data.table = FALSE)
  if (ncol(df) < 4) {
   warning("Unexpected format in: ", basename(file))
   next
  }

  # Optional filtering by score column
  if (!is.null(min_score) && "score" %in% tolower(names(df))) {
   score_col <- grep("score", names(df), ignore.case = TRUE, value = TRUE)[1]
   df <- subset(df, df[[score_col]] >= min_score)
  }

  # Filter for known Ensembl proteins
  df <- df[df[, 1] %in% names(ensp2ensg), ]
  colnames(df)[1:4] <- c("ensp", "taxid", "disease_id", "disease_name")

  #--------------------------------------------------
  # Disease → ENSP
  #--------------------------------------------------
  disease2ensp <- split(df$ensp, f = as.factor(df$disease_name))

  #--------------------------------------------------
  # Disease → ENSG (gene-level)
  #--------------------------------------------------
  message("⚙️ Mapping proteins to genes")
  disease2ensg <- lapply(disease2ensp, function(x) unique(unlist(ensp2ensg[x])))
  disease2ensg <- disease2ensg[sapply(disease2ensg, length) > 0]

  #--------------------------------------------------
  # ENSG → disease
  #--------------------------------------------------
  ensg <- unlist(disease2ensg, use.names = FALSE)
  diseases <- rep(names(disease2ensg), times = sapply(disease2ensg, length))
  ensg2disease <- split(diseases, f = as.factor(ensg))

  #--------------------------------------------------
  # Save results
  #--------------------------------------------------
  saveRDS(disease2ensp, out_ensp)
  saveRDS(disease2ensg, out_ensg)
  saveRDS(ensg2disease, out_ensg2)
  message("✅ Saved mappings for ", basename(file), " (",
          length(disease2ensg), " diseases).")
 }

 #--------------------------------------------------
 # 3️⃣ Update metadata
 #--------------------------------------------------
 GAlist$sets_metadata$modules <- unique(c(GAlist$sets_metadata$modules, "disease"))
 GAlist
}

createOMIMSets <- function(GAlist, overwrite = FALSE) {
 file <- file.path(GAlist$dirs["omim"], "mim2gene.txt.gz")
 if (!file.exists(file)) return(GAlist)
 dt <- data.table::fread(file, data.table = FALSE)
 sets <- split(dt$Ensembl_Gene, f = dt$Phenotype)
 saveRDS(sets, file.path(GAlist$dirs["gsets"], "omim2ensg.rds"))
 GAlist$sets_metadata$modules <- unique(c(GAlist$sets_metadata$modules, "omim"))
 GAlist
}

#' Create GWAS Catalog trait-to-marker sets
#'
#' @description
#' Builds SNP-level sets linking GWAS Catalog traits to rsIDs.
#' Requires the annotated GWAS Catalog TSV file.
#'
#' @param GAlist GACT database list
#' @param overwrite Logical; if TRUE, rebuild existing outputs
#' @return Updated GAlist with GWAS trait-level marker sets
#' @export
createGWASSets <- function(GAlist, overwrite = FALSE) {

 gwas_dir <- normalizePath(GAlist$dirs[["gwascatalog"]], winslash = "/", mustWork = FALSE)
 gsets_dir <- normalizePath(GAlist$dirs[["gsets"]], winslash = "/", mustWork = FALSE)

 gwas_file <- file.path(gwas_dir, "gwas-catalog-associations_ontology-annotated.tsv")
 out_file  <- file.path(gsets_dir, "gwas2rsids.rds")

 message("🧠 Creating GWAS Catalog trait-to-SNP sets")

 if (!file.exists(gwas_file)) {
  warning("Missing GWAS Catalog file: ", gwas_file)
  return(GAlist)
 }

 if (!overwrite && file.exists(out_file)) {
  message("✅ gwas2rsids.rds already exists — skipping.")
  return(GAlist)
 }

 message("📖 Reading GWAS Catalog associations")
 gwas <- data.table::fread(gwas_file, data.table = FALSE, quote = "")

 if (!("MAPPED_TRAIT" %in% names(gwas)) || !("SNPS" %in% names(gwas))) {
  stop("GWAS Catalog file does not contain expected columns: MAPPED_TRAIT, SNPS")
 }

 sets <- split(gwas$SNPS, f = gwas$MAPPED_TRAIT)
 sets <- lapply(sets, function(x) unique(na.omit(x)))

 saveRDS(sets, out_file)
 message("✅ Saved GWAS sets: ", length(sets), " traits.")

 GAlist$sets_metadata$modules <- unique(c(GAlist$sets_metadata$modules, "gwas"))
 return(GAlist)
}

#' Create PubMed–gene mapping sets
#'
#' @description
#' Creates mappings between Entrez Gene IDs and PubMed IDs using NCBI gene2pubmed data.
#'
#' @param GAlist GACT database list
#' @param overwrite Logical; if TRUE, rebuild existing outputs
#' @return Updated GAlist with PubMed–gene association sets
#' @export
createPubMedSets <- function(GAlist, overwrite = FALSE) {

 gsets_dir <- normalizePath(GAlist$dirs[["gsets"]], winslash = "/", mustWork = FALSE)
 pubmed_file <- file.path(gsets_dir, "gene2pubmed.gz")

 out_eg2pmid <- file.path(gsets_dir, "eg2pmid.rds")
 out_pmid2eg <- file.path(gsets_dir, "pmid2eg.rds")

 message("📚 Creating PubMed–gene mappings")

 if (!file.exists(pubmed_file)) {
  warning("Missing gene2pubmed.gz file: ", pubmed_file)
  return(GAlist)
 }

 if (!overwrite && all(file.exists(c(out_eg2pmid, out_pmid2eg)))) {
  message("✅ PubMed sets already exist — skipping.")
  return(GAlist)
 }

 message("📖 Reading NCBI gene2pubmed file")
 pubmed <- data.table::fread(pubmed_file, data.table = FALSE)
 if (ncol(pubmed) < 3) stop("Unexpected gene2pubmed format.")

 pubmed <- subset(pubmed, pubmed[, 1] == 9606)  # keep human
 pubmed <- pubmed[, -1, drop = FALSE]
 colnames(pubmed) <- c("GeneID", "PubMed_ID")

 eg2pmid <- split(pubmed$PubMed_ID, f = as.factor(pubmed$GeneID))
 pmid2eg <- split(pubmed$GeneID, f = as.factor(pubmed$PubMed_ID))

 saveRDS(eg2pmid, out_eg2pmid)
 saveRDS(pmid2eg, out_pmid2eg)

 message("✅ Saved PubMed–gene mappings (",
         length(eg2pmid), " Entrez genes, ",
         length(pmid2eg), " PubMed IDs).")

 GAlist$sets_metadata$modules <- unique(c(GAlist$sets_metadata$modules, "pubmed"))
 return(GAlist)
}

#' Create gene coordinate annotation summary
#'
#' @description
#' Summarizes chromosomal positions (chr, start, end, length) for Ensembl gene sets.
#' Equivalent to `*_genesplus_annotation.rds` in the old pipeline.
#'
#' @param GAlist GACT database list
#' @param overwrite Logical; if TRUE, rebuild even if existing outputs are found
#' @return Updated GAlist with coordinate summary annotation
#' @export
createGeneAnnotationSets <- function(GAlist, overwrite = FALSE) {

 gsets_dir <- normalizePath(GAlist$dirs[["gsets"]], winslash = "/", mustWork = FALSE)
 marker_file <- file.path(GAlist$dirs[["marker"]], "markers.txt.gz")
 ensg2rsids_file <- file.path(gsets_dir, "ensg2rsids.rds")
 out_file <- file.path(gsets_dir, "genesplus_annotation.rds")

 message("🧭 Creating gene coordinate annotation summary")

 if (!file.exists(marker_file) || !file.exists(ensg2rsids_file)) {
  warning("Missing required input files for gene annotation.")
  return(GAlist)
 }

 if (!overwrite && file.exists(out_file)) {
  message("✅ Gene annotation already exists — skipping.")
  return(GAlist)
 }

 markers <- data.table::fread(marker_file, data.table = FALSE)
 ensg2rsids <- readRDS(ensg2rsids_file)

 sets <- mapSetsDB(sets = ensg2rsids, featureID = markers$rsids, index = TRUE)

 chr   <- sapply(sets, function(x) unique(markers$chr[x]))
 start <- sapply(sets, function(x) min(markers$pos[x]))
 stop  <- sapply(sets, function(x) max(markers$pos[x]))

 df <- data.frame(EnsemblID = names(sets), chr, start, stop)
 df$length <- df$stop - df$start

 saveRDS(df, out_file)
 message("✅ Saved gene coordinate annotation (", nrow(df), " genes).")

 GAlist$sets_metadata$modules <- unique(c(GAlist$sets_metadata$modules, "gene_annotation"))
 return(GAlist)
}


#' Create STRING protein–protein interaction sets
#'
#' @description
#' Builds both forward and reverse gene-level sets from STRING protein–protein
#' interaction data.
#'
#' @param GAlist GACT database list
#' @param min_combined_score Minimum combined score threshold (default: 900)
#' @param min_interactions Minimum number of genes per set (default: 5)
#' @param overwrite Logical; if TRUE, overwrite existing sets
#' @return Updated GAlist with STRING-based sets
#' @export
createStringSets <- function(GAlist,
                             min_combined_score = 900,
                             min_interactions = 5,
                             overwrite = FALSE) {

 dirs <- lapply(GAlist$dirs, normalizePath, winslash = "/", mustWork = FALSE)
 map  <- GAlist$map

 string_file <- file.path(dirs[["string"]], "9606.protein.links.v12.0.txt.gz")
 out_forward <- file.path(dirs[["gsets"]], "string2ensg.rds")
 out_reverse <- file.path(dirs[["gsets"]], "ensg2string.rds")

 message("🧠 Creating STRING protein–protein sets")

 if (!file.exists(string_file)) {
  warning("❌ Missing STRING file: ", string_file)
  return(GAlist)
 }

 if (!overwrite && all(file.exists(c(out_forward, out_reverse)))) {
  message("✅ STRING sets already exist — skipping.")
  return(GAlist)
 }

 dt <- fread(string_file, data.table=FALSE)
 dt$protein1 <- gsub("9606.","",dt$protein1)
 dt$protein2 <- gsub("9606.","",dt$protein2)
 dt  <- dt[dt$combined_score>=min_combined_score,]
 dt <- dt[dt$protein2%in%names(GAlist$map$ensp2ensg),]
 dt <- split( dt$protein2,f=as.factor(dt$protein1))
 string2ensg <- dt[sapply(dt ,length)>=min_interactions]
 string2ensg <- lapply(string2ensg,function(x){na.omit(unlist(GAlist$map$ensp2ensg[x]))})

 saveRDS(string2ensg, out_forward)

 message("✅ STRING sets saved :",length(string2ensg))

 GAlist$sets_metadata$modules <- unique(c(GAlist$sets_metadata$modules, "string"))
 return(GAlist)
}

#' Create STITCH protein–chemical interaction sets
#'
#' @description
#' Builds both forward (chemical → gene) and reverse (gene → chemical)
#' mappings from STITCH interaction data.
#'
#' @param GAlist GACT database list
#' @param min_combined_score Minimum combined score to filter STITCH interactions
#' @param min_interactions Minimum number of proteins per chemical (default: 5)
#' @param overwrite Logical; if TRUE, overwrite existing sets
#' @return Updated GAlist with STITCH sets
#' @export
createStitchSets <- function(GAlist,
                             min_combined_score = 900,
                             min_interactions = 5,
                             overwrite = FALSE) {

 dirs <- lapply(GAlist$dirs, normalizePath, winslash = "/", mustWork = FALSE)
 map  <- GAlist$map

 stitch_file <- file.path(dirs[["stitch"]], "9606.protein_chemical.links.v5.0.tsv.gz")
 out_forward <- file.path(dirs[["gsets"]], "stitch2ensg.rds")
 out_reverse <- file.path(dirs[["gsets"]], "ensg2stitch.rds")

 message("💊 Creating STITCH chemical–protein sets")

 if (!file.exists(stitch_file)) {
  warning("❌ Missing STITCH file: ", stitch_file)
  return(GAlist)
 }

 if (!overwrite && all(file.exists(c(out_forward, out_reverse)))) {
  message("✅ STITCH sets already exist — skipping.")
  return(GAlist)
 }

 stitch <- fread(stitch_file, data.table=FALSE)
 stitch$protein <- gsub("9606.","",stitch$protein)
 stitch <- stitch[stitch$protein%in%names(GAlist$map$ensp2ensg),]
 stitch  <- stitch[stitch$combined_score>=min_combined_score,]
 stitch <- split( stitch$protein,f=as.factor(stitch$chemical))
 sets  <- stitch[sapply(stitch ,length)>=min_interactions]
 stitch2ensg <- lapply(sets,function(x){na.omit(unlist(GAlist$map$ensp2ensg[x]))})

 saveRDS(stitch2ensg, out_forward)
 message("✅ STITCH:",length(stitch2ensg))

 GAlist$sets_metadata$modules <- unique(c(GAlist$sets_metadata$modules, "stitch"))
 return(GAlist)
}



#' Create marker-level (SNP-level) sets from gene-level annotations
#'
#' @description
#' Converts gene-based functional sets (GO, Reactome, STRING, STITCH, Drug, etc.)
#' into marker-level sets using `ensg2rsids` mappings.
#'
#' @param GAlist GACT database list (output of [createDB()])
#' @param what Character vector indicating which sets to process:
#'   `"GO"`, `"reactome"`, `"string"`, `"stitch"`, `"drug"`, or others.
#' @param overwrite Logical; if TRUE, recompute existing sets.
#' @return Updated GAlist with SNP-level sets added to `gsets/`
#' @export
createMarkerSetsDB <- function(GAlist,
                               what = NULL,
                               overwrite = FALSE) {

 dirs <- GAlist$dirs
 dirs <- lapply(dirs, normalizePath, winslash = "/", mustWork = FALSE)
 gsets_dir <- dirs[["gsets"]]

 #--------------------------------------------------
 # Helper function: map gene-based sets to marker-based
 #--------------------------------------------------
 mapGeneSetsToMarkers <- function(file_in, file_out, ensg2rsids) {
  if (!file.exists(file_in)) {
   warning("Missing input file: ", file_in)
   return(NULL)
  }
  if (!overwrite && file.exists(file_out)) {
   message("✅ ", basename(file_out), " already exists — skipping.")
   return(NULL)
  }
  message("📖 Reading: ", basename(file_in))
  fset <- readRDS(file_in)
  fset <- mapSetsDB(fset, featureID = names(ensg2rsids))
  sets <- lapply(fset, function(x) unique(unlist(ensg2rsids[x])))
  sets <- sets[!sapply(sets, is.null)]
  saveRDS(sets, file_out)
  message("✅ Saved ", basename(file_out), " (", length(sets), " sets).")
  invisible(sets)
 }

 #--------------------------------------------------
 # Load ensg2rsids (required)
 #--------------------------------------------------
 ensg2rsids_file <- file.path(gsets_dir, "ensg2rsids.rds")
 if (!file.exists(ensg2rsids_file))
  stop("Missing ensg2rsids file — run createGeneMarkerSets() first.")
 ensg2rsids <- readRDS(ensg2rsids_file)

 #--------------------------------------------------
 # 1️⃣ GO → rsIDs
 #--------------------------------------------------
 if ("GO" %in% what) {
  mapGeneSetsToMarkers(
   file_in  = file.path(gsets_dir, "go.rds"),
   file_out = file.path(gsets_dir, "go2rsids.rds"),
   ensg2rsids = ensg2rsids
  )
 }

 #--------------------------------------------------
 # 2️⃣ Reactome → rsIDs
 #--------------------------------------------------
 if ("reactome" %in% what) {
  mapGeneSetsToMarkers(
   file_in  = file.path(gsets_dir, "reactome2ensg.rds"),
   file_out = file.path(gsets_dir, "reactome2rsids.rds"),
   ensg2rsids = ensg2rsids
  )
 }

 #--------------------------------------------------
 # 3️⃣ STRING → rsIDs
 #--------------------------------------------------
 if ("string" %in% what) {
  mapGeneSetsToMarkers(
   file_in  = file.path(gsets_dir, "string2ensg.rds"),
   file_out = file.path(gsets_dir, "string2rsids.rds"),
   ensg2rsids = ensg2rsids
  )
 }

 #--------------------------------------------------
 # 4️⃣ STITCH → rsIDs
 #--------------------------------------------------
 if ("stitch" %in% what) {
  mapGeneSetsToMarkers(
   file_in  = file.path(gsets_dir, "stitch2ensg.rds"),
   file_out = file.path(gsets_dir, "stitch2rsids.rds"),
   ensg2rsids = ensg2rsids
  )
 }

 #--------------------------------------------------
 # 5️⃣ Drug sets → rsIDs
 #--------------------------------------------------
 if ("drug" %in% what) {

  message("💊 Creating drug-level SNP sets")

  drug2ensg_file <- file.path(gsets_dir, "drug2ensg.rds")
  string2ensg_file <- file.path(gsets_dir, "string2ensg.rds")

  if (!file.exists(drug2ensg_file) || !file.exists(string2ensg_file)) {
   warning("Missing drug or STRING mapping files.")
   return(GAlist)
  }

  drug2ensg <- readRDS(drug2ensg_file)
  string2ensg <- readRDS(string2ensg_file)

  # Extend drug targets with STRING network neighbours
  drug2ensp <- lapply(drug2ensg, function(x)
   na.omit(unlist(GAlist$map$ensg2ensp[x]))
  )
  drug2string2ensg <- lapply(drug2ensp, function(x)
   unique(na.omit(unlist(string2ensg[x])))
  )
  for (i in seq_along(drug2string2ensg)) {
   drug2string2ensg[[i]] <- unique(c(drug2ensg[[i]], drug2string2ensg[[i]]))
  }

  # Save intermediate
  saveRDS(drug2string2ensg, file.path(gsets_dir, "drug2string2ensg.rds"))

  # Drug → SNP
  drug2rsids <- lapply(drug2ensg, function(x) unique(unlist(ensg2rsids[x])))
  drug2rsids <- drug2rsids[!sapply(drug2rsids, is.null)]
  saveRDS(drug2rsids, file.path(gsets_dir, "drug2rsids.rds"))

  # Drug + STRING neighbours → SNP
  drug2string2rsids <- lapply(drug2string2ensg, function(x)
   unique(unlist(ensg2rsids[x]))
  )
  drug2string2rsids <- drug2string2rsids[!sapply(drug2string2rsids, is.null)]
  saveRDS(drug2string2rsids, file.path(gsets_dir, "drug2string2rsids.rds"))

  message("✅ Saved drug2rsids and drug2string2rsids")
 }

 #--------------------------------------------------
 # Metadata and return
 #--------------------------------------------------
 GAlist$sets_metadata$modules <- unique(c(GAlist$sets_metadata$modules, "marker"))
 GAlist
}

#' Create regulatory element mappings (Ensembl Regulatory Build)
#'
#' @description
#' Builds mappings between regulatory elements (ENSRs), genes (ENSGs), and SNPs (rsIDs)
#' from Ensembl Regulatory Build annotations.
#'
#' @param GAlist GACT database list (output of [createDB()])
#' @param build Genome build: "GRCh37" or "GRCh38" (default: "GRCh37")
#' @param overwrite Logical; if TRUE, recompute even if outputs exist
#' @return Updated GAlist with regulatory mappings added to gsets/
#' @export
createRegulatorySets <- function(GAlist,
                                 build = "GRCh37",
                                 overwrite = FALSE) {

  dirs <- GAlist$dirs
  dirs <- lapply(dirs, normalizePath, winslash = "/", mustWork = FALSE)
  gsets_dir <- dirs[["gsets"]]
  ensembl_dir <- dirs[["ensembl"]]
  marker_dir <- dirs[["marker"]]

  message("🧬 Creating regulatory element mappings (", build, ")")

  #--------------------------------------------------
  # Input checks
  #--------------------------------------------------
  gff_file <- file.path(ensembl_dir,
    paste0(build, ".Regulatory_Build.regulatory_features.gff.gz"))
  if (!file.exists(gff_file))
    stop("Missing regulatory GFF file: ", gff_file)

  markers_file <- file.path(marker_dir, "markers.txt.gz")
  if (!file.exists(markers_file))
    stop("Missing marker file: ", markers_file)

  ensg2rsids_file <- file.path(gsets_dir, "ensg2rsids.rds")
  if (!file.exists(ensg2rsids_file))
    stop("Missing ensg2rsids.rds — run createGeneMarkerSets() first.")
  ensg2rsids <- readRDS(ensg2rsids_file)

  #--------------------------------------------------
  # Read regulatory annotation
  #--------------------------------------------------
  df <- data.table::fread(gff_file, data.table = FALSE)
  colnames(df) <- c("chr","source","type","start","end","score","strand","phase","attributes")
  df <- subset(df, chr %in% as.character(1:22))
  df$chr <- as.integer(df$chr)
  df <- df[!is.na(df$chr), ]

  # Extract regulatory IDs (ENSR)
  att <- strsplit(df$attributes, ";")
  att <- lapply(att, function(x) gsub("\"", "", x))
  att <- sapply(att, function(x) x[grep("ID=", x)])
  att <- strsplit(att, ":")
  df$reg_id <- sapply(att, function(x) x[2])
  rownames(df) <- df$reg_id

  # Save base annotation
  saveRDS(df[, c("reg_id", "type", "chr", "start", "end")],
          file = file.path(gsets_dir, "regulatory_annotation.rds"))

  #--------------------------------------------------
  # Build type-based and SNP-based mappings
  #--------------------------------------------------
  reg2ensr <- split(df$reg_id, f = as.factor(df$type))
  saveRDS(reg2ensr, file = file.path(gsets_dir, "reg2ensr.rds"))

  message("⚙️ Mapping regulatory elements to SNPs ...")
  markers <- data.table::fread(markers_file, data.table = FALSE)
  maxpos <- max(markers$pos, df$end)
  ensr2rsids <- vector("list", nrow(df))

  #for (chr in 1:22) {
  #  chr_idx <- which(df$chr == chr)
  #  if (length(chr_idx) == 0) next
  #  rsids_chr <- rep(NA, maxpos)
  #  pos_chr <- markers$pos[markers$chr == chr]
  #  rs_chr <- markers$rsids[markers$chr == chr]
  #  rsids_chr[pos_chr] <- rs_chr
  #  for (i in chr_idx) {
  #    region <- rsids_chr[df$start[i]:df$end[i]]
  #    ensr2rsids[[i]] <- region[!is.na(region)]
  #  }
  #}

  for (chr in 1:22) {
   chr_df <- df[df$chr == chr, ]
   chr_markers <- markers[markers$chr == chr, ]
   if (nrow(chr_df) == 0L || nrow(chr_markers) == 0L) next

   for (i in seq_len(nrow(chr_df))) {
    range_idx <- which(chr_markers$pos >= chr_df$start[i] &
                        chr_markers$pos <= chr_df$end[i])
    if (length(range_idx)) {
     ensr2rsids[[chr_df$reg_id[i]]] <- chr_markers$rsids[range_idx]
    }
   }
   message("Completed creating regulatory element mappings for chromosome(", chr, ")")
  }


  names(ensr2rsids) <- df$reg_id
  ensr2rsids <- ensr2rsids[!sapply(ensr2rsids, function(x) identical(x, character(0)))]
  saveRDS(ensr2rsids, file = file.path(gsets_dir, "ensr2rsids.rds"))

  # Type-level sets: reg2rsids
  reg2rsids <- sapply(reg2ensr, function(x) unique(unlist(ensr2rsids[x])))
  saveRDS(reg2rsids, file = file.path(gsets_dir, "reg2rsids.rds"))

  #--------------------------------------------------
  # Map ENSR ↔ ENSG via shared SNPs
  #--------------------------------------------------
  message("🔗 Linking regulatory features ↔ genes")
  ensr_df <- data.frame(
    rsids = unlist(ensr2rsids, use.names = FALSE),
    ensr = rep(names(ensr2rsids), times = sapply(ensr2rsids, length))
  )
  ensg_df <- data.frame(
    rsids = unlist(ensg2rsids, use.names = FALSE),
    ensg = rep(names(ensg2rsids), times = sapply(ensg2rsids, length))
  )
  merged_df <- merge(ensr_df, ensg_df, by = "rsids")

  ensr2ensg <- split(merged_df$ensr, merged_df$ensg)
  ensr2ensg <- lapply(ensr2ensg, unique)
  ensr2ensg <- ensr2ensg[!sapply(ensr2ensg, function(x) identical(x, character(0)))]
  saveRDS(ensr2ensg, file = file.path(gsets_dir, "ensr2ensg.rds"))

  ensg2ensr <- split(merged_df$ensg, merged_df$ensr)
  ensg2ensr <- lapply(ensg2ensr, unique)
  ensg2ensr <- ensg2ensr[!sapply(ensg2ensr, function(x) identical(x, character(0)))]
  saveRDS(ensg2ensr, file = file.path(gsets_dir, "ensg2ensr.rds"))

  GAlist$sets_metadata$modules <-
    unique(c(GAlist$sets_metadata$modules, "regulatory"))

  message("✅ Regulatory sets created (",
          length(ensr2rsids), " ENSR elements).")
  GAlist
}


#' Create SNP-level regulatory marker sets
#'
#' @description
#' Converts regulatory element mappings (ENSR → ENSG) into SNP-level sets (ENSR → rsIDs).
#' Works within the same architecture as [createMarkerSetsDB()].
#'
#' @param GAlist GACT database list
#' @param overwrite Logical; if TRUE, recompute if existing .rds files exist
#' @return Updated GAlist with marker-level regulatory sets
#' @export
createRegulatoryMarkerSets <- function(GAlist, overwrite = FALSE) {

 gsets_dir <- normalizePath(GAlist$dirs["gsets"], winslash = "/", mustWork = FALSE)
 message("🧩 Creating SNP-level marker sets for regulatory elements")

 ensr2ensg_file <- file.path(gsets_dir, "ensr2ensg.rds")
 ensg2rsids_file <- file.path(gsets_dir, "ensg2rsids.rds")
 out_file <- file.path(gsets_dir, "ensr2rsids.rds")

 if (!file.exists(ensr2ensg_file) || !file.exists(ensg2rsids_file)) {
  stop("Missing ENSR or ENSG mapping files — run createRegulatorySets() first.")
 }

 if (!overwrite && file.exists(out_file)) {
  message("✅ ensr2rsids.rds already exists — skipping.")
  return(GAlist)
 }

 ensr2ensg <- readRDS(ensr2ensg_file)
 ensg2rsids <- readRDS(ensg2rsids_file)

 message("⚙️ Mapping regulatory elements to SNPs")
 ensr2rsids <- lapply(ensr2ensg, function(x) unique(unlist(ensg2rsids[x])))
 ensr2rsids <- ensr2rsids[!sapply(ensr2rsids, is.null)]
 saveRDS(ensr2rsids, out_file)

 message("✅ Saved regulatory SNP-level sets: ", basename(out_file))
 GAlist$sets_metadata$modules <-
  unique(c(GAlist$sets_metadata$modules, "regulatory_marker"))
 GAlist
}


#' @export
#'
createSetsDB0 <- function(GAlist = NULL, what="ensembl",
                         upstream=35, downstream=10,
                         min_combined_score=900, min_interactions=5) {

 # default sets (order of execution is important)

 GAlist$map <- vector(mode = "list", length = length(GAlist$mapfiles))
 for(i in 1:length(GAlist$mapfiles)) {
  #GAlist$map[[i]] <- readRDS(file.path(GAlist$dirs["gsets"],GAlist$mapfiles[i]))
  GAlist$map[[i]] <- readRDS(GAlist$mapfiles[i])
 }
 names(GAlist$map) <- names(GAlist$mapfiles)
 GAlist$map <- GAlist$map[c("eg2ensg", "ensg2eg", "ensg2sym", "ensp2ensg")]


 # Ensembl genes, proteins, transcripts
 #file <- file.path(GAlist$dirs["gsets"],"GRCh38.109.entrez.tsv.gz")
 file <- file.path(GAlist$dirs["ensembl"],"GRCh38.109.entrez.tsv.gz")
 ensembl <- fread(file, data.table=FALSE)
 ensembl <- ensembl[!ensembl$protein_stable_id=="-",]

 GAlist$map$eg2ensg <- split( ensembl$gene_stable_id, f=as.factor(ensembl$xref) )
 GAlist$map$eg2ensp <- split( ensembl$protein_stable_id, f=as.factor(ensembl$xref) )
 GAlist$map$eg2enst <- split( ensembl$transcript_stable_id, f=as.factor(ensembl$xref) )

 GAlist$map$ensp2eg <- split( ensembl$xref, f=as.factor(ensembl$protein_stable_id) )
 GAlist$map$ensp2ensg <- split( ensembl$gene_stable_id, f=as.factor(ensembl$protein_stable_id) )
 GAlist$map$ensp2enst <- split( ensembl$transcript_stable_id, f=as.factor(ensembl$protein_stable_id) )

 GAlist$map$ensg2eg <- split( ensembl$xref, f=as.factor(ensembl$gene_stable_id) )
 GAlist$map$ensg2ensp <- split( ensembl$protein_stable_id, f=as.factor(ensembl$gene_stable_id) )
 GAlist$map$ensg2enst <- split( ensembl$transcript_stable_id, f=as.factor(ensembl$gene_stable_id) )

 saveRDS(GAlist$map$eg2ensg, file = file.path(GAlist$dirs["gsets"], "eg2ensg.rds"))
 saveRDS(GAlist$map$eg2ensp, file = file.path(GAlist$dirs["gsets"], "eg2ensp.rds"))
 saveRDS(GAlist$map$eg2enst, file = file.path(GAlist$dirs["gsets"], "eg2enst.rds"))

 saveRDS(GAlist$map$ensg2eg, file = file.path(GAlist$dirs["gsets"], "ensg2eg.rds"))
 saveRDS(GAlist$map$ensg2ensp, file = file.path(GAlist$dirs["gsets"], "ensg2ensp.rds"))
 saveRDS(GAlist$map$ensg2enst, file = file.path(GAlist$dirs["gsets"], "ensg2enst.rds"))

 saveRDS(GAlist$map$ensp2eg, file = file.path(GAlist$dirs["gsets"], "ensp2eg.rds"))
 saveRDS(GAlist$map$ensp2ensg, file = file.path(GAlist$dirs["gsets"], "ensp2ensg.rds"))
 saveRDS(GAlist$map$ensp2enst, file = file.path(GAlist$dirs["gsets"], "ensp2enst.rds"))

 # Gene marker sets
 # Specify parameters
 #upstream <- 35
 #downstream <- 10

 markers <- fread(file.path(GAlist$dirs["marker"],"markers.txt.gz"),
                  data.table=FALSE)


 upstream <- upstream*1000
 downstream <- downstream*1000

 #filesGRC <-  c(file.path(GAlist$dirs["gsets"],"Homo_sapiens.GRCh37.87.gtf.gz"),
 #               file.path(GAlist$dirs["gsets"],"Homo_sapiens.GRCh38.109.gtf.gz"))
 filesGRC <-  c(file.path(GAlist$dirs["ensembl"],"Homo_sapiens.GRCh37.87.gtf.gz"),
                file.path(GAlist$dirs["ensembl"],"Homo_sapiens.GRCh38.109.gtf.gz"))
 setsfilesGRC <-  c(file.path(GAlist$dirs["gsets"],"GRCh37.ensg2rsids.rds"),
                    file.path(GAlist$dirs["gsets"],"GRCh38.ensg2rsids.rds"))

 for (j in seq_along(filesGRC)) {
  df <- fread(filesGRC[j], data.table=FALSE,skip=1, header=FALSE)
  colnames(df) <- c("chr","source","type","start","end","score","strand","phase","attributes")
  df <- df[df$type=="gene" & df$source=="ensembl_havana",]
  att <- strsplit(df$attributes, ";")
  att <- lapply(att, function(x){gsub("\"","",x)})
  gene_id <- sapply(att, function(x){ x[grep("gene_id",x)]})
  df$gene_id <- gsub("gene_id ","",gene_id)
  df <- df[,c("gene_id","chr","source", "type", "start", "end","strand")]
  df <- df[!df$chr=="X",]
  df <- df[!df$chr=="Y",]
  df <- df[as.character(df$chr)%in%as.character(1:22),]
  df$chr <- as.integer(df$chr)

  ensg2rsids <- vector("list", nrow(df))
  ensg2cpra <- vector("list", nrow(df))

  start <- df$start-upstream
  start[start<1] <- 1
  end <- df$end+downstream
  maxpos <- max(markers$pos,end)
  pos <- 1:maxpos
  ensg2rsids <- vector("list", nrow(df))
  for (chr in 1:22) {
   message(paste("Processing chr:",chr))
   rsids <- rep(NA, maxpos)
   rsids[as.integer(markers[markers$chr==chr,"pos"])] <- markers[markers$chr==chr,"rsids"]
   for (i in 1:nrow(df)) {
    if(df$chr[i]==chr) {
     grsids <- rsids[start[i]:end[i]]
     ensg2rsids[[i]] <- grsids[!is.na(grsids)]
    }
   }
  }
  names(ensg2rsids) <- df$gene_id
  empty <- sapply(ensg2rsids, function(x){ identical(x, character(0))})
  ensg2rsids <- ensg2rsids[!empty]
  saveRDS(ensg2rsids, file = setsfilesGRC[j])
 }

 ensg2rsids <- readRDS(file = file.path(GAlist$dirs["gsets"], "ensg2rsids.rds"))
 sets <- mapSetsDB(sets=ensg2rsids, featureID=markers$rsids, index=TRUE)
 chr <- sapply(sets, function(x) {unique(markers$chr[x])})
 start <- sapply(sets, function(x) {min(markers$pos[x])})
 stop <- sapply(sets, function(x) {max(markers$pos[x])})
 df <- data.frame(EnsemblID=names(sets),chr=chr,start=start,stop=stop)
 df$length <- df$stop-df$start
 saveRDS(df, file = file.path(GAlist$dirs["gsets"], "genesplus_annotation.rds"))

 ensg2rsids <- readRDS(file = file.path(GAlist$dirs["gsets"], "GRCh37.ensg2rsids.rds"))
 sets <- mapSetsDB(sets=ensg2rsids, featureID=markers$rsids, index=TRUE)
 chr <- sapply(sets, function(x) {unique(markers$chr[x])})
 start <- sapply(sets, function(x) {min(markers$pos[x])})
 stop <- sapply(sets, function(x) {max(markers$pos[x])})
 df <- data.frame(EnsemblID=names(sets),chr=chr,start=start,stop=stop)
 df$length <- df$stop-df$start
 saveRDS(df, file = file.path(GAlist$dirs["gsets"], "GRCh37_genesplus_annotation.rds"))

 ensg2rsids <- readRDS(file = file.path(GAlist$dirs["gsets"], "GRCh38.ensg2rsids.rds"))
 sets <- mapSetsDB(sets=ensg2rsids, featureID=markers$rsids, index=TRUE)
 chr <- sapply(sets, function(x) {unique(markers$chr[x])})
 start <- sapply(sets, function(x) {min(markers$pos[x])})
 stop <- sapply(sets, function(x) {max(markers$pos[x])})
 df <- data.frame(EnsemblID=names(sets),chr=chr,start=start,stop=stop)
 df$length <- df$stop-df$start
 saveRDS(df, file = file.path(GAlist$dirs["gsets"], "GRCh38_genesplus_annotation.rds"))


 # GWAS catalog
 #dbdir <- file.path(GAlist$dbdir, "gwas")
 dbdir <- file.path(GAlist$dbdir, "gwascatalog")
 gwasfile <- file.path(dbdir, "gwas-catalog-associations_ontology-annotated.tsv")
 gwas <- fread(gwasfile, data.table=FALSE, quote="")
 sets <- split(gwas$SNPS,f=gwas$MAPPED_TRAIT)
 saveRDS(sets, file = file.path(GAlist$dirs["gsets"], "gwas2rsids.rds"))


 # Pubmed to genes
 file <- file.path(GAlist$dirs["gsets"],"gene2pubmed.gz")
 pubmed <- fread(file, data.table = FALSE)
 pubmed <- pubmed[pubmed[,1]%in%9606,-1]
 eg2pmid <- split(pubmed$PubMed_ID,f=as.factor(pubmed$GeneID))
 pmid2eg <- split(pubmed$GeneID,f=as.factor(pubmed$PubMed_ID))
 saveRDS(eg2pmid,file=file.path(GAlist$dirs["gsets"],"eg2pmid.rds"))
 saveRDS(pmid2eg,file=file.path(GAlist$dirs["gsets"],"pmid2eg.rds"))

 #String
 file <- file.path(GAlist$dirs["string"],"9606.protein.links.v12.0.txt.gz")
 #file <- file.path(GAlist$dirs["gsets"],"9606.protein.links.v12.0.txt.gz")
 string <- fread(file, data.table=FALSE)
 string$protein1 <- gsub("9606.","",string$protein1)
 string$protein2 <- gsub("9606.","",string$protein2)
 string  <- string[string$combined_score>=min_combined_score,]
 string <- string[string$protein2%in%names(GAlist$map$ensp2ensg),]
 string <- split( string$protein2,f=as.factor(string$protein1))
 sets <- string[sapply(string ,length)>=min_interactions]
 sets <- lapply(sets,function(x){na.omit(unlist(GAlist$map$ensp2ensg[x]))})
 saveRDS(sets,file=file.path(GAlist$dirs["gsets"],"string2ensg.rds"))

 #Stitch
 file <- file.path(GAlist$dirs["stitch"],"9606.protein_chemical.links.v5.0.tsv.gz")
 #file <- file.path(GAlist$dirs["gsets"],"9606.protein_chemical.links.v5.0.tsv.gz")
 stitch <- fread(file, data.table=FALSE)
 stitch$protein <- gsub("9606.","",stitch$protein)
 stitch <- stitch[stitch$protein%in%names(GAlist$map$ensp2ensg),]
 stitch  <- stitch[stitch$combined_score>=min_combined_score,]
 stitch <- split( stitch$protein,f=as.factor(stitch$chemical))
 sets  <- stitch[sapply(stitch ,length)>=min_interactions]
 sets <- lapply(sets,function(x){na.omit(unlist(GAlist$map$ensp2ensg[x]))})
 saveRDS(sets,file=file.path(GAlist$dirs["gsets"],"stitch2ensg.rds"))

 # Reactome
 file <- file.path(GAlist$dirs["reactome"],"Ensembl2Reactome.txt")
 #file <- file.path(GAlist$dirs["gsets"],"Ensembl2Reactome.txt")
 if(file.exists(file)) {
  reactome <- fread(file, data.table=FALSE, header=FALSE)
  isHSA <- grep("R-HSA",reactome[,2])
  reactome <- reactome[isHSA,]
  reac2ensg <- split(reactome[,1],f=reactome[,2])
  ensg2reac <- split(reactome[,2],f=reactome[,1])
  GAlist$map$reac2ensg <- reac2ensg
  GAlist$map$ensg2reac <- ensg2reac
  saveRDS(reac2ensg,file=file.path(GAlist$dirs["gsets"],"reactome2ensg.rds"))
  saveRDS(ensg2reac,file=file.path(GAlist$dirs["gsets"],"ensg2reactome.rds"))
 }
 file  <- file.path(GAlist$dirs["reactome"],"ReactomePathways.txt")
 #file  <- file.path(GAlist$dirs["gsets"],"ReactomePathways.txt")
 if(file.exists(file)) {
  pathway <- fread(file, data.table=FALSE, header=FALSE)
  isHSA <- grep("R-HSA",pathway[,1])
  pathway <- pathway[isHSA,]
  reac2names <- pathway[,2]
  names(reac2names) <- pathway[,1]
  GAlist$map$reac2names <- reac2names
  saveRDS(reac2names,file=file.path(GAlist$dirs["gsets"],"reactome2names.rds"))
 }

 # Regulatory elements
 if("regulatory"%in%what) {

  build <- "GRCh37"
  #if(build=="GRCh37") file <- file.path(GAlist$dirs["gsets"],"GRCh37.Regulatory_Build.regulatory_features.gff.gz")
  #if(build=="GRCh38") file <- file.path(GAlist$dirs["gsets"],"GRCh38.Regulatory_Build.regulatory_features.gff.gz")
  if(build=="GRCh37") file <- file.path(GAlist$dirs["ensembl"],"GRCh37.Regulatory_Build.regulatory_features.gff.gz")
  if(build=="GRCh38") file <- file.path(GAlist$dirs["ensembl"],"GRCh38.Regulatory_Build.regulatory_features.gff.gz")
  df <- fread(file, data.table=FALSE)
  colnames(df) <- c("chr","source","type","start","end","score","strand","phase","attributes")
  df <- df[df$chr%in%as.character(1:22),]
  df$chr <- as.integer(df$chr)
  df <- df[!is.na(df$chr),]
  att <- strsplit(df$attributes, ";")
  att <- lapply(att, function(x){gsub("\"","",x)})
  att <- sapply(att, function(x){ x[grep("ID=",x)]})
  att <- strsplit(att, ":")
  df$reg_id <- sapply(att, function(x){x[2]})
  rownames(df) <- df$reg_id
  saveRDS(df[,c("reg_id", "type", "chr", "start", "end")],
          file = file.path(GAlist$dirs["gsets"],"regulatory_annotation.rds"))
  #file = file.path(GAlist$dirs["gsets"], paste0(build,"_regulatory_annotation.rds")))

  reg2ensr <- split(df$reg_id, f=as.factor(df$type))
  saveRDS(reg2ensr, file = file.path(GAlist$dirs["gsets"], "reg2ensr.rds"))
  #saveRDS(reg2ensr, file = file.path(GAlist$dirs["gsets"], paste0(build,"_reg2ensr.rds")))

  #markers <- fread(GAlist$markerfiles, data.table=FALSE)
  # this is currently build GRC37
  markers <- fread(file.path(GAlist$dirs["marker"],"markers.txt.gz"),
                   data.table=FALSE)
  start <- df$start
  start[start<1] <- 1
  end <- df$end
  maxpos <- max(markers$pos,end)
  pos <- 1:maxpos
  ensr2rsids <- vector("list", nrow(df))
  for (chr in 1:22) {
   message(paste("Processing chr:",chr))
   rsids <- rep(NA, maxpos)
   rsids[as.integer(markers[markers$chr==chr,"pos"])] <- markers[markers$chr==chr,"rsids"]
   for (i in 1:nrow(df)) {
    if(df$chr[i]==chr) {
     grsids <- rsids[start[i]:end[i]]
     ensr2rsids[[i]] <- grsids[!is.na(grsids)]
    }
   }
  }
  names(ensr2rsids) <- df$reg_id
  empty <- sapply(ensr2rsids, function(x){ identical(x, character(0))})
  ensr2rsids <- ensr2rsids[!empty]
  saveRDS(ensr2rsids, file = file.path(GAlist$dirs["gsets"], "ensr2rsids.rds"))

  reg2rsids <- sapply(reg2ensr, function(x){unique(unlist(ensr2rsids[x]))})
  saveRDS(reg2rsids, file = file.path(GAlist$dirs["gsets"], "reg2rsids.rds"))

  ensg2rsids <- readRDS(file.path(GAlist$dirs["gsets"], "ensg2rsids.rds"))
  ensr2rsids <- readRDS(file.path(GAlist$dirs["gsets"], "ensr2rsids.rds"))

  # Create a data frame from ensr2rsids
  nsets <- sapply(ensr2rsids, length)
  ensr <- rep(names(ensr2rsids), times = nsets)
  rsids <- unlist(ensr2rsids, use.names = FALSE)
  ensr_df <- data.frame(rsids, ensr)

  # Create a data frame from ensg2rsids
  nsets <- sapply(ensg2rsids, length)
  ensg <- rep(names(ensg2rsids), times = nsets)
  rsids <- unlist(ensg2rsids, use.names = FALSE)
  ensg_df <- data.frame(rsids, ensg)

  # Merge the data frames on rsids
  merged_df <- merge(ensr_df, ensg_df, by = "rsids")

  # Create the final list structure for ensr2ensg
  ensr2ensg <- split(merged_df$ensr, merged_df$ensg)
  ensr2ensg <- lapply(ensr2ensg, unique)
  empty <- sapply(ensr2ensg, function(x){ identical(x, character(0))})
  ensr2ensg <- ensr2ensg[!empty]
  saveRDS(ensr2ensg, file = file.path(GAlist$dirs["gsets"], "ensr2ensg.rds"))

  # Create the final list structure for ensg2ensr
  ensg2ensr <- split(merged_df$ensg, merged_df$ensr)
  ensg2ensr <- lapply(ensg2ensr, unique)
  empty <- sapply(ensg2ensr, function(x){ identical(x, character(0))})
  ensg2ensr <- ensg2ensr[!empty]
  saveRDS(ensg2ensr, file = file.path(GAlist$dirs["gsets"], "ensg2ensr.rds"))
 }

 # Drug databases
 drugdb <- fread(file.path(GAlist$dirs["drugdb"], "interactions.tsv"),
                 quote = "\"", data.table = FALSE)
 #hgnc <- fread("https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt", data.table=FALSE)
 hgnc <- fread("https://ftp.ebi.ac.uk/pub/databases/genenames/out_of_date_hgnc/tsv/hgnc_complete_set.txt", data.table=FALSE)
 hgnc2ensg <- hgnc$ensembl_gene_id
 names(hgnc2ensg) <- tolower(hgnc$hgnc_id)
 drug2hgnc <- split( drugdb$gene_concept_id, f=as.factor(drugdb$drug_name) )
 drug2ensg <- lapply(drug2hgnc,function(x){unique(na.omit(hgnc2ensg[x]))})
 str(drug2ensg)
 length(drug2ensg)
 drug2ensg <- drug2ensg[sapply(drug2ensg,length)>0]
 length(drug2ensg)
 saveRDS(drug2ensg,file=file.path(GAlist$dirs["gsets"],"drug2ensg.rds"))

 # Add atc codes
 #GAlist <- downloadDB(GAlist=GAlist, what="atc")
 # Add drug target data frame with ATC information to GAlist
 drug2ensg <- readRDS(file=file.path(GAlist$dirs["gsets"],"drug2ensg.rds"))
 #drug2ensg <- getSetsDB(GAlist = GAlist, feature = "DrugGenes")
 nreps <- sapply(drug2ensg,length)
 drugs <- rep(names(drug2ensg), times=nreps)
 ensg <- unlist(drug2ensg)
 ensg2drug <- split(drugs, f=as.factor(ensg))
 df <- data.frame(Drug=drugs, Target=ensg)
 df$ATC <- rep("Unknown",nrow(df))
 atc <- readRDS(file.path(GAlist$dirs["drugdb"],"atc.rds"))
 #has_atc <- match(tolower(df$Drug),tolower(GAlist$atc$name))
 #df$ATC[!is.na(has_atc)] <- as.character(GAlist$atc$code[has_atc[!is.na(has_atc)]])
 has_atc <- match(tolower(df$Drug),tolower(atc$name))
 df$ATC[!is.na(has_atc)] <- as.character(atc$code[has_atc[!is.na(has_atc)]])
 targets <- df
 saveRDS(targets, file=file.path(GAlist$dirs["gsets"],"targets.rds"))

 #GAlist$targets <- df
 #target <- GAlist$targets
 targets <- readRDS(file=file.path(GAlist$dirs["gsets"],"targets.rds"))
 targets <- targets[!duplicated(targets$Drug),]
 drug2atc <- targets$ATC
 names(drug2atc) <- targets$Drug
 saveRDS(drug2atc, file=file.path(GAlist$dirs["gsets"],"drug2atc.rds"))
 #GAlist$drug2atc <- drug2atc

 saveRDS(GAlist$map, file=file.path(GAlist$dirs["gsets"],"map.rds"))

 if("diseases"%in%what) {
  ensp <- names(GAlist$map$ensp2ensg)
  filenames <- c("human_disease_textmining_full.tsv.gz",
                 "human_disease_textmining_filtered.tsv.gz",
                 "human_disease_knowledge_full.tsv.gz",
                 "human_disease_knowledge_filtered.tsv.gz",
                 "human_disease_experiments_full.tsv.gz",
                 "human_disease_experiments_filtered.tsv.gz",
                 "human_disease_integrated_full.tsv.gz")
  for (file in filenames) {
   df <- fread(file.path(GAlist$dirs["gsets"],file), data.table=FALSE)
   df <- df[df[,1]%in%ensp,]
   disease2ensp <- split(df[,1],f=as.factor(as.character(df[,4])))
   sets <- mapSetsDB(disease2ensp,featureID=ensp,index=TRUE)
   disease2ensg <- lapply(sets,function(x) {unique(unlist(GAlist$map$ensp2ensg[x]))})
   ensg <- unlist(disease2ensg)
   disease <- rep(names(disease2ensg),times=sapply(disease2ensg,length))
   ensg2disease <- split(disease,f=as.factor(ensg))
   saveRDS(disease2ensp, file = file.path(GAlist$dirs["gsets"], paste0("disease2ensp_",gsub(".tsv.gz",".rds",file))))
   saveRDS(disease2ensg, file = file.path(GAlist$dirs["gsets"], paste0("disease2ensg_",gsub(".tsv.gz",".rds",file))))
   saveRDS(ensg2disease, file = file.path(GAlist$dirs["gsets"], paste0("ensg2disease_",gsub(".tsv.gz",".rds",file))))
  }
  return(GAlist)
 }

 if("atc"%in%what) {
  targets <- readRDS(file=file.path(GAlist$dirs["gsets"],"targets.rds"))
  targets <- targets[!duplicated(targets$Drug),]
  drug2atc <- targets$ATC
  names(drug2atc) <- targets$Drug
  #GAlist$drug2atc <- drug2atc
  drug2atc <- drug2atc[!drug2atc=="Unknown"]
  atcL1 <- substr(drug2atc, 1,1)
  atcL2 <- substr(drug2atc, 1,3)
  atcL3 <- substr(drug2atc, 1,4)
  atcL4 <- substr(drug2atc, 1,5)
  atcSets1 <- split(names(atcL1),f=atcL1)
  atcSets2 <- split(names(atcL2),f=atcL2)
  atcSets3 <- split(names(atcL3),f=atcL3)
  atcSets4 <- split(names(atcL4),f=atcL4)
  saveRDS(atcSets1, file = file.path(GAlist$dirs["gsets"], "atcSets1.rds"))
  saveRDS(atcSets2, file = file.path(GAlist$dirs["gsets"], "atcSets2.rds"))
  saveRDS(atcSets3, file = file.path(GAlist$dirs["gsets"], "atcSets3.rds"))
  saveRDS(atcSets4, file = file.path(GAlist$dirs["gsets"], "atcSets4.rds"))
  return(GAlist)
 }


 return(GAlist)
}


#' @export
#'
createMarkerSetsDB0 <- function(GAlist = NULL, what=NULL,
                               upstream=35, downstream=10,
                               min_combined_score=900, min_interactions=5) {

 ensg2rsids <- readRDS(file=file.path(GAlist$dirs["gsets"],"ensg2rsids.rds"))

 if("GO"%in%what) {
  fset <- readRDS(file=file.path(GAlist$dirs["gsets"],"go.rds"))
  fset <- mapSetsDB(fset,featureID=names(ensg2rsids))
  sets <- lapply(fset, function(x){unique(unlist(ensg2rsids[x]))})
  sets <- sets[!sapply(sets, is.null)]
  saveRDS(sets, file = file.path(GAlist$dirs["gsets"], "go2rsids.rds"))
 }

 if("reactome"%in%what) {
  fset <- readRDS(file=file.path(GAlist$dirs["gsets"],"reactome2ensg.rds"))
  fset <- mapSetsDB(fset,featureID=names(ensg2rsids))
  sets <- lapply(fset, function(x){unique(unlist(ensg2rsids[x]))})
  sets <- sets[!sapply(sets, is.null)]
  saveRDS(sets, file = file.path(GAlist$dirs["gsets"], "reactome2rsids.rds"))
 }

 if("string"%in%what) {
  fset <- readRDS(file=file.path(GAlist$dirs["gsets"],"string2ensg.rds"))
  fset <- mapSetsDB(fset,featureID=names(ensg2rsids))
  sets <- lapply(fset, function(x){unique(unlist(ensg2rsids[x]))})
  sets <- sets[!sapply(sets, is.null)]
  saveRDS(sets, file = file.path(GAlist$dirs["gsets"], "string2rsids.rds"))
 }

 if("stitch"%in%what) {
  fset <- readRDS(file=file.path(GAlist$dirs["gsets"],"stitch2ensg.rds"))
  fset <- mapSetsDB(fset,featureID=names(ensg2rsids))
  sets <- lapply(fset, function(x){unique(unlist(ensg2rsids[x]))})
  sets <- sets[!sapply(sets, is.null)]
  saveRDS(sets, file = file.path(GAlist$dirs["gsets"], "stitch2rsids.rds"))
 }

 if("drug"%in%what) {

  drug2ensg <- readRDS(file=file.path(GAlist$dirs["gsets"],"drug2ensg.rds"))
  string2ensg <- readRDS(file=file.path(GAlist$dirs["gsets"],"string2ensg.rds"))
  drug2ensp <- lapply(drug2ensg,function(x){na.omit(unlist(GAlist$map$ensg2ensp[x]))})
  drug2string2ensg <- lapply(drug2ensp,function(x){na.omit(unlist(string2ensg[x]))})
  drug2string2ensg <- lapply(drug2string2ensg, function(x){unique(x)})
  for(i in 1:length(drug2string2ensg)) {
   drug2string2ensg[[i]] <- unique(c(drug2ensg[[i]], drug2string2ensg[[i]]))
  }
  saveRDS(drug2string2ensg,file=file.path(GAlist$dirs["gsets"],"drug2string2ensg.rds"))

  ensg2rsids <- readRDS(file.path(GAlist$dirs["gsets"], "ensg2rsids.rds"))

  drug2rsids <- lapply(drug2ensg,function(x){unique(unlist(ensg2rsids[x]))})
  drug2rsids <- drug2rsids[!sapply(drug2rsids,is.null)]
  saveRDS(drug2rsids,file=file.path(GAlist$dirs["gsets"],"drug2rsids.rds"))

  drug2string2rsids <- lapply(drug2string2ensg,function(x){unique(unlist(ensg2rsids[x]))})
  drug2string2rsids <- drug2string2rsids[!sapply(drug2string2rsids,is.null)]
  saveRDS(drug2string2rsids,file=file.path(GAlist$dirs["gsets"],"drug2string2rsids.rds"))
 }

 return(GAlist)
}

