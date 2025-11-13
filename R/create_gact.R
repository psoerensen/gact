#######################################################################################
# Function used for creating the GACT database
#######################################################################################
#'
#' Create the database Genomic Association of Complex Traits (GACT)
#'
#' @description
#' The `gact` function is designed for establishing and populating a comprehensive
#' database focused on genomic associations with complex traits. This function serves
#' two primary roles: infrastructure creation and data acquisition. It facilitates the
#' assembly of a structured repository encompassing single marker associations,
#' rigorously curated to ensure high-quality data. Beyond individual genetic markers,
#' the function integrates a broad spectrum of genomic entities, encompassing genes,
#' proteins, and an array of complexes (chemical, protein) as well as various
#' biological pathways. This integration aims to provide a holistic view of genomic
#' associations and their multifaceted relationships in the context of complex traits.
#'
#' @param GAlist A list object providing information and infrastructure of the gact database.
#' @param version The version of the gact database to use.
#' @param task The task to perform, either "download" or "process". By default "download".
#' @param dbdir The directory where the database should be stored, by default the current working directory.
#'
#' @return A list providing information and infrastructure of the gact database.
#'
#' @import data.table
#'
#' @examples
#' \dontrun{
#' GAlist <- gact()
#' }
#'
#' @importFrom R.utils gzip
#' @export
#'
gact <- function(GAlist=NULL, version=NULL, task="download",
                 dbdir=NULL, timeout=3000) {

 if (is.null(dbdir)) dbdir <- getwd()
 if (is.null(version))
  stop("Please provide a database version (e.g. version = 'hsa.0.0.1').\n",
       "Tip: Run list_gact_versions() to see available versions.")

 # ------------------------------------------------------------------
 # ✅ Check version exists before creating anything
 # ------------------------------------------------------------------
 if (!exists("gact_versions", inherits = TRUE)) {
  stop("No version metadata ('gact_versions') found in your environment.\n",
       "This usually means the gact package is not fully loaded.")
 }

 available_versions <- names(gact_versions)
 if (!(version %in% available_versions)) {
  stop(
   paste0(
    "Version '", version, "' is not defined in gact_versions.\n",
    "Available versions are: ", paste(available_versions, collapse = ", "), "\n",
    "Tip: Run list_gact_versions() to view supported database versions."
   )
  )
 }

 if(task=="download") {
  GAlist <- createDB(version=version, dbdir=dbdir)

  options(download.file.method="libcurl", url.method="libcurl", timeout=timeout)

  # Step 2: Download data from database:

  # default version of gact
  GAlist <- downloadDB(GAlist=GAlist, what="annotation")
  GAlist <- downloadDB(GAlist=GAlist, what="marker")
  GAlist <- downloadDB(GAlist=GAlist, what="gsea")
  GAlist <- downloadDB(GAlist=GAlist, what="gstat")
  GAlist <- downloadDB(GAlist=GAlist, what="gbayes")
  GAlist <- downloadDB(GAlist=GAlist, what="atc")
  GAlist <- downloadDB(GAlist=GAlist, what="dgi")

  # these are ressources used for building gene and marker sets
  # (not necessary to download unless custom sets is to be build)
  #GAlist <- downloadDB(GAlist=GAlist, what="ensembl")
  #GAlist <- downloadDB(GAlist=GAlist, what="reactome")
  #GAlist <- downloadDB(GAlist=GAlist, what="string")
  #GAlist <- downloadDB(GAlist=GAlist, what="stitch")

  # these can be downloaded by the user if they want to do some specialized analyses
  #GAlist <- downloadDB(GAlist=GAlist, what="gtex")
  #GAlist <- downloadDB(GAlist=GAlist, what="gwascatalog")
  #GAlist <- downloadDB(GAlist=GAlist, what="1000G")
  #GAlist <- downloadDB(GAlist=GAlist, what="diseases")

  # these may be ignored but potentially used at a later stage
  #GAlist <- downloadDB(GAlist=GAlist, what="gsets")
  #GAlist <- downloadDB(GAlist=GAlist, what="pubmed")
  #GAlist <- downloadDB(GAlist=GAlist, what="vep")
  #GAlist <- downloadDB(GAlist=GAlist, what="tiga")
  #GAlist <- downloadDB(GAlist=GAlist, what="pubchem")
  #GAlist <- downloadDB(GAlist=GAlist, what="pharmgkb")
  #GAlist <- downloadDB(GAlist=GAlist, what="opentargets")
  #GAlist <- downloadDB(GAlist=GAlist, what="alphamissense")

  #message("Creating full marker sets - this may take some time")
  #GAlist <- createSetsDB(GAlist=GAlist)
  #GAlist <- createSetsDB(GAlist=GAlist, what="diseases")
  #GAlist <- createMarkerSetsDB(GAlist=GAlist, what="GO")
  #GAlist <- createMarkerSetsDB(GAlist=GAlist, what="string")
  #GAlist <- createMarkerSetsDB(GAlist=GAlist, what="stitch")
  #GAlist <- createMarkerSetsDB(GAlist=GAlist, what="drug")
  #summaryDB(GAlist=GAlist)

 }
 # Step 3: Create marker sets from database:
 if(task=="createSets") {
  message("Creating full marker sets - this may take some time")
  GAlist <- createSetsDB(GAlist=GAlist)
  GAlist <- createMarkerSetsDB(GAlist=GAlist, what="reactome")
  #GAlist <- createMarkerSetsDB(GAlist=GAlist, what="GO")
  GAlist <- createMarkerSetsDB(GAlist=GAlist, what="string")
  #GAlist <- createMarkerSetsDB(GAlist=GAlist, what="stitch")
 }
 return(GAlist)
}

#' Create GACT database directory structure
#'
#' Initializes the directory tree used by downloadDB() and other GACT functions.
#'
#' @param version Character; name or version of the database (e.g., "v0.1").
#' @param dbdir Character; parent directory path in which the version folder will be created.
#'
#' @return A list (GAlist) containing the database version, base directory,
#'         directory paths (GAlist$dirs), and feature names.
#' @export
createDB <- function(version = NULL, dbdir = NULL) {
 if (is.null(version)) stop("Please include a database name using the version argument")

 # Normalize and check directory path
 dbdir <- normalizePath(file.path(dbdir, version), winslash = "/", mustWork = FALSE)
 if (dir.exists(dbdir)) stop(paste("Directory:", dbdir, "already exists"))

 # Create database directory structure
 dir.create(dbdir, recursive = TRUE)
 dirnames <- c(
  "glist", "gstat", "gsets", "gsea", "gbayes", "gtex", "gwas", "ldsc",
  "marker", "drugdb", "download", "script", "vep",
  # optional resources
  "ensembl", "string", "stitch", "reactome", "dgi", "gwascatalog", "diseases"
 )
 dirs <- setNames(file.path(dbdir, dirnames), dirnames)
 lapply(dirs, dir.create)

 # Initialize GAlist
 features <- c("Markers", "Genes", "Proteins", "GO", "Pathways",
               "ProteinComplexes", "ChemicalComplexes", "GTEx", "GWAS")
 GAlist <- list(version = version, dbdir = dbdir, dirs = dirs, features = features)

 # Attach version metadata
 if (exists("gact_versions", inherits = TRUE) && version %in% names(gact_versions)) {
  GAlist$metadata <- gact_versions[[version]]
 } else {
  warning("Version metadata not found in gact_versions for version: ", version)
 }

 message("✅ Database directory created at: ", dbdir)
 return(GAlist)
}



#' Metadata describing available GACT database versions and sources
#'
#' This list defines all database versions available within GACT, including
#' species, description, creation date, and the corresponding data sources
#' (DOIs or URLs) for each component. These are used by [downloadDB()] to
#' resolve download targets.
#'
#' @keywords internal
gact_versions <- list(
 "hsa.0.0.1" = list(
  species = "Homo sapiens",
  description = "Default GACT human database release 0.0.1",
  created = "2025-11-06",
  sources = list(
   # Core Zenodo resources
   annotation = "10.5281/zenodo.14234857",
   marker     = "10.5281/zenodo.10462385",
   gsea       = "10.5281/zenodo.10462484",
   gstat      = "10.5281/zenodo.10462495",
   gbayes     = "10.5281/zenodo.10462421",

   # Drug classification
   atc = "https://www.dropbox.com/s/n5ehglmhhs0kcue/WHO%20ATC-DDD%202023-03-28.csv?dl=1",

   # Optional reference data
   ensembl = paste(
    "https://ftp.ensembl.org/pub/release-109/tsv/homo_sapiens/Homo_sapiens.GRCh38.109.entrez.tsv.gz",
    "https://ftp.ensembl.org/pub/release-109/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz",
    "https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz",
    "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz",
    "https://ftp.ensembl.org/pub/grch37/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz",
    "https://ftp.ensembl.org/pub/grch37/release-113/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff.gz",
    sep = ";"
   ),

   string = "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz",

   stitch = "http://stitch.embl.de/download/protein_chemical.links.v5.0/9606.protein_chemical.links.v5.0.tsv.gz",

   reactome = paste(
    "https://reactome.org/download/current/ReactomePathways.txt",
    "https://reactome.org/download/current/Ensembl2Reactome.txt",
    sep = ";"
   ),

   dgi = paste(
    "https://www.dgidb.org/data/2023-Dec/interactions.tsv",
    "https://www.dgidb.org/data/2023-Dec/genes.tsv",
    "https://www.dgidb.org/data/2023-Dec/drugs.tsv",
    "https://www.dgidb.org/data/2023-Dec/categories.tsv",
    sep = ";"
   ),

   gtex = paste(
    "https://storage.googleapis.com/adult-gtex/bulk-qtl/v7/single-tissue-cis-qtl/GTEx_Analysis_v7_eQTL.tar.gz",
    "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar",
    sep = ";"
   ),

   gwascatalog = paste(
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-studies.tsv",
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip",
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-ancestry.tsv",
    sep = ";"
   ),

   diseases = paste(
    "https://download.jensenlab.org/human_disease_textmining_full.tsv",
    "https://download.jensenlab.org/human_disease_textmining_filtered.tsv",
    "https://download.jensenlab.org/human_disease_knowledge_full.tsv",
    "https://download.jensenlab.org/human_disease_knowledge_filtered.tsv",
    "https://download.jensenlab.org/human_disease_experiments_full.tsv",
    "https://download.jensenlab.org/human_disease_experiments_filtered.tsv",
    "https://download.jensenlab.org/human_disease_integrated_full.tsv",
    sep = ";"
   ),

   # 1000 Genomes reference panels (Zenodo)
   `1000G` = paste(
    "10.5281/zenodo.10462403",  # 1000G_EUR_Phase3_plink.zip
    "10.5281/zenodo.14141438",  # g1000_eur.zip
    "10.5281/zenodo.14141523",  # g1000_eas.zip
    "10.5281/zenodo.14141597",  # g1000_sas.zip
    "10.5281/zenodo.14141774",  # g1000_afr.zip
    "10.5281/zenodo.14141885",  # g1000_amr.zip
    sep = ";"
   )


  )
 )
)

#' List available GACT database versions
#'
#' @description
#' Displays all database versions shipped with the GACT package, including
#' species, creation date, and a short description. These correspond to
#' versioned data bundles that can be installed using [gact()].
#'
#' @return A data frame with columns: `version`, `species`, `created`, and `description`.
#' @export
#'
#' @examples
#' list_gact_versions()
list_gact_versions <- function() {
 versions <- lapply(names(gact_versions), function(v) {
  list(
   version = v,
   species = gact_versions[[v]]$species,
   created = gact_versions[[v]]$created,
   description = gact_versions[[v]]$description
  )
 })

 df <- do.call(rbind, lapply(versions, as.data.frame, stringsAsFactors = FALSE))
 rownames(df) <- NULL

 cat("📦 Available GACT database versions:\n\n")
 print(df, row.names = FALSE)
 invisible(df)
}



#' Download database components for genomic associations (with logging)
#'
#' This function downloads the various components of a genomic associations
#' database and logs metadata for each operation via an external `log_download()`
#' helper function. It automatically retrieves DOIs or URLs from version metadata
#' defined in `gact_versions`.
#'
#' @param GAlist A list object containing the structure and paths needed for the database.
#' @param what A character string specifying which data type to download.
#' @param min_combined_score Minimum combined score for filtering STRING/STITCH data.
#' @param min_interactions Minimum number of interactions for filtering STRING/STITCH data.
#'
#' @return The updated GAlist object with file paths and GAlist$downloads updated.
#' @export
downloadDB <- function(GAlist = NULL, what = NULL,
                       min_combined_score = 900, min_interactions = 5) {

 if (is.null(what)) stop("Please specify what to download, e.g. what='gsets'")

 # Ensure timeout and libcurl settings
 options(download.file.method = "libcurl", url.method = "libcurl", timeout = 3000)

 # --------------------------------------------------------------------
 # Version-aware source lookup
 # --------------------------------------------------------------------
 source_map <- NULL
 if (!is.null(GAlist$metadata) && !is.null(GAlist$metadata$sources))
  source_map <- GAlist$metadata$sources

 get_source <- function(what) {
  if (!is.null(source_map) && what %in% names(source_map)) {
   return(source_map[[what]])
  } else {
   return(NULL)
  }
 }

 # Helper for safe downloading + logging
 #safe_download <- function(expr, what, source, urls, dest, success_msg, fail_msg) {
 # tryCatch({
 #  force(expr)
 #  GAlist <<- log_download(GAlist, what, source, paste(urls, collapse = "; "), dest, "success", success_msg)
 # }, error = function(e) {
 #  GAlist <<- log_download(GAlist, what, source, paste(urls, collapse = "; "), dest, "failed", paste(fail_msg, e$message))
 # })
 #}
 #safe_download <- function(expr, what, source, urls, dest, success_msg, fail_msg) {
 # dest <- normalizePath(dest, winslash = "/", mustWork = FALSE)  # ✅ normalize here
 # tryCatch({
 #  force(expr)
 #  GAlist <<- log_download(GAlist, what, source, paste(urls, collapse = "; "), dest, "success", success_msg)
 # }, error = function(e) {
 #  GAlist <<- log_download(GAlist, what, source, paste(urls, collapse = "; "), dest, "failed", paste(fail_msg, e$message))
 # })
 #}
 safe_download <- function(expr, what, source, urls, dest, success_msg, fail_msg) {
  dest <- normalizePath(dest, winslash = "/", mustWork = FALSE)
  tryCatch({
   expr_result <- eval(substitute(expr), envir = parent.frame())
   GAlist <- log_download(GAlist, what, source, paste(urls, collapse = "; "),
                          dest, "success", success_msg)
   return(GAlist)
  }, error = function(e) {
   GAlist <- log_download(GAlist, what, source, paste(urls, collapse = "; "),
                          dest, "failed", paste(fail_msg, e$message))
   return(GAlist)
  })
 }


 # --------------------------------------------------------------------
 # Annotation (Zenodo)
 # --------------------------------------------------------------------
 if (what == "annotation") {
  dest <- GAlist$dirs["gsets"]
  doi <- get_source("annotation")
  if (is.null(doi)) stop("No DOI or URL found for 'annotation' in version metadata.")
  message("Downloading annotation sets")
  GAlist <- safe_download({
   download_zenodo(doi = doi, path = dest, parallel = TRUE)
   zipfile <- file.path(dest, "hsa.0.0.1.zip")
   Sys.sleep(1)
   unzip(zipfile, exdir = dest)
   GAlist$map <- readRDS(file.path(dest, "map.rds"))
  }, what, "Zenodo", doi, dest, "Annotation sets downloaded", "Annotation sets failed: ")
 }

 # --------------------------------------------------------------------
 # Marker information (Zenodo)
 # --------------------------------------------------------------------
 if (what == "marker") {
  dest <- GAlist$dirs["marker"]
  doi <- get_source("marker")
  if (is.null(doi)) stop("No DOI or URL found for 'marker' in version metadata.")
  message("Downloading marker information")
  GAlist <- safe_download({
   download_zenodo(doi = doi, path = dest, parallel = TRUE)
   markers <- fread(file.path(dest, "markers.txt.gz"), data.table = FALSE)
   GAlist$rsids <- markers$rsids
   cpra <- paste(markers$chr, markers$pos, markers$ea, markers$nea, sep = "_")
   fwrite(data.frame(cpra = cpra), file = file.path(dest, "cpra.txt.gz"), col.names = FALSE)
  }, what, "Zenodo", doi, dest, "Marker data downloaded", "Marker data failed: ")
 }


 # --------------------------------------------------------------------
 # GSEA (Zenodo)
 # --------------------------------------------------------------------
 if (what == "gsea") {
  dest <- GAlist$dirs["gsea"]
  doi <- get_source("gsea")
  if (is.null(doi)) stop("No DOI or URL found for 'gsea' in version metadata.")
  message("Downloading GSEA summary statistics")
  GAlist <- safe_download({
   download_zenodo(doi = doi, path = dest, parallel = TRUE)
  }, what, "Zenodo", doi, dest, "GSEA statistics downloaded", "GSEA download failed: ")
 }

 # --------------------------------------------------------------------
 # GWAS summary statistics (Zenodo)
 # --------------------------------------------------------------------
 if (what == "gstat") {
  dest <- GAlist$dirs["gstat"]
  doi <- get_source("gstat")
  if (is.null(doi)) stop("No DOI or URL found for 'gstat' in version metadata.")
  message("Downloading GWAS summary statistics")
  GAlist <- safe_download({
   download_zenodo(doi = doi, path = dest, parallel = TRUE)
   meta <- file.path(dest, "GWAS_information.csv")
   GAlist$study <- as.list(read.csv2(meta))
   GAlist$studies <- as.data.frame(GAlist$study)
   rownames(GAlist$studies) <- GAlist$studies$id
   gz <- list.files(dest, pattern = ".gz", full.names = TRUE)
   GAlist$studyfiles <- gz
   names(GAlist$studyfiles) <- gsub(".txt.gz", "", basename(gz))
  }, what, "Zenodo", doi, dest, "GWAS summary statistics downloaded", "GWAS download failed: ")
 }

 # --------------------------------------------------------------------
 # Bayesian statistics (Zenodo)
 # --------------------------------------------------------------------
 if (what == "gbayes") {
  dest <- GAlist$dirs["gbayes"]
  doi <- get_source("gbayes")
  if (is.null(doi)) stop("No DOI or URL found for 'gbayes' in version metadata.")
  message("Downloading Bayesian fine-mapping results")
  GAlist <- safe_download({
   download_zenodo(doi = doi, path = dest, parallel = TRUE)
   files <- list.files(dest, pattern = "^fit_blr_pruned_GWAS[0-9]+\\.rds$", full.names = TRUE)
   old_names <- basename(files)
   new_names <- gsub("^fit_blr_pruned_GWAS", "BLR", old_names)
   gwas_id <- sub("^fit_blr_pruned_(.*?)\\.rds$", "\\1", old_names)

   ok <- file.rename(from = files, to = file.path(dest, new_names))
   if (any(!ok)) warning("Some files could not be renamed.")

   meta <- read.csv2(file.path(GAlist$dirs["gstat"], "GWAS_information.csv"), row.names = 1)
   meta <- meta[intersect(rownames(meta), gwas_id), ]
   rownames(meta) <- gsub("GWAS", "BLR", rownames(meta))
   meta$file <- gsub("GWAS", "BLR", meta$file)
   write.csv2(meta, file.path(GAlist$dirs["gbayes"], "BLR_information.csv"), quote = TRUE)
  }, what, "Zenodo", doi, dest, "Bayesian results downloaded", "Bayesian download failed: ")
 }

 # --------------------------------------------------------------------
 # 1000 Genomes reference data (Zenodo)
 # --------------------------------------------------------------------
 if (what == "1000G") {
  dest <- GAlist$dirs["marker"]
  doi <- get_source("1000G")
  if (is.null(doi))
   stop("No DOI or URL found for '1000G' in version metadata.")

  message("Downloading 1000 Genomes reference panels (≈ several GB total)")

  GAlist <- safe_download({
   urls <- unlist(strsplit(doi, ";"))

   # Iterate through each DOI and download sequentially
   for (u in urls) {
    message("Downloading from Zenodo DOI: ", u)
    download_zenodo(doi = u, path = dest, parallel = FALSE)
   }

   # Find all ZIP files and unzip them
   zip_files <- list.files(dest, pattern = "\\.zip$", full.names = TRUE)
   for (z in zip_files) {
    message("Unzipping: ", basename(z))
    unzip(z, exdir = dest)
   }

   message("✅ 1000 Genomes reference data successfully extracted")

  }, what, "Zenodo", doi, dest,
  "1000 Genomes reference data downloaded",
  "1000 Genomes download failed: ")
 }


 # --------------------------------------------------------------------
 # GTEx eQTL data (Google Cloud Storage)
 # --------------------------------------------------------------------
 if (what == "gtex") {
  dest <- GAlist$dirs["gtex"]
  urls <- get_source("gtex")
  if (is.null(urls))
   stop("No URL found for 'gtex' in version metadata.")

  message("🧬 Downloading GTEx eQTL datasets (≈ several GB total)")

  GAlist <- safe_download({
   urls <- unlist(strsplit(urls, ";"))
   dir.create(dest, recursive = TRUE, showWarnings = FALSE)

   # Iterate through each GTEx archive (v7, v8, etc.)
   for (u in urls) {
    fname <- basename(u)
    destfile <- file.path(dest, fname)

    message("⬇️  Downloading: ", fname)
    download.file(u, destfile, mode = "wb")

    # Extract if it's a tar or tar.gz file
    if (grepl("\\.tar(\\.gz)?$", fname, ignore.case = TRUE)) {
     message("📦 Extracting: ", fname)
     untar(destfile, exdir = dest)
     # Optional cleanup of large tar files
     file.remove(destfile)
    } else {
     message("ℹ️ Skipping extraction (not a tar file): ", fname)
    }
   }

   message("✅ GTEx eQTL data successfully extracted to: ", dest)

  }, what, "GTEx (Google Cloud Storage)", urls, dest,
  "GTEx eQTL data downloaded",
  "GTEx eQTL download failed: ")
 }


 # --------------------------------------------------------------------
 # Optional data sources (STRING, STITCH, Ensembl, Reactome, DGIdb, GTEx, GWAS, DISEASES)
 # --------------------------------------------------------------------
 optional_resources <- list(
  ensembl      = get_source("ensembl"),
  string       = get_source("string"),
  stitch       = get_source("stitch"),
  reactome     = get_source("reactome"),
  dgi          = get_source("dgi"),
  #gtex         = get_source("gtex"),
  gwascatalog  = get_source("gwascatalog"),
  diseases     = get_source("diseases")
 )

 # Handle direct FTP/HTTP downloads for optional resources
 for (res in names(optional_resources)) {
  if (what == res && !is.null(optional_resources[[res]])) {

   urls <- unlist(strsplit(optional_resources[[res]], ";"))

   # ✅ Explicitly store DGIdb and ATC in drugdb/
   if (res %in% c("dgi", "atc")) {
    dest <- GAlist$dirs["drugdb"]
   } else if (!is.null(GAlist$dirs[[res]])) {
    dest <- GAlist$dirs[[res]]
   } else {
    dest <- file.path(GAlist$dbdir, res)
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    message("Created missing directory: ", dest)
   }

   message("Downloading resource: ", res)

   #safe_download({
   # for (u in urls) {
   #  destfile <- file.path(dest, basename(u))
   #  download.file(u, destfile, mode = "wb")
   # }
   #}, res, res, urls, dest,
   #paste(res, "downloaded"),
   #paste(res, "failed: "))
  GAlist <- safe_download({
  for (u in urls) {
    destfile <- file.path(dest, basename(u))
    message("Downloading: ", basename(u))
    download.file(u, destfile, mode = "wb")
  }

  # 🧠 Post-processing for Ensembl files
  if (res == "ensembl") {
    message("🔧 Standardizing Ensembl filenames ...")

    files <- list.files(dest, pattern = "Homo_sapiens\\.|homo_sapiens\\.", full.names = TRUE)
    for (f in files) {
      base <- basename(f)
      new <- sub("^(H|h)omo_sapiens\\.", "", base)                 # remove prefix
      new <- sub("\\.\\d{8}(?=\\.gff\\.gz$)", "", new, perl = TRUE) # remove date suffix
      new <- sub("homo_sapiens\\.", "", new)                        # cleanup duplicate
      newpath <- file.path(dest, new)

      if (base != new) {
        message("Renaming: ", base, " → ", new)
        file.rename(f, newpath)
      }
    }
  }

}, res, res, urls, dest,
paste(res, "downloaded"),
paste(res, "failed: "))

  }
 }

 # --------------------------------------------------------------------
 # ATC classification (WHO)
 # --------------------------------------------------------------------
 if (what == "atc") {
  dest <- GAlist$dirs["drugdb"]
  url <- get_source("atc")
  if (is.null(url)) stop("No URL found for 'atc' in version metadata.")
  message("Downloading WHO ATC classification")
  GAlist <- safe_download({
   df <- fread(url, data.table = FALSE)
   atc <- list(code = df$atc_code, name = df$atc_name)
   names(atc$code) <- atc$name
   names(atc$name) <- atc$code
   saveRDS(atc, file.path(dest, "atc.rds"))
  }, what, "WHO ATC", url, dest, "ATC data downloaded", "ATC download failed: ")
 }

 return(GAlist)
}




#' Log downloaded files in the GACT database structure
#'
#' Records one row per successfully downloaded (or failed) file in
#' `GAlist$downloads`. Each entry includes file metadata such as URL,
#' destination path, file size, and timestamp.
#'
#' This function supports multiple URLs: if `urls` contains a semicolon-separated
#' string or a character vector of URLs, one entry is created for each file.
#' It appends new entries to any existing `GAlist$downloads` table.
#'
#' @param GAlist A GACT database list object, typically created by
#'   [createDB()] and updated by [downloadDB()].
#' @param what Character string indicating the dataset name
#'   (e.g. `"ensembl"`, `"string"`, `"gstat"`).
#' @param source Character string describing the data source
#'   (e.g. `"Zenodo"`, `"Ensembl FTP"`, `"DGIdb"`).
#' @param urls Character vector or semicolon-separated string of download URLs
#'   or DOIs corresponding to the files being logged.
#' @param dest_path Character string giving the directory where the files are
#'   downloaded.
#' @param status Character string describing the result of the download,
#'   typically `"success"` or `"failed"`.
#' @param message Optional character string with additional information or
#'   error message.
#'
#' @return The updated `GAlist` object with a data frame `GAlist$downloads`
#'   containing one row per file and the following columns:
#'   \itemize{
#'     \item `what` — dataset identifier
#'     \item `source` — data source
#'     \item `url` — download URL or DOI
#'     \item `dest_file` — destination file path
#'     \item `file_size_MB` — file size in megabytes (if found)
#'     \item `status` — `"success"` or `"failed"`
#'     \item `message` — optional descriptive message
#'     \item `date_downloaded` — timestamp of logging
#'   }
#'
#' @examples
#' \dontrun{
#' GAlist <- list(dbdir = tempdir())
#' GAlist <- log_download(
#'   GAlist,
#'   what = "ensembl",
#'   source = "Ensembl FTP",
#'   urls = c("https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz"),
#'   dest_path = file.path(GAlist$dbdir, "ensembl"),
#'   status = "success",
#'   message = "Downloaded GTF successfully"
#' )
#' }
#'
#' @export
log_download <- function(GAlist, what, source, urls, dest_path, status, message = "") {
 # Handle multiple URLs
 if (length(urls) == 1) urls <- strsplit(urls, ";")[[1]]
 urls <- trimws(urls)

 entries <- lapply(urls, function(u) {
  dest_file <- file.path(dest_path, basename(u))
  size <- if (file.exists(dest_file)) file.info(dest_file)$size / 1e6 else NA
  data.frame(
   what = what,
   source = source,
   url = u,
   dest_file = dest_file,
   file_size_MB = round(size, 3),
   status = status,
   message = message,
   date_downloaded = Sys.time(),
   stringsAsFactors = FALSE
  )
 })

 entry_df <- do.call(rbind, entries)

 if (!"downloads" %in% names(GAlist) || is.null(GAlist$downloads)) {
  GAlist$downloads <- entry_df
 } else {
  GAlist$downloads <- rbind(GAlist$downloads, entry_df)
 }

 return(GAlist)
}


#' Check which GACT database components are downloaded and complete
#'
#' Compares expected GACT directories with on-disk files and download logs.
#' Handles unpacked bundles like annotation stored in gsets/, and shared folders
#' like drugdb/ for dgi and atc datasets.
#'
#' @param GAlist A GACT database list created by [createDB()] and updated by [downloadDB()].
#' @param min_size_MB Minimum total size (in MB) to consider a dataset complete.
#' @param verbose Logical; if TRUE, prints a summary table.
#'
#' @return A data.frame summarizing download status per dataset.
#' @export
check_downloads <- function(GAlist, min_size_MB = 0.5, verbose = TRUE) {

 if (is.null(GAlist$dirs))
  stop("GAlist must contain a 'dirs' element from createDB().")

 # Expected datasets to check
 expected <- c(
  "annotation", "marker", "gsea", "gstat", "gbayes",
  "ensembl", "string", "stitch", "reactome",
  "dgi", "gtex", "gwascatalog", "atc"
 )

 # Map datasets to their actual storage folders
 folder_map <- c(
  annotation  = "gsets",
  dgi         = "drugdb",
  atc         = "drugdb",
  gwascatalog = "gwascatalog" # manually created in downloadDB
 )

 # Function to resolve correct folder path for each dataset
 get_path <- function(w) {
  if (w %in% names(folder_map)) {
   return(file.path(GAlist$dbdir, folder_map[[w]]))
  }
  if (w %in% names(GAlist$dirs)) {
   return(GAlist$dirs[[w]])
  }
  file.path(GAlist$dbdir, w)
 }

 logdf <- GAlist$downloads %||% NULL

 results <- data.frame(
  what = expected,
  path = vapply(expected, get_path, FUN.VALUE = character(1)),
  n_files = NA_integer_,
  total_size_MB = NA_real_,
  in_log = FALSE,
  exists = FALSE,
  large_enough = FALSE,
  status = "❌ missing",
  stringsAsFactors = FALSE
 )

 # Loop over each dataset
 for (i in seq_along(expected)) {
  w <- expected[i]
  path <- results$path[i]

  if (!dir.exists(path)) next

  # --- Special case: annotation unpacked ZIP in gsets/ ---
  if (w == "annotation") {
   core <- file.path(path, "map.rds")
   if (file.exists(core)) {
    files <- list.files(path, full.names = TRUE)
    results$n_files[i] <- length(files)
    results$total_size_MB[i] <- sum(file.info(files)$size, na.rm = TRUE) / 1e6
    results$exists[i] <- TRUE
    results$large_enough[i] <- results$total_size_MB[i] >= min_size_MB
   }
   next
  }

  # --- Standard case ---
  files <- list.files(path, full.names = TRUE, recursive = TRUE)
  n <- length(files)
  total <- if (n > 0) sum(file.info(files)$size, na.rm = TRUE) / 1e6 else 0

  results$n_files[i] <- n
  results$total_size_MB[i] <- round(total, 2)
  results$exists[i] <- n > 0
  results$large_enough[i] <- total >= min_size_MB
 }

 # --- Mark datasets with successful log entries ---
 if (!is.null(logdf)) {
  results$in_log <- vapply(
   results$what,
   function(w) any(logdf$what == w & grepl("success", logdf$status, ignore.case = TRUE)),
   FUN.VALUE = logical(1)
  )
 }

 # --- Final status assignment ---
 results$status <- ifelse(
  results$exists & results$large_enough & results$in_log, "✅ complete",
  ifelse(results$exists & results$large_enough, "⚠️ no log entry",
         ifelse(results$exists, "⚠️ too small", "❌ missing"))
 )

 # --- Print summary ---
 if (verbose) {
  cat("\n📦 GACT download status summary:\n")
  print(results[, c("what", "n_files", "total_size_MB", "in_log", "status")], row.names = FALSE)
 }

 invisible(results)
}

# Small helper (NULL coalescing)
`%||%` <- function(x, y) if (!is.null(x)) x else y


#' @export
#'
removeStatDB <- function(GAlist,studyID=NULL) {
 if (is.null(studyID)) stop("Please provide a valid studyID.")
 if (!all(studyID %in% GAlist$study$id))
  stop("One or more study IDs not found in the database.")
 # remove summary statistics files
 for (i in 1:length(studyID)) {
  file.remove(GAlist$studyfiles[studyID[i]])
 }
 # update GAlist
 #study <- as.data.frame(GAlist$study)
 #study <- study[!study$id%in%studyID,]
 #GAlist$study <- as.list(study)
 #studyfiles <- GAlist$studyfiles
 #studyfiles <- studyfiles[!names(studyfiles)%in%studyID]
 #GAlist$studyfiles <- studyfiles

 # Update GAlist metadata
 study <- as.data.frame(GAlist$study)
 study <- study[!study$id%in%studyID,]
 GAlist$study <- as.list(study)

 GAlist$studyfiles <- GAlist$studyfiles[!names(GAlist$studyfiles) %in% studyID]
 GAlist$studies <- study
 rownames(GAlist$studies) <- GAlist$studies$id

 # Update on-disk metadata
 file_stat_information <- file.path(GAlist$dirs["gstat"], "GWAS_information.csv")
 if (file.exists(file_stat_information))
  write.csv2(study, file_stat_information, row.names = FALSE)

 message("Removed studies: ", paste(studyID, collapse = ", "))
 return(GAlist)
}




#' Update and process genetic association study summary statistics
#'
#' This function updates and processes genetic association study summary statistics. The input to the function includes study information, summary statistics, and information about the sample.
#'
#' @param GAlist A list object providing information and infrastructure of the gact database.
#' @param stat A dataframe of summary statistics for each SNP
#' @param source Source of the summary statistics
#' @param trait Trait being studied
#' @param type Type of trait (binary or quantitative)
#' @param gender Gender of study participants
#' @param ancestry Ancestry of study participants
#' @param build Genome build used in the study
#' @param n Sample size
#' @param ncase Number of cases
#' @param ncontrol Number of controls
#' @param reference Reference for the study
#' @param comments Additional comments
#' @param writeStatDB A flag indicating whether to perform quality control and write the processed summary statistics to a file
#' @param excludeMAFDIFF A threshold used to exclude SNPs with a difference in minor allele frequency (MAF) between cases and controls that exceeds this threshold
#' @return A list of updated study information, including the study id, file name, trait, type, gender, ancestry, build, sample size, case/control count, reference, and comments.
#'
#' @export
#'
updateStatDB <- function(GAlist=NULL,
                         stat=NULL,
                         source=NULL,
                         trait="unknown",
                         type="unknown",
                         gender="unknown",
                         ancestry="unknown",
                         build="unknown",
                         n=NULL,
                         ncase=NULL,
                         ncontrol=NULL,
                         reference="unknown",
                         comments="none",
                         writeStatDB=TRUE) {

 # Update study information only - useful if we need to add extra information
 message("Collecting information on external summary statistics")
 study_number <- length(GAlist$study$id)+1

 #studyID <- paste0("GWAS",study_number)

 # Extract numeric part of all current IDs - works if some studies have been removed
 existing_ids <- GAlist$study$id
 existing_nums <- as.integer(gsub("\\D", "", existing_ids))
 next_id_num <- ifelse(length(existing_nums) == 0, 1, max(existing_nums, na.rm = TRUE) + 1)
 studyID <- paste0("GWAS", next_id_num)

 GAlist$study$id[study_number] <- studyID
 GAlist$study$file[study_number] <- paste0(studyID,".txt")
 GAlist$study$trait[study_number] <- trait
 GAlist$study$type[study_number] <- type
 GAlist$study$gender[study_number] <- gender
 GAlist$study$ancestry[study_number] <- ancestry
 GAlist$study$build[study_number] <- build
 GAlist$study$n[study_number] <- n
 GAlist$study$ncase[study_number] <- ncase
 GAlist$study$ncontrol[study_number] <- ncontrol
 GAlist$study$comments[study_number] <- comments

 if(type=="binary") {
  ntotal <- ncase + ncontrol
  pcase <- ncase/(ncase+ncontrol)
  GAlist$study$neff[study_number] <- ntotal*pcase*(1-pcase)
 }
 if(type=="quantitative") GAlist$study$neff[study_number] <- n

 GAlist$study$reference[study_number] <- reference
 GAlist$study$source[study_number] <- source

 file_stat_information <- file.path(GAlist$dirs["gstat"],"GWAS_information.csv")
 write.csv2(as.data.frame(GAlist$study),file=file_stat_information,row.names=FALSE)
 GAlist$studyfiles <- file.path(GAlist$dirs["gstat"],paste0(GAlist$study$file,".gz"))
 names(GAlist$studyfiles) <- GAlist$study$id


 # processing of summary statistics
 if(writeStatDB) {
  message("Perform quality control of external summary statistics")

  if(any(!sapply(stat, class)[c("eaf", "b", "seb", "p")]%in%"numeric")) warning("Column classes wrong")
  stat$p <- as.numeric(stat$p)
  stat$p[stat$p==0] <- .Machine$double.xmin

  stat <- qcStatDB(GAlist=GAlist,stat=stat, excludeMAF=0.5, excludeMAFDIFF=0.5,
                   excludeINFO=0.8, excludeCGAT=TRUE, excludeINDEL=TRUE,
                   excludeDUPS=TRUE, excludeMHC=FALSE, excludeMISS=0.05,
                   excludeHWE=1e-12)
  message(paste("Writing processed summary statistics to internal file:",
                GAlist$study$file[study_number]))
  file_stat <- GAlist$studyfiles[study_number]
  if(file.exists(file_stat)) stop(paste("GWAS summary statistics file already exists:",
                                        file_stat))
  fwrite(stat, file_stat, col.names=TRUE)
 }
 GAlist$studies <- as.data.frame(GAlist$study)
 rownames(GAlist$studies) <- GAlist$studies$id
 return(GAlist)
}

#' @export
#'
qcStatDB <- function(GAlist=NULL, stat=NULL, excludeMAF=0.01, excludeMAFDIFF=0.05,
                     excludeINFO=0.8, excludeCGAT=TRUE, excludeINDEL=TRUE,
                     excludeDUPS=TRUE, excludeMHC=FALSE, excludeMISS=0.05,
                     excludeHWE=1e-12) {

 # we use cpra to link sumstats and Glist
 cpra <- as.vector(fread(file = file.path(GAlist$dirs["marker"], "cpra.txt.gz"), header=FALSE, data.table=FALSE))[[1]]
 if(cpra[1]=="cpra") cpra <- cpra[-1]    # this is a fix
 rsids <- GAlist$rsids

 # stat is a data.frame
 if(!is.data.frame(stat)) stop("stat should be  a data frame")
 if(!is.null(stat$rsids)) rownames(stat) <- stat$rsids

 if (!"chr" %in% names(stat) || !"pos" %in% names(stat))
  stop("Input `stat` must contain 'chr' and 'pos' columns for mapping.")


 # Remove rows with missing values
 stat <- na.omit(stat)

 # Handle zero values in p column
 stat$p <- as.numeric(stat$p)
 stat$p[stat$p == 0] <- .Machine$double.xmin


 # Internal summary statistic column format
 # data.frame(rsids, chr, pos, a1, a2, af, b, seb, stat, p, n)     (single trait)
 # list(marker=(rsids, chr, pos, a1, a2, af), b, seb, stat, p, n)  (multiple trait)

 fm_internal <- c("marker","chr","pos","ea","nea","eaf","b","seb","p","n")

 format <- "unknown"

 if(all(fm_internal%in%colnames(stat))) format <- "internal"
 #if(all(fm_internal[1:9]%in%colnames(stat))) format <- "internal"
 if(all(fm_internal[1:5]%in%colnames(stat))) format <- "internal"

 if(format=="unknown") {
  message("Column headings for stat object not found")
  message("Column headings for stat object should be:")
  print(fm_internal)
  stop("please revised your stat object accordingly")
 }

 cpra1 <- paste(stat[,"chr"],stat[,"pos"],toupper(stat[,"ea"]),toupper(stat[,"nea"]),sep="_")
 cpra2 <- paste(stat[,"chr"],stat[,"pos"],toupper(stat[,"nea"]),toupper(stat[,"ea"]),sep="_")

 mapped <- cpra1%in%cpra | cpra2%in%cpra
 message("Map markers based on cpra")
 message(paste("Number of markers in stat mapped to marker ids in GAlist:",sum(mapped)))
 message(paste("Number of markers in stat not mapped to marker ids in GAlist:",sum(!mapped)))

 stat <- stat[mapped,]
 cpra1 <- cpra1[mapped]
 cpra2 <- cpra2[mapped]
 rws1 <- match(cpra1,cpra)
 rws2 <- match(cpra2,cpra)

 stat$marker[!is.na(rws1)] <- rsids[rws1[!is.na(rws1)]]
 stat$marker[!is.na(rws2)] <- rsids[rws2[!is.na(rws2)]]

 isdup <- duplicated(stat$marker)
 if(any(isdup)) message("Removing markers with duplicated ids")
 if(any(isdup)) message(paste("Number of markers duplicated in stat:",sum(isdup)))
 stat <- stat[!isdup,]
 rownames(stat) <- stat$marker

 marker <- fread(file.path(GAlist$dirs["marker"],"markers.txt.gz"),
                 data.table=FALSE)

 rownames(marker) <- marker$rsids

 #message("Filtering markers based on information in GAlist:")
 #message("")

 #message("Filtering markers based on information in stat:")
 #message("")

 if(!is.null(stat$marker)) marker_in_stat <- marker$rsids%in%stat$marker
 message(paste("Number of markers in stat also found in bimfiles:", sum(marker_in_stat)))
 message("")
 if(sum(marker_in_stat)==0) stop("No marker ids found in bimfiles")

 # align marker and stat object
 marker <- marker[marker_in_stat,]
 stat <- stat[marker$rsids,]

 aligned <- stat$ea==marker$ea
 message(paste("Number of effect alleles aligned with first allele in bimfiles:", sum(aligned)))
 message(paste("Number of effect alleles not aligned with first allele in bimfiles:", sum(!aligned)))
 message("")

 #original
 effect <- stat[,"b"]
 effect_allele <- stat[,"ea"]
 non_effect_allele <- stat[,"nea"]
 if(!is.null(stat$eaf)) effect_allele_freq <- stat[,"eaf"]

 # flip is not aligned
 stat[!aligned,"b"] <- -effect[!aligned]
 stat[!aligned,"ea"] <- non_effect_allele[!aligned]
 stat[!aligned,"nea"] <- effect_allele[!aligned]
 if(!is.null(stat$eaf)) stat[!aligned,"eaf"] <- 1-effect_allele_freq[!aligned]

 if(is.null(stat$eaf)) {
  message("No effect allele frequency (eaf) provided - using eaf in database")
  stat$eaf <- marker$eaf
 }


 # exclude based on maf
 excludeMAFDIFF <- abs(marker$eaf-stat$eaf) > excludeMAFDIFF
 message(paste("Number of markers excluded by large difference between MAF difference:", sum(excludeMAFDIFF)))
 message("")
 stat <- stat[!excludeMAFDIFF,]
 marker <- marker[!excludeMAFDIFF,]

 if(is.null(stat$n)) stat$n <- neff(seb=stat$seb,af=stat$eaf)
 colnames(stat)[1] <- "rsids"
 return(stat)
}



#' @export
#'
columnStatDB <- function(stat=NULL) {
 col_names <- colnames(stat)
 types <- sapply(stat,class)
 isInteger <- types=="integer"
 isNumeric <- types=="numeric"
 isCharacter <- types=="character"
 isChr <- isPos <- isMarker <- NULL

 if(any(isInteger)) {
  nlevs <- sapply(stat[,isInteger],function(x){length(unique(x))})
  maxVal <- sapply(stat[,isInteger],function(x){max(unique(x))})
  isChr <- col_names[isInteger][22 <= nlevs & nlevs <= 24]
  isPos <- col_names[isInteger][maxVal>100000000]
 }
 if(any(isCharacter)) {
  nlevs <- sapply(stat[,isCharacter],function(x){length(unique(x))})
  if(is.null(isChr)) isChr <- col_names[isCharacter][22 <= nlevs & nlevs <= 24]
  #isAllele <- col_names[isCharacter][nlevs<5]
  isAllele <- sapply(stat[,isCharacter],function(x){sum(c("A","T","C","G")%in%toupper(unique(x)))==4})
  isAllele <- col_names[isCharacter][isAllele]
  if(length(isAllele)==1) isEA <- isNEA <- isAllele
  if(length(isAllele)==2) {
   isEA <- isAllele[1]
   isNEA <- isAllele[2]
  }
  isMarker <- col_names[isCharacter][nlevs>10000]
 }
 if(any(isNumeric)) {
  minVal <- sapply(stat[,isNumeric],min, na.rm=TRUE)
  maxVal <- sapply(stat[,isNumeric],max, na.rm=TRUE)
  meanVal <- sapply(stat[,isNumeric],mean, na.rm=TRUE)
  medianVal <- sapply(stat[,isNumeric],median, na.rm=TRUE)
  isB <- col_names[isNumeric][minVal<0]
  isSEB <- col_names[isNumeric][0<minVal & maxVal<0.2]
  isP <- col_names[isNumeric][meanVal>0 & minVal<0.00001 & meanVal<0.4]
  isN <- col_names[isNumeric][minVal>1000]
  isNca <- isN[1]
  isNco <- isN[1]
  isN <- isN[1]
  if(max(maxVal)>100000000) isPos <- col_names[isNumeric][maxVal>100000000 && minVal>0]
  isEAF <- col_names[isNumeric][0<=minVal & maxVal<1 & medianVal>0.1 & medianVal<0.7]
  if(length(isEAF)==0) isEAF <- "Unknown"
  isInfo <- col_names[isNumeric][0.3<minVal & maxVal<=1 & medianVal>0.7]
  if(length(isInfo)==0) isInfo <- "Unknown"
 }
 best_guess <- c(isMarker,isChr,isPos, isEA, isNEA, isEAF,
                 isB, isSEB, isP, isN, isNca, isNco, isInfo)
 #names(best_guess) <- c("marker","chr","pos","ea","nea","eaf",
 #                      "b","seb","p","n","ncase","ncontrols","info")
 best_guess <- as.character(na.omit(best_guess))
 best_guess
}

#' Get data from a Zenodo archive

#' Download files from a Zenodo archive using a DOI
#' Original version authors see below - modified by Peter Sørensen
#'
#' @author Hans Van Calster, \email{hans.vancalster@@inbo.be}
#' @author Floris Vanderhaeghe, \email{floris.vanderhaeghe@@inbo.be}
#'
#' @param doi The DOI of the Zenodo record
#' @param path Destination path for downloaded files (default: current directory)
#' @param parallel Logical, whether to download files in parallel (default: TRUE)
#' @param quiet Logical, whether to suppress messages (default: FALSE)
#' @importFrom stringr str_remove str_split str_match
#' @importFrom curl curl_fetch_memory curl_download
#' @importFrom jsonlite fromJSON
#' @importFrom tools md5sum
#' @importFrom utils tail
#' @importFrom assertthat assert_that is.string is.flag noNA
#' @noRd
download_zenodo <- function(doi, path = ".", parallel = FALSE, quiet = FALSE, max_tries = 3) {
 if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Please install the 'jsonlite' package.")
 if (!requireNamespace("curl", quietly = TRUE)) stop("Please install the 'curl' package.")

 if (!is.character(doi) || length(doi) != 1) stop("Argument 'doi' must be a single character string.")
 if (!is.character(path) || length(path) != 1) stop("Argument 'path' must be a single character string.")
 if (!dir.exists(path)) stop("Directory does not exist: ", path)
 if (!is.logical(parallel) || length(parallel) != 1) stop("'parallel' must be TRUE or FALSE.")
 if (!is.logical(quiet) || length(quiet) != 1) stop("'quiet' must be TRUE or FALSE.")

 record_id <- sub("^10\\.5281/zenodo\\.", "", doi)
 zenodo_api_url <- paste0("https://zenodo.org/api/records/", record_id)
 response <- curl::curl_fetch_memory(zenodo_api_url)
 content <- jsonlite::fromJSON(rawToChar(response$content))

 file_urls <- content$files$links$self
 filenames <- basename(content$files$key)
 destfiles <- file.path(path, filenames)
 file_md5s <- vapply(strsplit(content$files$checksum, ":"), function(x) x[2], character(1))
 file_sizes <- content$files$size

 total_size <- sum(file_sizes)
 num_files <- length(filenames)

 message("Zenodo record contains ", num_files, " file(s), total size: ",
         format(structure(total_size, class = "object_size")), "\n")

 safe_download <- function(url, destfile, expected_md5, quiet) {
  for (i in seq_len(max_tries)) {
   try({
    curl::curl_download(url, destfile, quiet = quiet)
    md5 <- tolower(trimws(unname(tools::md5sum(destfile))))
    if (identical(md5, expected_md5)) {
     if (!quiet) message(basename(destfile), " download complete and verified (md5: ", md5, ")")
     return(TRUE)
    }
   }, silent = TRUE)
   if (file.exists(destfile)) file.remove(destfile)
   if (!quiet) message("Retrying ", basename(destfile), " (attempt ", i, " of ", max_tries, ")")
   Sys.sleep(2 ^ i)
  }
  stop("Download failed after ", max_tries, " attempts: ", basename(destfile))
 }

 if (parallel && length(file_urls) > 1) {
  message("Downloading in parallel...")
  for (i in seq_along(file_urls)) {
   safe_download(file_urls[i], destfiles[i], file_md5s[i], quiet)
  }
 } else {
  message("Downloading sequentially...")
  for (i in seq_along(file_urls)) {
   safe_download(file_urls[i], destfiles[i], file_md5s[i], quiet)
  }
 }

 message("\nAll files downloaded and verified successfully.")
}

# download_zenodo <- function(doi, path = ".", parallel = FALSE, quiet = FALSE) {
#  # Validate input arguments
#  assert_that(is.string(doi), is.string(path))
#  assert_that(is.flag(parallel), noNA(parallel), is.flag(quiet), noNA(quiet))
#
#  # Ensure the directory exists
#  stopifnot(dir.exists(path))
#
#  # Remove DOI prefix and fetch record details from Zenodo
#  record_id <- str_remove(doi, "10.5281/zenodo.")
#  zenodo_api_url <- paste0("https://zenodo.org/api/records/", record_id)
#  response <- curl_fetch_memory(zenodo_api_url)
#  content <- fromJSON(rawToChar(response$content))
#
#  # Extract file information
#  file_urls <- content$files$links$self
#  filenames <- basename(content$files$key)
#  destfiles <- file.path(path, filenames)
#  file_md5s <- content$files$checksum
#
#
#  # Calculate total file size and number of files
#  total_size <- sum(content$files$size)
#  num_files <- length(filenames)
#
#   # Download files (either in parallel or sequentially)
#  if (parallel && length(file_urls) > 1) {
#   curl::multi_download(urls = file_urls, destfiles = destfiles, progress = !quiet)
#  } else {
#   mapply(curl_download, file_urls, destfiles, MoreArgs = list(quiet = quiet))
#  }
#
#  # Verify file integrity
#  if (!quiet) message("\nVerifying file integrity...\n")
#  for (i in seq_along(file_urls)) {
#   expected_md5 <- str_split(file_md5s[i], ":")[[1]][2]
#   verify_file_integrity(filenames[i], destfiles[i], expected_md5, quiet)
#  }
# }

#' Verify the integrity of a downloaded file
#'
#' @param filename Name of the file
#' @param destfile Path to the downloaded file
#' @param expected_md5 Expected MD5 checksum
#' @param quiet Logical, whether to suppress messages
#' @return None
#' @noRd

verify_file_integrity <- function(filename, destfile, expected_md5, quiet) {

 computed_md5 <- unname(md5sum(destfile))

 # Ensure both MD5 strings are in the same format
 computed_md5 <- tolower(trimws(computed_md5))
 expected_md5 <- tolower(trimws(expected_md5))

 if (identical(computed_md5, expected_md5)) {
  if (!quiet) {
   message(filename, " was downloaded and its integrity verified (md5sum: ", computed_md5, ")")
  }
 } else {
  warning("Incorrect download! md5sum ", computed_md5, " for file ", filename, " does not match the expected md5sum ", expected_md5)
 }
}


#' Summarize GSEA and BLR results
#'
#' @param GAlist A list containing directory paths and potentially other information
#' @return None
#' @export

summaryDB <- function(GAlist = NULL) {
 if (is.null(GAlist)) {
  stop("GAlist is required.")
 }

 processGSEA(GAlist)
 processBLR(GAlist)
}

processGSEA <- function(GAlist) {
 message("Processing GSEA result")

 # Read GSEA results and remove NAs
 resEUR <- readRDS(file.path(GAlist$dirs["gsea"], "res_vegas_eur_filtered.rds"))
 resEAS <- na.omit(readRDS(file.path(GAlist$dirs["gsea"], "res_vegas_eas_filtered.rds")))
 resSAS <- na.omit(readRDS(file.path(GAlist$dirs["gsea"], "res_vegas_sas_filtered.rds")))

 # Determine common genes
 genes <- Reduce(intersect, list(rownames(resEUR$Z), rownames(resEAS), rownames(resSAS)))

 # Combine and save results
 saveCombinedResults(resEUR, resEAS, resSAS, genes, GAlist$dirs["gsea"])
}

processBLR <- function(GAlist) {
 message("Processing BLR result")

 # Process each study
 studyIDs <- getStudyIDs(GAlist$dirs["gbayes"])
 gbayesfiles <- list.files(GAlist$dirs["gbayes"], pattern = ".rds", full.names = TRUE)
 processStudies(gbayesfiles, studyIDs, GAlist$dirs["gbayes"])
}

getStudyIDs <- function(dirPath) {
 studyIDs <- list.files(dirPath, pattern = ".rds")
 studyIDs <- gsub("fit_blr_pruned_|.rds", "", studyIDs)
 return(studyIDs)
}

processStudies <- function(gbayesfiles, studyIDs, dirPath) {
 # Iterate over studies and save results
 for (i in seq_along(gbayesfiles)) {
  fit <- readRDS(gbayesfiles[i])
  saveStudyResults(fit, studyIDs[i], dirPath)
 }
}

saveCombinedResults <- function(resEUR, resEAS, resSAS, genes, dirPath) {
 # Combining Z and P values
 Z <- cbind(resEUR$Z[genes,], GWAS7 = resEAS[genes, "Z"], GWAS8 = resSAS[genes, "Z"])
 P <- cbind(resEUR$p[genes,], GWAS7 = resEAS[genes, "p"], GWAS8 = resSAS[genes, "p"])

 # Save the combined results
 saveRDS(Z, file = file.path(dirPath, "Z_vegas_filtered.rds"))
 saveRDS(P, file = file.path(dirPath, "P_vegas_filtered.rds"))

 # Additional processing and saving data if required...
}

saveStudyResults <- function(fit, studyID, dirPath) {
 # Saving different statistical results
 fnstat <- file.path(dirPath, paste0(studyID, "_stat_BayesC.txt.gz"))
 fwrite(fit$stat[!fit$stat$bm == 0,], file = fnstat)

 fnpost <- file.path(dirPath, paste0(studyID, "_post_BayesC.txt.gz"))
 fwrite(cbind(fit$post, fit$conv[, -4]), file = fnpost)

}



# https://www.medrxiv.org/content/10.1101/2020.09.08.20190561v1
# Promoter capture Hi-C
# Enhancer-promoter correlation
# ABC-Max


# if(what=="pharmgkb") {
#  cwd <- getwd()
#  setwd(GAlist$dirs["drugdb"])
#  url <- "https://api.pharmgkb.org/v1/download/file/data/drugLabels.zip"
#  output_file <-  file.path(GAlist$dirs["drugdb"],"drugLabels.zip")
#  download.file(url, destfile = output_file, mode = "wb")
#  unzip(output_file)
#  file.remove(output_file)
#  url <- "https://api.pharmgkb.org/v1/download/file/data/relationships.zip"
#  output_file <-  file.path(GAlist$dirs["drugdb"],"relationships.zip")
#  download.file(url, destfile = output_file, mode = "wb")
#  unzip(output_file)
#  file.remove(output_file)
#  url <- "https://api.pharmgkb.org/v1/download/file/data/clinicalVariants.zip"
#  output_file <-  file.path(GAlist$dirs["drugdb"],"clinicalVariants.zip")
#  download.file(url, destfile = output_file, mode = "wb")
#  unzip(output_file)
#  file.remove(output_file)
#  url <- "https://api.pharmgkb.org/v1/download/file/data/automated_annotations.zip"
#  output_file <-  file.path(GAlist$dirs["drugdb"],"automated_annotations.zip")
#  download.file(url, destfile = output_file, mode = "wb")
#  unzip(output_file)
#  file.remove(output_file)
#  setwd(cwd)
# }

# if(what=="opentargets") {
#  dir.create(file.path(GAlist$dirs["drugdb"], "targets"))
#  dir.create(file.path(GAlist$dirs["drugdb"], "diseases"))
#  dir.create(file.path(GAlist$dirs["drugdb"], "diseaseToPhenotype"))
#  dir.create(file.path(GAlist$dirs["drugdb"], "significantAdverseDrugReactions"))
#  dir.create(file.path(GAlist$dirs["drugdb"], "associationByDatasourceDirect"))
#  dir.create(file.path(GAlist$dirs["drugdb"], "associationByDatatypeDirect"))
#  dir.create(file.path(GAlist$dirs["drugdb"], "associationByOverallDirect"))
#
#  # Download opentarget
#  urls <- c("ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.11/output/etl/json/targets/",
#            "ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.11/output/etl/json/diseases/",
#            "ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.11/output/etl/json/diseaseToPhenotype/",
#            "ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.11/output/etl/json/fda/significantAdverseDrugReactions/",
#            "ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.11/output/etl/json/associationByDatasourceDirect/",
#            "ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.11/output/etl/json/associationByDatatypeDirect/",
#            "ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.11/output/etl/json/associationByOverallDirect/")
#
#  for (i in 1:length(urls)) {
#   urldir <- unlist(strsplit(urls[i],split="/"))
#   urldir <- urldir[length(urldir)]
#   files = getURL(urls[i], ftp.use.epsv = FALSE, dirlistonly = TRUE)
#   files <- stringr::str_c(urls[i], stringr::str_split(files, "\n")[[1]])
#   files <- stringr::str_trim(files)
#   is.json <- grep(".json",files, fixed=TRUE)
#   files <- files[is.json]
#   for (j in 1:length(files)) {
#    destfile <- unlist(strsplit(files[j],split="/"))
#    destfile <- destfile[length(destfile)]
#    destfile <- file.path(file.path(GAlist$dirs["drugdb"], urldir),destfile)
#    download.file(url=files[j], mode = "wb", dest=destfile)
#   }
#  }
#  files <- dir(file.path(GAlist$dirs["drugdb"], "associationByOverallDirect"), full.names = TRUE)
#  df <- NULL
#  for (i in 1:length(files)) {
#   con <- file(files[i],"r")
#   df <- rbind(df,jsonlite::stream_in(con))
#   close(con)
#  }
#  filename <- file.path(GAlist$dirs["drugdb"], "associationByOverallDirect.tsv")
#  fwrite(df, file=filename)
#
#  files <- dir(file.path(GAlist$dirs["drugdb"], "associationByDatatypeDirect"), full.names = TRUE)
#  df <- NULL
#  for (i in 1:length(files)) {
#   con <- file(files[i],"r")
#   df <- rbind(df,jsonlite::stream_in(con))
#   close(con)
#  }
#  filename <- file.path(GAlist$dirs["drugdb"], "associationByDatatypeDirect.tsv")
#  fwrite(df, file=filename)
#
#  files <- dir(file.path(GAlist$dirs["drugdb"], "associationByDatasourceDirect"), full.names = TRUE)
#  df <- NULL
#  for (i in 1:length(files)) {
#   con <- file(files[i],"r")
#   df <- rbind(df,jsonlite::stream_in(con))
#   close(con)
#  }
#  filename <- file.path(GAlist$dirs["drugdb"], "associationByDatasourceDirect.tsv")
#  fwrite(df, file=filename)
# }

# if(what=="alphamissense") {
#  url <- "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg19.tsv.gz"
#  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://storage.googleapis.com/dm_alphamissense/","",url))
#  download.file(url=url, mode = "wb", dest=destfile)
#  url <- "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
#  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://storage.googleapis.com/dm_alphamissense/","",url))
#  download.file(url=url, mode = "wb", dest=destfile)
# }

#
# drug2eg <- split( drugdb$entrez_id, f=as.factor(drugdb$drug_name) )
# drug2eg <- lapply(drug2eg,function(x){as.character(x)})
# drug2eg <- drug2eg[!names(drug2eg)==""]
#
# drug2ensg <- lapply(drug2eg,function(x){na.omit(unlist(GAlist$map$eg2ensg[x]))})
# drug2ensg <- lapply(drug2ensg, function(x){unique(x)})
# drug2ensg <- drug2ensg[sapply(drug2ensg, function(x){ !any(is.na(x)) } )]
# drug2ensg <- drug2ensg[ sapply(drug2ensg, length)>0]
#
# #saveRDS(drug2ensg,file=file.path(GAlist$dirs["gsets"],"drug2ensg.rds"))
# saveRDS(drug2ensg,file=file.path(GAlist$dirs["gsets"],"drugGenes.rds"))
#


#' Create database for genomic association of complex traits (gact)
#'
#' This function sets up a directory structure for a genomic database and initializes
#' a list to manage the database. It is intended for internal use.
#'
#' @param version The name/version of the database.
#' @param dbdir The root directory for the database.
#' @return A list providing information and infrastructure of the genomic database.
#' @noRd

createDB_old <- function(version = NULL, dbdir = NULL) {
 # Validate required inputs
 if (is.null(version)) stop("Please include a database name using the version argument")

 # Normalize and check directory path
 dbdir <- normalizePath(file.path(dbdir, version), winslash = "/", mustWork = FALSE)
 if (dir.exists(dbdir)) stop(paste("Directory:", dbdir, "already exists"))

 # Create database directory structure
 dir.create(dbdir)
 dirnames <- c("glist", "gstat", "gsets", "gsea", "gbayes", "gtex", "gwas",
               "ldsc", "marker", "drugdb", "download", "script", "vep")
 dirs <- setNames(file.path(dbdir, dirnames), dirnames)
 lapply(dirs, dir.create)

 # Initialize GAlist
 features <- c("Markers", "Genes", "Proteins", "GO", "Pathways", "ProteinComplexes", "ChemicalComplexes", "GTEx", "GWAS")
 GAlist <- list(version = version, dbdir = dbdir, dirs = dirs, features = features)

 return(GAlist)
}


#' Download database components for genomic associations
#'
#' This function facilitates downloading various components of a genomic
#' associations database. It allows selective downloading of different
#' sections such as genetic marker information, gene sets, GSEA summary
#' statistics, and other genomic data. The data is fetched from specified
#' sources and stored in a structured format.
#'
#' @param GAlist A list object containing the structure and paths needed for the database.
#' @param what A character string specifying the type of data to download.
#'             Possible values include "db", "gsets", "gsea", "gstat",
#'             "gbayes", "gwascatalog", "ensembl", "reactome", "string",
#'             "stitch", "pubmed", "dgi", "drugbank", "atc", and others.
#' @param min_combined_score Minimum combined score for filtering (applicable if relevant).
#' @param min_interactions Minimum number of interactions for filtering (applicable if relevant).
#'
#' @return The updated GAlist object, including paths to the downloaded files.
#'
#' @examples
#' \dontrun{
#'   GAlist <- downloadDB(GAlist, what = "gsets")
#' }
#'
#' @export
downloadDB_old <- function(GAlist=NULL, what=NULL, min_combined_score=900,  min_interactions=5) {

 #options(download.file.method="libcurl", url.method="libcurl", timeout=600)
 if(is.null(what)) stop("Please specify what to download e.g. what=gsets")

 if(what=="annotation") {
  message("Downloading annotation sets")
  download_zenodo(doi = "10.5281/zenodo.14234857", path=GAlist$dirs["gsets"], parallel = FALSE)
  dest <- file.path(GAlist$dirs["gsets"],"hsa.0.0.1.zip")
  Sys.sleep(1)
  unzip(dest, exdir=GAlist$dirs["gsets"])
  GAlist$map <- readRDS(file.path(GAlist$dirs["gsets"],"map.rds"))
 }

 # if(what=="gsets") {
 #  message("Downloading annotation sets")
 #  download_zenodo(doi = "10.5281/zenodo.10462983", path=GAlist$dirs["gsets"])
 #  gsetsNames <- list.files(file.path(GAlist$dirs["gsets"]), pattern=".rds")
 #  GAlist$mapfiles <- list.files(file.path(GAlist$dirs["gsets"]), pattern=".rds", full.names=TRUE)
 #  names(GAlist$mapfiles) <- gsub(".rds","", gsetsNames)
 # }

 if(what=="gsea") {
  message("Downloading gsea summary statistics")
  download_zenodo(doi = "10.5281/zenodo.10462484", path=GAlist$dirs["gsea"], parallel = FALSE)
 }

 if(what=="gstat") {
  message("Downloading GWAS summary statistics")
  download_zenodo(doi = "10.5281/zenodo.10462495", path=GAlist$dirs["gstat"], parallel = FALSE)
  destfile <- file.path(GAlist$dirs["gstat"],"GWAS_information.csv")
  GAlist$study <- as.list(read.csv2(destfile))
  GAlist$studies <- as.data.frame(GAlist$study)
  rownames(GAlist$studies) <- GAlist$studies$id
  GAlist$studyfiles <- list.files(file.path(GAlist$dirs["gstat"]), pattern=".gz", full.names = TRUE)
  studyfiles <- list.files(file.path(GAlist$dirs["gstat"]), pattern=".gz")
  names(GAlist$studyfiles) <- gsub(".txt.gz","",studyfiles)
 }

 if(what=="marker") {
  message("Downloading marker information")
  download_zenodo(doi = "10.5281/zenodo.10462385", path=GAlist$dirs["marker"], parallel = FALSE)
  markers <- fread(file.path(GAlist$dirs["marker"],"markers.txt.gz"),
                   data.table=FALSE)
  GAlist$rsids <- markers$rsids
  cpra <- paste(markers$chr,
                markers$pos,
                markers$ea,
                markers$nea,sep="_")
  fwrite(data.frame(cpra = cpra), col.names = FALSE,
         file = file.path(GAlist$dirs["marker"], "cpra.txt"))
 }

 if(what=="ensembl") {
  message("Downloading ensembl information")
  url <- "https://ftp.ensembl.org/pub/release-109/tsv/homo_sapiens/Homo_sapiens.GRCh38.109.entrez.tsv.gz"
  destfile <- file.path(GAlist$dirs["gsets"],"GRCh38.109.entrez.tsv.gz")
  download.file(url=url, mode = "wb", dest=destfile)

  url <- "https://ftp.ensembl.org/pub/release-109/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz"
  destfile <- file.path(GAlist$dirs["gsets"],"GRCh38.Regulatory_Build.regulatory_features.gff.gz")
  download.file(url=url, mode = "wb", dest=destfile)

  url <- "https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz"
  destfile <- file.path(GAlist$dirs["gsets"],"Homo_sapiens.GRCh38.109.gtf.gz")
  download.file(url=url, mode = "wb", dest=destfile)

  url <- "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
  destfile <- file.path(GAlist$dirs["gsets"],"Homo_sapiens.GRCh38.110.gtf.gz")
  download.file(url=url, mode = "wb", dest=destfile)

  url <- "https://ftp.ensembl.org/pub/grch37/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
  #url <- "https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
  destfile <- file.path(GAlist$dirs["gsets"],"Homo_sapiens.GRCh37.87.gtf.gz")
  download.file(url=url, mode = "wb", dest=destfile)

  #url <- "https://ftp.ensembl.org/pub/grch37/current/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff.gz"
  url <- "https://ftp.ensembl.org/pub/grch37/release-113/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff.gz"
  destfile <- file.path(GAlist$dirs["gsets"],"GRCh37.Regulatory_Build.regulatory_features.gff.gz")
  download.file(url=url, mode = "wb", dest=destfile)

 }

 if(what=="pubmed") {
  message("Downloading files from pubmed")
  url <- "https://ftp.ncbi.nih.gov/gene/DATA/gene2pubmed.gz"
  destfile <- file.path(GAlist$dirs["gsets"],"gene2pubmed.gz")
  download.file(url=url, mode = "wb", dest=destfile)
 }

 if(what=="diseases") {
  message("Downloading files from diseases")
  url <- "https://download.jensenlab.org/human_disease_textmining_full.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://download.jensenlab.org/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  Sys.sleep(1)
  gzip(destfile)

  url <- "https://download.jensenlab.org/human_disease_textmining_filtered.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://download.jensenlab.org/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  Sys.sleep(1)
  gzip(destfile)

  url <- "https://download.jensenlab.org/human_disease_knowledge_full.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://download.jensenlab.org/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  Sys.sleep(1)
  gzip(destfile)

  url <- "https://download.jensenlab.org/human_disease_knowledge_filtered.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://download.jensenlab.org/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  Sys.sleep(1)
  gzip(destfile)

  url <- "https://download.jensenlab.org/human_disease_experiments_full.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://download.jensenlab.org/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  Sys.sleep(1)
  gzip(destfile)

  url <- "https://download.jensenlab.org/human_disease_experiments_filtered.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://download.jensenlab.org/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  Sys.sleep(1)
  gzip(destfile)

  url <- "https://download.jensenlab.org/human_disease_integrated_full.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://download.jensenlab.org/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  Sys.sleep(1)
  gzip(destfile)
 }

 if(what=="tiga") {
  url <- "https://unmtid-shinyapps.net/download/TIGA/latest/tiga_gene-trait_provenance.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://unmtid-shinyapps.net/download/TIGA/latest/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  Sys.sleep(1)
  gzip(destfile)

  url <- "https://unmtid-shinyapps.net/download/TIGA/latest/tiga_gene-trait_stats.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://unmtid-shinyapps.net/download/TIGA/latest/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  Sys.sleep(1)
  gzip(destfile)

  url <- "https://unmtid-shinyapps.net/download/TIGA/latest/tiga_genes.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://unmtid-shinyapps.net/download/TIGA/latest/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  Sys.sleep(1)
  gzip(destfile)

  url <- "https://unmtid-shinyapps.net/download/TIGA/latest/tiga_traits.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://unmtid-shinyapps.net/download/TIGA/latest/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  Sys.sleep(1)
  gzip(destfile)
 }

 if(what=="pubchem") {
  #https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/xrefs/PubMedID/TXT
  #https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/
 }

 if(what=="string") {
  message("Downloading files from string")
  url <- "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
  destfile <- file.path(GAlist$dirs["gsets"],"9606.protein.links.v12.0.txt.gz")
  download.file(url=url, mode = "wb", dest=destfile)
 }

 if(what=="stitch") {
  message("Downloading files from stitch")
  url <- "http://stitch.embl.de/download/protein_chemical.links.v5.0/9606.protein_chemical.links.v5.0.tsv.gz"
  destfile <- file.path(GAlist$dirs["gsets"],"9606.protein_chemical.links.v5.0.tsv.gz")
  download.file(url=url, mode = "wb", dest=destfile)
 }

 if(what=="gbayes") {
  message("Downloading gbayes files")
  download_zenodo(doi = "10.5281/zenodo.10462421", path=GAlist$dirs["gbayes"], parallel = FALSE)
 }

 # if(what=="1000G-LDscores") {
 #
 #  # This is a gzipped copy of the European LD scores from 1000 Genomes provided by the Alkes Group,
 #  # originally available at https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2.
 #
 #  download_zenodo(doi = "10.5281/zenodo.8182036", path=GAlist$dirs["marker"])
 #  tarfile <- file.path(GAlist$dirs["marker"],"eur_w_ld_chr.tar.gz")
 #  exdir <- file.path(GAlist$dirs["marker"],"ldscores")
 #  untar(tarfile=tarfile, exdir = exdir)
 #
 #  markers <- NULL
 #  for ( chr in 1:22) {
 #   filename <- file.path(GAlist$dirs["marker"],"ldscores", "eur_w_ld_chr",paste0(chr,".l2.ldscore.gz"))
 #   markers <- rbind(markers,fread(filename, data.table=FALSE))
 #  }
 #  colnames(markers) <- c("chr","rsids","pos","map","eaf", "ldscores")
 #  fwrite(markers, file.path(GAlist$dirs["marker"], "markers_1000G_eur_w_ld.txt.gz"))
 #
 # }


 if(what=="1000G") {
  message("Downloading 1000G files")
  download_zenodo(doi = "10.5281/zenodo.10462403", path=GAlist$dirs["marker"], parallel = FALSE)
  dest <- file.path(GAlist$dirs["marker"],"1000G_EUR_Phase3_plink.zip")
  unzip(dest, exdir=GAlist$dirs["marker"])

  download_zenodo(doi = "10.5281/zenodo.14141438", path=GAlist$dirs["marker"], parallel = FALSE)
  dest <- file.path(GAlist$dirs["marker"],"g1000_eur.zip")
  unzip(dest, exdir=GAlist$dirs["marker"])

  download_zenodo(doi = "10.5281/zenodo.14141523", path=GAlist$dirs["marker"], parallel = FALSE)
  dest <- file.path(GAlist$dirs["marker"],"g1000_eas.zip")
  unzip(dest, exdir=GAlist$dirs["marker"])

  download_zenodo(doi = "10.5281/zenodo.14141597", path=GAlist$dirs["marker"], parallel = FALSE)
  dest <- file.path(GAlist$dirs["marker"],"g1000_sas.zip")
  unzip(dest, exdir=GAlist$dirs["marker"])

  download_zenodo(doi = "10.5281/zenodo.14141774", path=GAlist$dirs["marker"], parallel = FALSE)
  dest <- file.path(GAlist$dirs["marker"],"g1000_afr.zip")
  unzip(dest, exdir=GAlist$dirs["marker"])

  download_zenodo(doi = "10.5281/zenodo.14141885", path=GAlist$dirs["marker"], parallel = FALSE)
  dest <- file.path(GAlist$dirs["marker"],"g1000_amr.zip")
  unzip(dest, exdir=GAlist$dirs["marker"])

  #url <- "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip"
  #url <- "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_afr.zip"
  #url <- "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eas.zip"
  #url <- "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_sas.zip"
  #url <- "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_amr.zip"

 }

 if(what=="vep") {
  #message("Downloading vep files")
  message("Not implemented yet")

  # rsids <- GAlist$rsids
  # rs2vep_list <- list()  # Store data frames temporarily
  #
  # # Define columns to retain
  # cls <- c("#Uploaded_variation", "SIFT", "PolyPhen", "LOEUF", "CADD_PHRED")
  #
  # # Read chromosome-specific VEP files
  # files <- paste0(file.path(GAlist$dbdir, "vep"), "/vep", sprintf("%02d", 1:22), ".txt.gz")
  # for (i in seq_along(files)) {
  #  df <- fread(files[i], data.table = FALSE, na.strings = "-")
  #  rws <- df[,"#Uploaded_variation"] %in% rsids
  #  df <- df[rws,cls]
  #  df <- df[!duplicated(df[,1]),]
  #  rs2vep_list[[i]] <- df
  #  print(paste("Finished processing chromosome", i))
  # }
  #
  # # Combine all chromosome data into one data frame
  # rs2vep <- do.call(rbind, rs2vep_list)
  # colnames(rs2vep) <- c("rsids","sift","polyphen","loeuf","cadd")
  #
  # sift <- gsub(".*\\(([^)]+)\\).*", "\\1", rs2vep[,2])
  # rs2vep[,2] <- as.numeric(sift)
  #
  # polyphen <- gsub(".*\\(([^)]+)\\).*", "\\1", rs2vep[,3])
  # rs2vep[,3] <- as.numeric(polyphen)
  #
  # rownames(rs2vep) <- rs2vep$rsids
  # rs2vep <- rs2vep[,2:5]
  #
  # saveRDS(rs2vep, file.path(GAlist$dirs["gsets"],"rs2vep.rds"))
  #
  # options(download.file.method="libcurl", url.method="libcurl", timeout=3000)
  #
  # url <- "https://ftp.ensembl.org/pub/release-112/variation/vep/homo_sapiens_vep_112_GRCh37.tar.gz"
  # dbdir <- file.path(GAlist$dbdir, "vep")
  # if(!dir.exists(dbdir)) dir.create(dbdir)
  # dest <- file.path(GAlist$dbdir, "vep/homo_sapiens_vep_112_GRCh37.tar")
  # download.file(url=url,dest=dest, mode="wb")
  # untar(tarfile=dest,exdir = dbdir)

 }

 if(what=="clinvar") {
  message("Downloading files from clinvar")
  url <- "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
  destfile <- file.path(GAlist$dirs["gsets"],"variant_summary.txt.gz")
  download.file(url=url, mode = "wb", dest=destfile)
 }


 if(what=="gtex") {
  message("Downloading gtex files")

  options(download.file.method="libcurl", url.method="libcurl", timeout=600)

  url <- "https://storage.googleapis.com/adult-gtex/bulk-qtl/v7/single-tissue-cis-qtl/GTEx_Analysis_v7_eQTL.tar.gz"
  dbdir <- file.path(GAlist$dbdir, "gtex")
  dest <- file.path(GAlist$dbdir, "gtex/GTEx_Analysis_v7_eQTL.tar")
  download.file(url=url,dest=dest, mode="wb")
  untar(tarfile=dest,exdir = dbdir)

  url <- "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar"
  dbdir <- file.path(GAlist$dbdir, "gtex")
  dest <- file.path(GAlist$dbdir, "gtex/GTEx_Analysis_v8_eQTL.tar")
  download.file(url=url,dest=dest, mode="wb")
  untar(tarfile=dest,exdir = dbdir)
 }

 if(what=="gwascatalog") {
  message("Downloading gwascatalog files")
  options(download.file.method="libcurl", url.method="libcurl", timeout=600)
  dbdir <- file.path(GAlist$dbdir, "gwas")
  if(!dir.exists(dbdir)) dir.create(dbdir)
  #file_studies <- "http://ftp.ebi.ac.uk/pub/databases/gwas/releases/2023/04/25/gwas-catalog-studies_ontology-annotated.tsv"
  file_studies <- "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-studies.tsv"
  destfile <- file.path(dbdir, "gwas-catalog-studies.tsv")
  download.file(file_studies, destfile = destfile, mode = "wb")
  #file_gwas <- "http://ftp.ebi.ac.uk/pub/databases/gwas/releases/2023/04/25/gwas-catalog-associations_ontology-annotated.tsv"
  file_gwas <- "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated.tsv"

  destfile <- file.path(dbdir, "gwas-catalog-associations_ontology-annotated.tsv")
  download.file(file_gwas, destfile = destfile, mode = "wb")
  #GAlist$gwasfiles <- c("gwas-catalog-studies_ontology-annotated.tsv","gwas-catalog-associations_ontology-annotated.tsv")
 }

 if(what=="reactome") {
  message("Downloading reactome files")

  url <- "https://reactome.org/download/current/ReactomePathways.txt"
  destfile <- file.path(GAlist$dirs["gsets"],"ReactomePathways.txt")
  download.file(url=url, mode = "wb", dest=destfile)

  url <- "https://reactome.org/download/current/Ensembl2Reactome.txt"
  destfile <- file.path(GAlist$dirs["gsets"],"Ensembl2Reactome.txt")
  download.file(url=url, mode = "wb", dest=destfile)

 }


 if(what=="dgi") {
  message("Downloading files from dgidb")
  message("Downloading Drug Gene Interaction database")
  url_db <- "https://www.dgidb.org/data/2023-Dec/interactions.tsv"
  destfile <- file.path(GAlist$dirs["drugdb"],"interactions.tsv")
  download.file(url=url_db, mode = "wb", dest=destfile)
  url_db <- "https://www.dgidb.org/data/2023-Dec/genes.tsv"
  destfile <- file.path(GAlist$dirs["drugdb"],"genes.tsv")
  download.file(url=url_db, mode = "wb", dest=destfile)
  url_db <- "https://www.dgidb.org/data/2023-Dec/drugs.tsv"
  destfile <- file.path(GAlist$dirs["drugdb"],"drugs.tsv")
  download.file(url=url_db, mode = "wb", dest=destfile)
  url_db <- "https://www.dgidb.org/data/2023-Dec/categories.tsv"
  destfile <- file.path(GAlist$dirs["drugdb"],"categories.tsv")
  download.file(url=url_db, mode = "wb", dest=destfile)
 }

 if(what=="drugbank") {
  setwd(GAlist$dirs["drugdb"])
  url <- "https://www.dropbox.com/s/2rf7plojpkmy18z/drugid2drugname.rds?dl=1"
  file_drugbank <- file.path(GAlist$dirs["drugdb"],"drugid2drugname.rds")
  download.file(url, destfile = file_drugbank, mode = "wb")
  GAlist$drugbank <- readRDS(file=file_drugbank)
 }

 if(what=="atc") {
  message("Downloading atc information")
  # Add targets to GAlist
  df <- fread("https://www.dropbox.com/s/n5ehglmhhs0kcue/WHO%20ATC-DDD%202023-03-28.csv?dl=1",data.table = FALSE)
  atc <- NULL
  atc$code <- df$atc_code
  atc$name <- df$atc_name
  names(atc$code) <- atc$name
  names(atc$name) <- atc$code
  # GAlist$atc <- NULL
  # GAlist$atc$code <- df$atc_code
  # GAlist$atc$name <- df$atc_name
  # names(GAlist$atc$code) <- GAlist$atc$name
  # names(GAlist$atc$name) <- GAlist$atc$code
  saveRDS(atc, file.path(GAlist$dirs["drugdb"],"atc.rds"))

  # # Add drug target data frame with ATC information to GAlist
  # drug2ensg <- readRDS(file=file.path(GAlist$dirs["gsets"],"drug2ensg.rds"))
  # #drug2ensg <- getSetsDB(GAlist = GAlist, feature = "DrugGenes")
  # nreps <- sapply(drug2ensg,length)
  # drugs <- rep(names(drug2ensg), times=nreps)
  # ensg <- unlist(drug2ensg)
  # ensg2drug <- split(drugs, f=as.factor(ensg))
  # df <- data.frame(Drug=drugs, Target=ensg)
  #   df$ATC <- rep("Unknown",nrow(df))
  # has_atc <- match(tolower(df$Drug),tolower(GAlist$atc$name))
  # df$ATC[!is.na(has_atc)] <- as.character(GAlist$atc$code[has_atc[!is.na(has_atc)]])
  # GAlist$targets <- df
  # target <- GAlist$targets
  # target <- target[!duplicated(target$Drug),]
  # drug2atc <- target$ATC
  # names(drug2atc) <- target$Drug
  # GAlist$drug2atc <- drug2atc
 }

 return(GAlist)
}
#' Download database components for genomic associations (with logging)
#'
#' This function downloads the various components of a genomic associations
#' database and logs metadata for each operation via an external `log_download()`
#' helper function.
#'
#' @param GAlist A list object containing the structure and paths needed for the database.
#' @param what A character string specifying which data type to download.
#' @param min_combined_score Minimum combined score for filtering STRING/STITCH data.
#' @param min_interactions Minimum number of interactions for filtering STRING/STITCH data.
#'
#' @return The updated GAlist object with file paths and GAlist$downloads updated.
#' @export
downloadDB_v0 <- function(GAlist = NULL, what = NULL,
                          min_combined_score = 900, min_interactions = 5) {

 if (is.null(what)) stop("Please specify what to download, e.g. what='gsets'")

 # Ensure timeout and libcurl settings
 options(download.file.method = "libcurl", url.method = "libcurl", timeout = 3000)

 # --------------------------------------------------------------------
 # Annotation (Zenodo)
 # --------------------------------------------------------------------
 if (what == "annotation") {
  dest <- GAlist$dirs["gsets"]
  doi <- "10.5281/zenodo.14234857"
  message("Downloading annotation sets")
  tryCatch({
   download_zenodo(doi = doi, path = dest, parallel = TRUE)
   zipfile <- file.path(dest, "hsa.0.0.1.zip")
   Sys.sleep(1)
   unzip(zipfile, exdir = dest)
   GAlist$map <- readRDS(file.path(dest, "map.rds"))
   GAlist <- log_download(GAlist, what, "Zenodo", doi, dest, "success", "Annotation sets downloaded")
  }, error = function(e) {
   GAlist <- log_download(GAlist, what, "Zenodo", doi, dest, "failed", e$message)
  })
 }

 # --------------------------------------------------------------------
 # Marker information (Zenodo)
 # --------------------------------------------------------------------
 if (what == "marker") {
  dest <- GAlist$dirs["marker"]
  doi <- "10.5281/zenodo.10462385"
  message("Downloading marker information")
  tryCatch({
   download_zenodo(doi = doi, path = dest, parallel = TRUE)
   markers <- fread(file.path(dest, "markers.txt.gz"), data.table = FALSE)
   GAlist$rsids <- markers$rsids
   cpra <- paste(markers$chr, markers$pos, markers$ea, markers$nea, sep = "_")
   fwrite(data.frame(cpra = cpra), col.names = FALSE, file = file.path(dest, "cpra.txt"))
   GAlist <- log_download(GAlist, what, "Zenodo", doi, dest, "success", "Marker data downloaded")
  }, error = function(e) {
   GAlist <- log_download(GAlist, what, "Zenodo", doi, dest, "failed", e$message)
  })
 }

 # --------------------------------------------------------------------
 # GSEA (Zenodo)
 # --------------------------------------------------------------------
 if (what == "gsea") {
  dest <- GAlist$dirs["gsea"]
  doi <- "10.5281/zenodo.10462484"
  message("Downloading GSEA summary statistics")
  tryCatch({
   download_zenodo(doi = doi, path = dest, parallel = TRUE)
   GAlist <- log_download(GAlist, what, "Zenodo", doi, dest, "success", "GSEA statistics downloaded")
  }, error = function(e) {
   GAlist <- log_download(GAlist, what, "Zenodo", doi, dest, "failed", e$message)
  })
 }

 # --------------------------------------------------------------------
 # GWAS summary statistics (Zenodo)
 # --------------------------------------------------------------------
 if (what == "gstat") {
  dest <- GAlist$dirs["gstat"]
  doi <- "10.5281/zenodo.10462495"
  message("Downloading GWAS summary statistics")
  tryCatch({
   download_zenodo(doi = doi, path = dest, parallel = TRUE)
   meta <- file.path(dest, "GWAS_information.csv")
   GAlist$study <- as.list(read.csv2(meta))
   GAlist$studies <- as.data.frame(GAlist$study)
   rownames(GAlist$studies) <- GAlist$studies$id
   gz <- list.files(dest, pattern = ".gz", full.names = TRUE)
   GAlist$studyfiles <- gz
   names(GAlist$studyfiles) <- gsub(".txt.gz", "", basename(gz))
   GAlist <- log_download(GAlist, what, "Zenodo", doi, dest, "success", "GWAS summary statistics downloaded")
  }, error = function(e) {
   GAlist <- log_download(GAlist, what, "Zenodo", doi, dest, "failed", e$message)
  })
 }

 # --------------------------------------------------------------------
 # Bayesian statistics results (Zenodo)
 # --------------------------------------------------------------------
 if (what == "gbayes") {
  dest <- GAlist$dirs["gbayes"]
  doi <- "10.5281/zenodo.10462421"
  message("Downloading gbayes files")
  tryCatch({
   download_zenodo(doi = doi, path = dest, parallel = TRUE)
   GAlist <- log_download(GAlist, what, "Zenodo", doi, dest, "success", "gbayes results downloaded")
  }, error = function(e) {
   GAlist <- log_download(GAlist, what, "Zenodo", doi, dest, "failed", e$message)
  })
 }

 # --------------------------------------------------------------------
 # Ensembl
 # --------------------------------------------------------------------
 if (what == "ensembl") {
  #dest <- GAlist$dirs["gsets"]
  dest <- GAlist$dirs["ensembl"]
  if (!dir.exists(dest)) dir.create(dest)
  base <- "https://ftp.ensembl.org/pub"
  urls <- c(
   paste0(base, "/release-109/tsv/homo_sapiens/Homo_sapiens.GRCh38.109.entrez.tsv.gz"),
   paste0(base, "/release-109/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz"),
   paste0(base, "/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz"),
   paste0(base, "/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"),
   paste0(base, "/grch37/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"),
   paste0(base, "/grch37/release-113/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff.gz")
  )
  message("Downloading Ensembl GTF and regulatory data")
  tryCatch({
   for (u in urls) download.file(u, file.path(dest, basename(u)), mode = "wb")
   GAlist <- log_download(GAlist, what, "Ensembl FTP", paste(urls, collapse = "; "), dest, "success", "Ensembl files downloaded")
  }, error = function(e) {
   GAlist <- log_download(GAlist, what, "Ensembl FTP", paste(urls, collapse = "; "), dest, "failed", e$message)
  })
 }

 # --------------------------------------------------------------------
 # STRING
 # --------------------------------------------------------------------
 if (what == "string") {
  #dest <- GAlist$dirs["gsets"]
  dest <- GAlist$dirs["string"]
  if (!dir.exists(dest)) dir.create(dest)
  url <- "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
  message("Downloading STRING v12.0 data")
  tryCatch({
   download.file(url, file.path(dest, basename(url)), mode = "wb")
   GAlist <- log_download(GAlist, what, "STRING", url, dest, "success", "STRING v12.0 network downloaded")
  }, error = function(e) {
   GAlist <- log_download(GAlist, what, "STRING", url, dest, "failed", e$message)
  })
 }

 # --------------------------------------------------------------------
 # STITCH
 # --------------------------------------------------------------------
 if (what == "stitch") {
  #dest <- GAlist$dirs["gsets"]
  dest <- GAlist$dirs["stitch"]
  if (!dir.exists(dest)) dir.create(dest)
  url <- "http://stitch.embl.de/download/protein_chemical.links.v5.0/9606.protein_chemical.links.v5.0.tsv.gz"
  message("Downloading STITCH v5.0 data")
  tryCatch({
   download.file(url, file.path(dest, basename(url)), mode = "wb")
   GAlist <- log_download(GAlist, what, "STITCH", url, dest, "success", "STITCH protein–chemical links downloaded")
  }, error = function(e) {
   GAlist <- log_download(GAlist, what, "STITCH", url, dest, "failed", e$message)
  })
 }

 # --------------------------------------------------------------------
 # Reactome
 # --------------------------------------------------------------------
 if (what == "reactome") {
  #dest <- GAlist$dirs["gsets"]
  dest <- GAlist$dirs["reactome"]
  if (!dir.exists(dest)) dir.create(dest)
  urls <- c(
   "https://reactome.org/download/current/ReactomePathways.txt",
   "https://reactome.org/download/current/Ensembl2Reactome.txt"
  )
  message("Downloading Reactome pathway mappings")
  tryCatch({
   for (u in urls) download.file(u, file.path(dest, basename(u)), mode = "wb")
   GAlist <- log_download(GAlist, what, "Reactome", paste(urls, collapse = "; "), dest, "success", "Reactome pathway data downloaded")
  }, error = function(e) {
   GAlist <- log_download(GAlist, what, "Reactome", paste(urls, collapse = "; "), dest, "failed", e$message)
  })
 }

 # --------------------------------------------------------------------
 # DGIdb
 # --------------------------------------------------------------------
 if (what == "dgi") {
  dest <- GAlist$dirs["drugdb"]
  if (!dir.exists(dest)) dir.create(dest)
  base <- "https://www.dgidb.org/data/2023-Dec"
  urls <- paste0(base, "/", c("interactions.tsv", "genes.tsv", "drugs.tsv", "categories.tsv"))
  message("Downloading DGIdb drug–gene interaction data")
  tryCatch({
   for (u in urls) download.file(u, file.path(dest, basename(u)), mode = "wb")
   GAlist <- log_download(GAlist, what, "DGIdb", paste(urls, collapse = "; "), dest, "success", "DGIdb files downloaded")
  }, error = function(e) {
   GAlist <- log_download(GAlist, what, "DGIdb", paste(urls, collapse = "; "), dest, "failed", e$message)
  })
 }

 # --------------------------------------------------------------------
 # GTEx (v7 + v8)
 # --------------------------------------------------------------------
 if (what == "gtex") {
  dest <- GAlist$dirs["gtex"]
  if (!dir.exists(dest)) dir.create(dest)
  urls <- c(
   "https://storage.googleapis.com/adult-gtex/bulk-qtl/v7/single-tissue-cis-qtl/GTEx_Analysis_v7_eQTL.tar.gz",
   "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar"
  )
  message("Downloading GTEx eQTL data (v7 and v8)")
  tryCatch({
   for (u in urls) {
    destfile <- file.path(dest, basename(u))
    download.file(u, destfile, mode = "wb")
    untar(destfile, exdir = dest)
   }
   GAlist <- log_download(GAlist, what, "GTEx", paste(urls, collapse = "; "), dest, "success", "GTEx eQTLs downloaded and extracted")
  }, error = function(e) {
   GAlist <- log_download(GAlist, what, "GTEx", paste(urls, collapse = "; "), dest, "failed", e$message)
  })
 }

 # --------------------------------------------------------------------
 # GWAS Catalog (EBI)
 # --------------------------------------------------------------------
 if (what == "gwascatalog") {
  #dest <- file.path(GAlist$dbdir, "gwas")
  dest <- file.path(GAlist$dbdir, "gwascatalog")
  if (!dir.exists(dest)) dir.create(dest)
  urls <- c(
   "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-studies.tsv",
   #"https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated.tsv"
   "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip",
   "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-ancestry.tsv"
  )
  message("Downloading GWAS Catalog data")
  tryCatch({
   for (u in urls) download.file(u, destfile = file.path(dest, basename(u)), mode = "wb")
   GAlist <- log_download(GAlist, what, "EBI GWAS Catalog", paste(urls, collapse = "; "), dest, "success", "GWAS Catalog data downloaded")
  }, error = function(e) {
   GAlist <- log_download(GAlist, what, "EBI GWAS Catalog", paste(urls, collapse = "; "), dest, "failed", e$message)
  })
 }

 # --------------------------------------------------------------------
 # ATC classification
 # --------------------------------------------------------------------
 if (what == "atc") {
  dest <- GAlist$dirs["drugdb"]
  url <- "https://www.dropbox.com/s/n5ehglmhhs0kcue/WHO%20ATC-DDD%202023-03-28.csv?dl=1"
  message("Downloading WHO ATC classification")
  tryCatch({
   df <- fread(url, data.table = FALSE)
   atc <- list(code = df$atc_code, name = df$atc_name)
   names(atc$code) <- atc$name
   names(atc$name) <- atc$code
   saveRDS(atc, file.path(dest, "atc.rds"))
   GAlist <- log_download(GAlist, what, "WHO ATC", url, dest, "success", "ATC data downloaded")
  }, error = function(e) {
   GAlist <- log_download(GAlist, what, "WHO ATC", url, dest, "failed", e$message)
  })
 }

 # --------------------------------------------------------------------
 # 1000 Genomes reference panels (Zenodo)
 # --------------------------------------------------------------------
 if (what == "1000G") {
  message("Downloading 1000 Genomes reference data")

  dest <- GAlist$dirs["marker"]

  # Define DOIs and expected zip names explicitly
  refs <- data.frame(
   name = c("1000G_EUR_Phase3_plink", "g1000_eur", "g1000_eas",
            "g1000_sas", "g1000_afr", "g1000_amr"),
   doi  = c("10.5281/zenodo.10462403",
            "10.5281/zenodo.14141438",
            "10.5281/zenodo.14141523",
            "10.5281/zenodo.14141597",
            "10.5281/zenodo.14141774",
            "10.5281/zenodo.14141885"),
   stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(refs))) {
   nm  <- refs$name[i]
   doi <- refs$doi[i]

   message(" - Downloading ", nm, " (Zenodo ", doi, ")")

   tryCatch({
    # Download from Zenodo
    files_downloaded <- download_zenodo(doi = doi, path = dest, parallel = FALSE)

    # If download_zenodo() returns files, check what we got
    if (is.character(files_downloaded) && length(files_downloaded) > 0) {
     zipfile <- files_downloaded[grep("\\.zip$", files_downloaded)]
    } else {
     # Fallback: assume file name based on known pattern
     zipfile <- file.path(dest, paste0(nm, ".zip"))
    }

    # Ensure file exists
    if (length(zipfile) > 0 && file.exists(zipfile)) {
     unzip(zipfile, exdir = dest)
     message("   ✓ Unzipped: ", basename(zipfile))
    } else {
     warning("   ⚠ No zip file found for ", nm)
    }

    GAlist <- log_download(
     GAlist,
     what = "1000G",
     source = "Zenodo",
     urls = doi,
     dest_path = dest,
     status = "success",
     message = paste(nm, "reference panel downloaded and unpacked")
    )

   }, error = function(e) {
    warning("Failed to download ", nm, ": ", e$message)
    GAlist <- log_download(
     GAlist,
     what = "1000G",
     source = "Zenodo",
     urls = doi,
     dest_path = dest,
     status = "failed",
     message = paste("Error downloading", nm, ":", e$message)
    )
   })
  }
 }



 # --------------------------------------------------------------------
 # DISEASES database (Jensen Lab)
 # --------------------------------------------------------------------
 if (what == "diseases") {
  message("Downloading DISEASES database files from JensenLab")

  dest <- GAlist$dirs["gsets"]

  urls <- c(
   "https://download.jensenlab.org/human_disease_textmining_full.tsv",
   "https://download.jensenlab.org/human_disease_textmining_filtered.tsv",
   "https://download.jensenlab.org/human_disease_knowledge_full.tsv",
   "https://download.jensenlab.org/human_disease_knowledge_filtered.tsv",
   "https://download.jensenlab.org/human_disease_experiments_full.tsv",
   "https://download.jensenlab.org/human_disease_experiments_filtered.tsv",
   "https://download.jensenlab.org/human_disease_integrated_full.tsv"
  )

  for (url in urls) {
   fname <- basename(url)
   destfile <- file.path(dest, fname)

   message(" - Downloading ", fname)

   tryCatch({
    download.file(url = url, destfile = destfile, mode = "wb")
    Sys.sleep(1)

    # Compress (if not already gzipped)
    if (!grepl("\\.gz$", destfile)) {
     R.utils::gzip(destfile, overwrite = TRUE)
    }

    GAlist <- log_download(
     GAlist,
     what = "diseases",
     source = "JensenLab DISEASES",
     urls = url,
     dest_path = dest,
     status = "success",
     message = paste(fname, "downloaded and compressed")
    )

   }, error = function(e) {
    GAlist <- log_download(
     GAlist,
     what = "diseases",
     source = "JensenLab DISEASES",
     urls = url,
     dest_path = dest,
     status = "failed",
     message = paste("Failed to download", fname, ":", e$message)
    )
   })
  }
 }

 return(GAlist)
}

