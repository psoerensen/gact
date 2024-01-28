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
                 dbdir=NULL) {

 if(is.null(dbdir)) dbdir <- getwd()
 if(is.null(version)) stop("Please provide database version")

 if(task=="download") {
  GAlist <- createDB(version=version, dbdir=dbdir)

  options(download.file.method="libcurl", url.method="libcurl", timeout=1200)

  # Step 2: Download data from database:
  GAlist <- downloadDB(GAlist=GAlist, what="marker")
  GAlist <- downloadDB(GAlist=GAlist, what="gsets")
  GAlist <- downloadDB(GAlist=GAlist, what="gsea")
  GAlist <- downloadDB(GAlist=GAlist, what="gstat")
  GAlist <- downloadDB(GAlist=GAlist, what="gbayes")
  GAlist <- downloadDB(GAlist=GAlist, what="gwascatalog")
  GAlist <- downloadDB(GAlist=GAlist, what="ensembl")
  GAlist <- downloadDB(GAlist=GAlist, what="reactome")
  GAlist <- downloadDB(GAlist=GAlist, what="string")
  GAlist <- downloadDB(GAlist=GAlist, what="stitch")
  GAlist <- downloadDB(GAlist=GAlist, what="pubmed")
  GAlist <- downloadDB(GAlist=GAlist, what="dgi")
  #GAlist <- downloadDB(GAlist=GAlist, what="gtex")
  #GAlist <- downloadDB(GAlist=GAlist, what="1000G")
  #GAlist <- downloadDB(GAlist=GAlist, what="diseases")
  #GAlist <- downloadDB(GAlist=GAlist, what="tiga")
  #GAlist <- downloadDB(GAlist=GAlist, what="pubchem")
  #GAlist <- downloadDB(GAlist=GAlist, what="pharmgkb")
  #GAlist <- downloadDB(GAlist=GAlist, what="opentargets")
  #GAlist <- downloadDB(GAlist=GAlist, what="atc")
  #GAlist <- downloadDB(GAlist=GAlist, what="alphamissense")

  message("Creating full marker sets - this may take some time")
  GAlist <- createSetsDB(GAlist=GAlist)
  #GAlist <- createSetsDB(GAlist=GAlist, what="diseases")
  #GAlist <- createMarkerSetsDB(GAlist=GAlist, what="GO")
  #summaryDB(GAlist=GAlist)

 }
 # Step 3: Create marker sets from database:
 if(task=="createSets") {
  message("Creating full marker sets - this may take some time")
  GAlist <- createSetsDB(GAlist=GAlist)
  GAlist <- createMarkerSetsDB(GAlist=GAlist, what="reactome")
  #GAlist <- createMarkerSetsDB(GAlist=GAlist, what="GO")
  #GAlist <- createMarkerSetsDB(GAlist=GAlist, what="string")
  #GAlist <- createMarkerSetsDB(GAlist=GAlist, what="stitch")
 }
 return(GAlist)
}



#' Create database for genomic association of complex traits (gact)
#'
#' This function sets up a directory structure for a genomic database and initializes
#' a list to manage the database. It is intended for internal use.
#'
#' @param version The name/version of the database.
#' @param dbdir The root directory for the database.
#' @return A list providing information and infrastructure of the genomic database.
#' @noRd

createDB <- function(version = NULL, dbdir = NULL) {
 # Validate required inputs
 if (is.null(version)) stop("Please include a database name using the version argument")

 # Normalize and check directory path
 dbdir <- normalizePath(file.path(dbdir, version), winslash = "/", mustWork = FALSE)
 if (dir.exists(dbdir)) stop(paste("Directory:", dbdir, "already exists"))

 # Create database directory structure
 dir.create(dbdir)
 dirnames <- c("glist", "gstat", "gsets", "gsea", "gbayes", "gtex", "gwas", "ldsc", "marker", "drugdb", "download")
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
downloadDB <- function(GAlist=NULL, what=NULL, min_combined_score=900,  min_interactions=5) {

 options(download.file.method="libcurl", url.method="libcurl", timeout=600)
 if(is.null(what)) stop("Please specify what to download e.g. what=gsets")

 if(what=="db") {
  message("Downloading hsa.0.01.zip")
  download_zenodo(doi = "10.5281/zenodo.10464396", path=GAlist$dbdir)
  dest <- file.path(GAlist$dbdir,"hsa.0.01.zip")
  unzip(dest, exdir=dest)
 }

 if(what=="gsets") {
  message("Downloading marker sets")
  download_zenodo(doi = "10.5281/zenodo.10462983", path=GAlist$dirs["gsets"])
  gsetsNames <- list.files(file.path(GAlist$dirs["gsets"]), pattern=".rds")
  GAlist$gsetsfiles <- list.files(file.path(GAlist$dirs["gsets"]), pattern=".rds", full.names=TRUE)
  names(GAlist$gsetsfiles) <- gsub(".rds","", gsetsNames)
 }

 if(what=="gsea") {
  message("Downloading gsea summary statistics")
  download_zenodo(doi = "10.5281/zenodo.10462485", path=GAlist$dirs["gsea"])
  gseaNames <- list.files(file.path(GAlist$dirs["gsea"]), pattern=".rds")
  GAlist$gseafiles <- list.files(file.path(GAlist$dirs["gsea"]), pattern=".rds", full.names=TRUE)
  names(GAlist$gseafiles) <- gsub(".rds","",gseaNames)
 }

 if(what=="gstat") {
  message("Downloading GWAS summary statistics")
  download_zenodo(doi = "10.5281/zenodo.10462496", path=GAlist$dirs["gstat"])
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
  download_zenodo(doi = "10.5281/zenodo.10467174", path=GAlist$dirs["marker"])
  GAlist$markerfiles <- file.path(GAlist$dirs["marker"],"markers.txt.gz")
  markers <- fread(file.path(GAlist$dirs["marker"],"markers.txt.gz"),
                   data.table=FALSE)
  GAlist$rsids <- markers$rsids
  GAlist$cpra <- paste(markers$chr,
                       markers$pos,
                       markers$ea,
                       markers$nea,sep="_")
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
  R.utils::gzip(destfile)

  url <- "https://download.jensenlab.org/human_disease_textmining_filtered.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://download.jensenlab.org/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  gzip(destfile)

  url <- "https://download.jensenlab.org/human_disease_knowledge_full.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://download.jensenlab.org/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  gzip(destfile)

  url <- "https://download.jensenlab.org/human_disease_knowledge_filtered.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://download.jensenlab.org/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  gzip(destfile)

  url <- "https://download.jensenlab.org/human_disease_experiments_full.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://download.jensenlab.org/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  gzip(destfile)

  url <- "https://download.jensenlab.org/human_disease_experiments_filtered.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://download.jensenlab.org/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  gzip(destfile)

  url <- "https://download.jensenlab.org/human_disease_integrated_full.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://download.jensenlab.org/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  gzip(destfile)
 }

 if(what=="tiga") {
  url <- "https://unmtid-shinyapps.net/download/TIGA/latest/tiga_gene-trait_provenance.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://unmtid-shinyapps.net/download/TIGA/latest/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  gzip(destfile)

  url <- "https://unmtid-shinyapps.net/download/TIGA/latest/tiga_gene-trait_stats.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://unmtid-shinyapps.net/download/TIGA/latest/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  gzip(destfile)

  url <- "https://unmtid-shinyapps.net/download/TIGA/latest/tiga_genes.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://unmtid-shinyapps.net/download/TIGA/latest/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
  gzip(destfile)

  url <- "https://unmtid-shinyapps.net/download/TIGA/latest/tiga_traits.tsv"
  destfile <- file.path(GAlist$dirs["gsets"],gsub("https://unmtid-shinyapps.net/download/TIGA/latest/","",url))
  download.file(url=url, mode = "wb", dest=destfile)
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
  download_zenodo(doi = "10.5281/zenodo.10462421", path=GAlist$dirs["gbayes"])
 }



 if(what=="1000G") {
  message("Downloading 1000G files")
  url <- "https://zenodo.org/api/records/10462403"
  download_zenodo(doi = "10.5281/zenodo.10462403", path=GAlist$dirs["marker"])
  dest <- file.path(GAlist$dirs["marker"],"1000G_EUR_Phase3_plink.zip")
  unzip(dest, exdir=GAlist$dirs["marker"])

  url <- "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip"
  dest <- file.path(GAlist$dirs["marker"],
                    gsub("https://ctg.cncr.nl/software/MAGMA/ref_data/","",url))
  download.file(url=url, mode = "wb", dest=dest)
  unzip(dest, exdir=GAlist$dirs["marker"])

  url <- "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_afr.zip"
  dest <- file.path(GAlist$dirs["marker"],
                    gsub("https://ctg.cncr.nl/software/MAGMA/ref_data/","",url))
  download.file(url=url, mode = "wb", dest=dest)
  unzip(dest, exdir=GAlist$dirs["marker"])

  url <- "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eas.zip"
  dest <- file.path(GAlist$dirs["marker"],
                    gsub("https://ctg.cncr.nl/software/MAGMA/ref_data/","",url))
  download.file(url=url, mode = "wb", dest=dest)
  unzip(dest, exdir=GAlist$dirs["marker"])

  url <- "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_sas.zip"
  dest <- file.path(GAlist$dirs["marker"],
                    gsub("https://ctg.cncr.nl/software/MAGMA/ref_data/","",url))
  download.file(url=url, mode = "wb", dest=dest)
  unzip(dest, exdir=GAlist$dirs["marker"])

  url <- "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_amr.zip"
  dest <- file.path(GAlist$dirs["marker"],
                    gsub("https://ctg.cncr.nl/software/MAGMA/ref_data/","",url))
  download.file(url=url, mode = "wb", dest=dest)
  unzip(dest, exdir=GAlist$dirs["marker"])

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
  file_studies <- "http://ftp.ebi.ac.uk/pub/databases/gwas/releases/2023/04/25/gwas-catalog-studies_ontology-annotated.tsv"
  destfile <- file.path(dbdir, "gwas-catalog-studies_ontology-annotated.tsv")
  download.file(file_studies, destfile = destfile, mode = "wb")
  file_gwas <- "http://ftp.ebi.ac.uk/pub/databases/gwas/releases/2023/04/25/gwas-catalog-associations_ontology-annotated.tsv"
  destfile <- file.path(dbdir, "gwas-catalog-associations_ontology-annotated.tsv")
  download.file(file_gwas, destfile = destfile, mode = "wb")
  GAlist$gwasfiles <- c("gwas-catalog-studies_ontology-annotated.tsv","gwas-catalog-associations_ontology-annotated.tsv")
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
  GAlist$atc <- NULL
  GAlist$atc$code <- df$atc_code
  GAlist$atc$name <- df$atc_name
  names(GAlist$atc$code) <- GAlist$atc$name
  names(GAlist$atc$name) <- GAlist$atc$code

  # Add drug target data frame with ATC information to GAlist
  drug2ensg <- getSetsDB(GAlist = GAlist, feature = "DrugGenes")
  nreps <- sapply(drug2ensg,length)
  drugs <- rep(names(drug2ensg), times=nreps)
  ensg <- unlist(drug2ensg)
  ensg2drug <- split(drugs, f=as.factor(ensg))
  df <- data.frame(Drug=drugs, Target=ensg)
    df$ATC <- rep("Unknown",nrow(df))
  has_atc <- match(tolower(df$Drug),tolower(GAlist$atc$name))
  df$ATC[!is.na(has_atc)] <- as.character(GAlist$atc$code[has_atc[!is.na(has_atc)]])
  GAlist$targets <- df
  target <- GAlist$targets
  target <- target[!duplicated(target$Drug),]
  drug2atc <- target$ATC
  names(drug2atc) <- target$Drug
  GAlist$drug2atc <- drug2atc
 }

 return(GAlist)
}



#' @export
#'
createSetsDB <- function(GAlist = NULL, what="ensembl",
                         upstream=35, downstream=10,
                         min_combined_score=900, min_interactions=5) {

 # default sets (order of execution is important)

 GAlist$gsets <- vector(mode = "list", length = length(GAlist$gsetsfiles))
 for(i in 1:length(GAlist$gsetsfiles)) {
  #GAlist$gsets[[i]] <- readRDS(file.path(GAlist$dirs["gsets"],GAlist$gsetsfiles[i]))
  GAlist$gsets[[i]] <- readRDS(GAlist$gsetsfiles[i])
 }
 names(GAlist$gsets) <- names(GAlist$gsetsfiles)
 GAlist$gsets <- GAlist$gsets[c("eg2ensg", "ensg2eg", "ensg2sym", "ensp2ensg")]

 # Ensembl genes, proteins, transcripts
 file <- file.path(GAlist$dirs["gsets"],"GRCh38.109.entrez.tsv.gz")
 ensembl <- fread(file, data.table=FALSE)
 ensembl <- ensembl[!ensembl$protein_stable_id=="-",]

 GAlist$gsets$eg2ensg <- split( ensembl$gene_stable_id, f=as.factor(ensembl$xref) )
 GAlist$gsets$eg2ensp <- split( ensembl$protein_stable_id, f=as.factor(ensembl$xref) )
 GAlist$gsets$eg2enst <- split( ensembl$transcript_stable_id, f=as.factor(ensembl$xref) )

 GAlist$gsets$ensp2eg <- split( ensembl$xref, f=as.factor(ensembl$protein_stable_id) )
 GAlist$gsets$ensp2ensg <- split( ensembl$gene_stable_id, f=as.factor(ensembl$protein_stable_id) )
 GAlist$gsets$ensp2enst <- split( ensembl$transcript_stable_id, f=as.factor(ensembl$protein_stable_id) )

 GAlist$gsets$ensg2eg <- split( ensembl$xref, f=as.factor(ensembl$gene_stable_id) )
 GAlist$gsets$ensg2ensp <- split( ensembl$protein_stable_id, f=as.factor(ensembl$gene_stable_id) )
 GAlist$gsets$ensg2enst <- split( ensembl$transcript_stable_id, f=as.factor(ensembl$gene_stable_id) )

 saveRDS(GAlist$gsets$eg2ensg, file = file.path(GAlist$dirs["gsets"], "eg2ensg.rds"))
 saveRDS(GAlist$gsets$eg2ensp, file = file.path(GAlist$dirs["gsets"], "eg2ensp.rds"))
 saveRDS(GAlist$gsets$eg2enst, file = file.path(GAlist$dirs["gsets"], "eg2enst.rds"))

 saveRDS(GAlist$gsets$ensg2eg, file = file.path(GAlist$dirs["gsets"], "ensg2eg.rds"))
 saveRDS(GAlist$gsets$ensg2ensp, file = file.path(GAlist$dirs["gsets"], "ensg2ensp.rds"))
 saveRDS(GAlist$gsets$ensg2enst, file = file.path(GAlist$dirs["gsets"], "ensg2enst.rds"))

 saveRDS(GAlist$gsets$ensp2eg, file = file.path(GAlist$dirs["gsets"], "ensp2eg.rds"))
 saveRDS(GAlist$gsets$ensp2ensg, file = file.path(GAlist$dirs["gsets"], "ensp2ensg.rds"))
 saveRDS(GAlist$gsets$ensp2enst, file = file.path(GAlist$dirs["gsets"], "ensp2enst.rds"))

 # Gene marker sets
 # Specify parameters
 #upstream <- 35
 #downstream <- 10

 markers <- fread(file.path(GAlist$dirs["marker"],"markers.txt.gz"),
                  data.table=FALSE)

 # file <- file.path(GAlist$dirs["gsets"],"Homo_sapiens.GRCh38.109.gtf.gz")
 # #file <- file.path(GAlist$dirs["gsets"],"Homo_sapiens.GRCh38.110.gtf.gz")
 # df <- fread(file, data.table=FALSE,skip=1, header=FALSE)
 # colnames(df) <- c("chr","source","type","start","end","score","strand","phase","attributes")
 # df <- df[df$type=="gene" & df$source=="ensembl_havana",]
 # att <- strsplit(df$attributes, ";")
 # att <- lapply(att, function(x){gsub("\"","",x)})
 # gene_id <- sapply(att, function(x){ x[grep("gene_id",x)]})
 # df$gene_id <- gsub("gene_id ","",gene_id)
 # df <- df[,c("gene_id","chr","source", "type", "start", "end","strand")]
 # df <- df[!df$chr=="X",]
 # df <- df[!df$chr=="Y",]
 # df$chr <- as.integer(df$chr)
 #
 # upstream <- upstream*1000
 # downstream <- downstream*1000
 # ensg2rsids <- vector("list", nrow(df))
 # ensg2cpra <- vector("list", nrow(df))
 #
 # start <- df$start-upstream
 # start[start<1] <- 1
 # end <- df$end+downstream
 # maxpos <- max(markers$pos,end)
 # pos <- 1:maxpos
 # ensg2rsids <- vector("list", nrow(df))
 # for (chr in 1:22) {
 #  message(paste("Processing chr:",chr))
 #  rsids <- rep(NA, maxpos)
 #  rsids[as.integer(markers[markers$chr==chr,"pos"])] <- markers[markers$chr==chr,"rsids"]
 #  for (i in 1:nrow(df)) {
 #   if(df$chr[i]==chr) {
 #    grsids <- rsids[start[i]:end[i]]
 #    ensg2rsids[[i]] <- grsids[!is.na(grsids)]
 #   }
 #  }
 # }
 # names(ensg2rsids) <- df$gene_id
 # empty <- sapply(ensg2rsids, function(x){ identical(x, character(0))})
 # ensg2rsids <- ensg2rsids[!empty]
 # #setsfile <- file.path(GAlist$dirs["gsets"], "GRCh38.110.ensg2rsids.rds")
 # setsfile <- file.path(GAlist$dirs["gsets"], "ensg2rsids.rds")
 # saveRDS(ensg2rsids, file = setsfile)

 ensg2rsids <- readRDS(file = file.path(GAlist$dirs["gsets"], "ensg2rsids.rds"))
 sets <- mapSetsDB(sets=ensg2rsids, featureID=markers$rsids, index=TRUE)
 chr <- sapply(sets, function(x) {unique(markers$chr[x])})
 start <- sapply(sets, function(x) {min(markers$pos[x])})
 stop <- sapply(sets, function(x) {max(markers$pos[x])})
 df <- data.frame(EnsemblID=names(sets),chr=chr,start=start,stop=stop)
 saveRDS(df, file = file.path(GAlist$dirs["gsets"], "genesplus_annotation.rds"))

 # GWAS catalog
 dbdir <- file.path(GAlist$dbdir, "gwas")
 gwasfile <- file.path(dbdir, "gwas-catalog-associations_ontology-annotated.tsv")
 gwas <- fread(gwasfile, data.table=FALSE, quote="")
 sets <- split(gwas$SNPS,f=gwas$MAPPED_TRAIT)
 saveRDS(sets, file = file.path(GAlist$dirs["gsets"], "gwasSets.rds"))


 # Pubmed to genes
 file <- file.path(GAlist$dirs["gsets"],"gene2pubmed.gz")
 pubmed <- fread(file, data.table = FALSE)
 pubmed <- pubmed[pubmed[,1]%in%9606,-1]
 eg2pmid <- split(pubmed$PubMed_ID,f=as.factor(pubmed$GeneID))
 pmid2eg <- split(pubmed$GeneID,f=as.factor(pubmed$PubMed_ID))
 saveRDS(eg2pmid,file=file.path(GAlist$dirs["gsets"],"eg2pmid.rds"))
 saveRDS(pmid2eg,file=file.path(GAlist$dirs["gsets"],"pmid2eg.rds"))

 #String
 file <- file.path(GAlist$dirs["gsets"],"9606.protein.links.v12.0.txt.gz")
 string <- fread(file, data.table=FALSE)
 string$protein1 <- gsub("9606.","",string$protein1)
 string$protein2 <- gsub("9606.","",string$protein2)
 string  <- string[string$combined_score>=min_combined_score,]
 string <- string[string$protein2%in%names(GAlist$gsets$ensp2ensg),]
 string <- split( string$protein2,f=as.factor(string$protein1))
 sets <- string[sapply(string ,length)>=min_interactions]
 sets <- lapply(sets,function(x){na.omit(unlist(GAlist$gsets$ensp2ensg[x]))})
 saveRDS(sets,file=file.path(GAlist$dirs["gsets"],"string2ensg.rds"))

 #Stitch
 file <- file.path(GAlist$dirs["gsets"],"9606.protein_chemical.links.v5.0.tsv.gz")
 stitch <- fread(file, data.table=FALSE)
 stitch$protein <- gsub("9606.","",stitch$protein)
 stitch <- stitch[stitch$protein%in%names(GAlist$gsets$ensp2ensg),]
 stitch  <- stitch[stitch$combined_score>=min_combined_score,]
 stitch <- split( stitch$protein,f=as.factor(stitch$chemical))
 sets  <- stitch[sapply(stitch ,length)>=min_interactions]
 sets <- lapply(sets,function(x){na.omit(unlist(GAlist$gsets$ensp2ensg[x]))})
 saveRDS(sets,file=file.path(GAlist$dirs["gsets"],"stitch2ensg.rds"))

 # Reactome
 file <- file.path(GAlist$dirs["gsets"],"Ensembl2Reactome.txt")
 if(file.exists(file)) {
  reactome <- fread(file, data.table=FALSE, header=FALSE)
  isHSA <- grep("R-HSA",reactome[,2])
  reactome <- reactome[isHSA,]
  reac2ensg <- split(reactome[,1],f=reactome[,2])
  ensg2reac <- split(reactome[,2],f=reactome[,1])
  GAlist$gsets$reac2ensg <- reac2ensg
  GAlist$gsets$ensg2reac <- ensg2reac
  saveRDS(reac2ensg,file=file.path(GAlist$dirs["gsets"],"reactome2ensg.rds"))
  saveRDS(ensg2reac,file=file.path(GAlist$dirs["gsets"],"ensg2reactome.rds"))
 }
 file <- file.path(GAlist$dirs["gsets"],"ReactomePathways.txt")
 if(file.exists(file)) {
  pathway <- fread(file, data.table=FALSE, header=FALSE)
  isHSA <- grep("R-HSA",pathway[,1])
  pathway <- pathway[isHSA,]
  reac2names <- pathway[,2]
  names(reac2names) <- pathway[,1]
  GAlist$gsets$reac2names <- reac2names
  saveRDS(reac2names,file=file.path(GAlist$dirs["gsets"],"reactome2names.rds"))
 }

 # Regulatory elements
 file <- file.path(GAlist$dirs["gsets"],"GRCh38.Regulatory_Build.regulatory_features.gff.gz")
 df <- fread(file, data.table=FALSE)
 colnames(df) <- c("chr","source","type","start","end","score","strand","phase","attributes")
 df <- df[df$chr%in%as.character(1:22),]
 #df <- df[!df$chr=="X",]
 #df <- df[!df$chr=="Y",]
 df$chr <- as.integer(df$chr)
 df <- df[!is.na(df$chr),]
 att <- strsplit(df$attributes, ";")
 att <- lapply(att, function(x){gsub("\"","",x)})
 att <- sapply(att, function(x){ x[grep("ID=",x)]})
 att <- strsplit(att, ":")
 df$reg_id <- sapply(att, function(x){x[2]})
 rownames(df) <- df$reg_id
 saveRDS(df[,c("reg_id", "type", "chr", "start", "end")], file = file.path(GAlist$dirs["gsets"], "regulatory_annotation.rds"))

 regSets <- split(df$reg_id, f=as.factor(df$type))
 saveRDS(regSets, file = file.path(GAlist$dirs["gsets"], "regSets.rds"))

 #markers <- fread(GAlist$markerfiles, data.table=FALSE)
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

 reg2rsids <- sapply(regSets, function(x){unique(unlist(ensr2rsids[x]))})
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


 # Drug databases
 drugdb <- fread(file.path(GAlist$dirs["drugdb"], "interactions.tsv"),
                 quote = "\"", data.table = FALSE)
 hgnc <- fread("https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt", data.table=FALSE)
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
 GAlist <- downloadDB(GAlist=GAlist, what="atc")

 string2ensg <- readRDS(file=file.path(GAlist$dirs["gsets"],"string2ensg.rds"))
 drug2ensp <- lapply(drug2ensg,function(x){na.omit(unlist(GAlist$gsets$ensg2ensp[x]))})
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

 # drug2string2rsids <- lapply(drug2string2ensg,function(x){unique(unlist(ensg2rsids[x]))})
 # drug2string2rsids <- drug2string2rsids[!sapply(drug2string2rsids,is.null)]
 # saveRDS(drug2string2rsids,file=file.path(GAlist$dirs["gsets"],"drug2string2rsids.rds"))

 if("diseases"%in%what) {
  ensp <- names(GAlist$gsets$ensp2ensg)
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
   disease2ensg <- lapply(sets,function(x) {unique(unlist(GAlist$gsets$ensp2ensg[x]))})
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
  target <- GAlist$targets
  target <- target[!duplicated(target$Drug),]
  drug2atc <- target$ATC
  names(drug2atc) <- target$Drug
  GAlist$drug2atc <- drug2atc
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
createMarkerSetsDB <- function(GAlist = NULL, what=NULL,
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
 return(GAlist)
}

#' @export
#'
removeStatDB <- function(GAlist,studyID=NULL) {
 if (!any(studyID%in%GAlist$study$id)) stop("studyID not data base")
 # remove summary statistics files
 for (i in 1:length(studyID)) {
  file.remove(GAlist$studyfiles[studyID[i]])
 }
 # update GAlist
 study <- as.data.frame(GAlist$study)
 study <- study[!study$id%in%studyID,]
 GAlist$study <- as.list(study)
 studyfiles <- GAlist$studyfiles
 studyfiles <- studyfiles[!names(studyfiles)%in%studyID]
 GAlist$studyfiles <- studyfiles
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
                         writeStatDB=TRUE,
                         excludeMAFDIFF=0.05) {

 # Update study information only - useful if we need to add extra information
 message("Collecting information on external summary statistics")
 study_number <- length(GAlist$study$id)+1
 studyID <- paste0("GWAS",study_number)
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

  stat <- qcStatDB(GAlist=GAlist,stat=stat, excludeMAFDIFF=excludeMAFDIFF)
  message(paste("Writing processed summary statistics to internal file:",
                GAlist$study$file[study_number]))
  file_stat <- GAlist$studyfiles[study_number]
  if(file.exists(file_stat)) stop(paste("GWAS summary statistics file allready exists:",
                                        file_stat))
  fwrite(stat, file_stat)
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
 cpra <- GAlist$cpra
 rsids <- GAlist$rsids

 # stat is a data.frame
 if(!is.data.frame(stat)) stop("stat should be  a data frame")
 if(!is.null(stat$rsids)) rownames(stat) <- stat$rsids

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

download_zenodo <- function(doi, path = ".", parallel = TRUE, quiet = FALSE) {
 # Validate input arguments
 assert_that(is.string(doi), is.string(path))
 assert_that(is.flag(parallel), noNA(parallel), is.flag(quiet), noNA(quiet))

 # Ensure the directory exists
 stopifnot(dir.exists(path))

 # Remove DOI prefix and fetch record details from Zenodo
 record_id <- str_remove(doi, "10.5281/zenodo.")
 zenodo_api_url <- paste0("https://zenodo.org/api/records/", record_id)
 response <- curl_fetch_memory(zenodo_api_url)
 content <- fromJSON(rawToChar(response$content))

 # Extract file information
 file_urls <- content$files$links$self
 filenames <- basename(content$files$key)
 destfiles <- file.path(path, filenames)
 file_md5s <- content$files$checksum


 # Calculate total file size and number of files
 total_size <- sum(content$files$size)
 num_files <- length(filenames)

  # Download files (either in parallel or sequentially)
 if (parallel && length(file_urls) > 1) {
  curl::multi_download(urls = file_urls, destfiles = destfiles, progress = !quiet)
 } else {
  mapply(curl_download, file_urls, destfiles, MoreArgs = list(quiet = quiet))
 }

 # Verify file integrity
 if (!quiet) message("\nVerifying file integrity...\n")
 for (i in seq_along(file_urls)) {
  expected_md5 <- str_split(file_md5s[i], ":")[[1]][2]
  verify_file_integrity(filenames[i], destfiles[i], expected_md5, quiet)
 }
}

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


#' summaryDB <- function(GAlist=NULL) {
#'
#'
#'  ################################################################################
#'  # Summarize gsea results
#'  ################################################################################
#'
#'  message("Processing GSEA result")
#'  resEUR <- readRDS(file=file.path(GAlist$dirs["gsea"],"res_vegas_eur_filtered.rds"))
#'  resEAS <- readRDS(file=file.path(GAlist$dirs["gsea"],"res_vegas_eas_filtered.rds"))
#'  resEAS <- na.omit(resEAS)
#'  resSAS <- readRDS(file=file.path(GAlist$dirs["gsea"],"res_vegas_sas_filtered.rds"))
#'  resSAS <- na.omit(resSAS)
#'  genes <- intersect(intersect(rownames(resEUR$Z),rownames(resEAS)),rownames(resSAS))
#'
#'  Z <- cbind(resEUR$Z[genes,],GWAS7=resEAS[genes,"Z"],GWAS8=resSAS[genes,"Z"])
#'  P <- cbind(resEUR$p[genes,],GWAS7=resEAS[genes,"p"],GWAS8=resSAS[genes,"p"])
#'
#'  saveRDS(Z,file=file.path(GAlist$dirs["gsea"],"Z_vegas_filtered.rds"))
#'  saveRDS(P,file=file.path(GAlist$dirs["gsea"],"P_vegas_filtered.rds"))
#'
#'
#'  df <- addAnnotationDB(df=P, hyperlinkEXCEL=TRUE)
#'  file.csv <- file.path(GAlist$dirs["gsea"],"summary_vegas_filtered.csv")
#'  write.csv2(df,file=file.csv,row.names=FALSE)
#'
#'  df <- addAnnotationDB(df=P)
#'  filename <- file.path(GAlist$dirs["gsea"],"summary_vegas_filtered.txt.gz")
#'  fwrite(df,file=filename,row.names=FALSE)
#'
#'
#'  ################################################################################
#'  # Summarize gbayes results
#'  ################################################################################
#'
#'
#'  message("Processing BLR result")
#'  studyIDs <- list.files(file.path(GAlist$dirs["gbayes"]), pattern=".rds")
#'  studyIDs <- gsub("fit_blr_pruned_","",studyIDs)
#'  studyIDs <- gsub(".rds","",studyIDs)
#'  gbayesfiles <- list.files(file.path(GAlist$dirs["gbayes"]), pattern=".rds", full.names=TRUE)
#'
#'  for (i in 1:length(gbayesfiles)) {
#'   fit <- readRDS(file=gbayesfiles[i])
#'   fnstat <- file.path(GAlist$dirs["gbayes"], paste0(studyIDs[i],"_stat_BayesC.txt.gz"))
#'   fwrite(fit$stat[!fit$stat$bm==0,], file=fnstat)
#'   fnpost <- file.path(GAlist$dirs["gbayes"], paste0(studyIDs[i],"_post_BayesC.txt.gz"))
#'   fwrite(cbind(fit$post,fit$conv[,-4]), file=fnpost)
#'  }
#'
#'  # Extract gene-marker sets
#'  sets <- getMarkerSetsDB(GAlist = GAlist, feature = "Genesplus")
#'  bm <- dm <- vm <- matrix(0,nrow=length(sets),ncol=length(gbayesfiles))
#'  rownames(bm) <- rownames(dm) <- rownames(vm) <- names(sets)
#'  colnames(bm) <- colnames(dm) <- colnames(vm) <- studyIDs
#'
#'  for (i in 1:length(gbayesfiles)) {
#'   fit <- readRDS(file=gbayesfiles[i])
#'   isets <- qgg:::mapSets(sets,fit$stat$rsids,index=TRUE)
#'   b <- sapply(isets,function(x){sum(abs(fit$stat$bm[x]))})
#'   d <- sapply(isets,function(x){sum(abs(fit$stat$dm[x]))})
#'   vb <- sapply(isets,function(x){sum(abs(fit$stat$vm[x]))})
#'   bm[names(b),i] <- b
#'   dm[names(d),i] <- d
#'   vm[names(vb),i] <- vb
#'  }
#'
#'  df <- addAnnotationDB(df=bm, hyperlinkEXCEL=TRUE)
#'  file.csv <- file.path(GAlist$dirs["gbayes"],"summary_b.csv")
#'  write.csv2(df,file=file.csv,row.names=FALSE)
#'
#'  df <- addAnnotationDB(df=dm, hyperlinkEXCEL=TRUE)
#'  file.csv <- file.path(GAlist$dirs["gbayes"],"summary_d.csv")
#'  write.csv2(df,file=file.csv,row.names=FALSE)
#'
#'  df <- addAnnotationDB(df=vm, hyperlinkEXCEL=TRUE)
#'  file.csv <- file.path(GAlist$dirs["gbayes"],"summary_vb.csv")
#'  write.csv2(df,file=file.csv,row.names=FALSE)
#'
#'
#'  df <- addAnnotationDB(df=bm)
#'  filename <- file.path(GAlist$dirs["gbayes"],"summary_b.txt.gz")
#'  fwrite(df,file=filename,row.names=FALSE)
#'
#'  df <- addAnnotationDB(df=dm)
#'  filename <- file.path(GAlist$dirs["gbayes"],"summary_d.txt.gz")
#'  fwrite(df,file=filename,row.names=FALSE)
#'
#'  df <- addAnnotationDB(df=vm)
#'  filename <- file.path(GAlist$dirs["gbayes"],"summary_vb.txt.gz")
#'  fwrite(df,file=filename,row.names=FALSE)
#'
#' }




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
# drug2ensg <- lapply(drug2eg,function(x){na.omit(unlist(GAlist$gsets$eg2ensg[x]))})
# drug2ensg <- lapply(drug2ensg, function(x){unique(x)})
# drug2ensg <- drug2ensg[sapply(drug2ensg, function(x){ !any(is.na(x)) } )]
# drug2ensg <- drug2ensg[ sapply(drug2ensg, length)>0]
#
# #saveRDS(drug2ensg,file=file.path(GAlist$dirs["gsets"],"drug2ensg.rds"))
# saveRDS(drug2ensg,file=file.path(GAlist$dirs["gsets"],"drugGenes.rds"))
#
