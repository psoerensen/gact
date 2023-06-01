#######################################################################################
# Create GACT database functions
#######################################################################################
#'
#' Create or download GACT data
#'
#' @description
#' The gact function is used to create or download a GACT database with information
#' on genomic associations for complex traits. The genomic associations
#' include single marker associations linked to genetic markers in a Glist.
#' This includes stringent quality control. In addition genomic associations
#' linked to different genomic features (e.g. genes, proteins, chemical-complexes,
#' protein-complexes, biological pathways) are provided.
#'
#' @param GAlist A list object providing information and infrastructure of the gact database.
#' @param version The version of the gact database to use.
#' @param task The task to perform, either "download" or "process". By default "download".
#' @param dbdir The directory where the database should be stored, by default the current working directory.
#' @param what Specifies the type of data to download, either "lite" or "full". By default "lite".
#'
#' @return A list providing information and infrastructure of the gact database.
#'
#' @import data.table
#'
#' @examples
#' \dontrun{
#' GAlist <- gact()
#' }


#' @export
#'
gact <- function(GAlist=NULL, version=NULL, task="download",
                 dbdir=NULL, what="lite") {

 if(is.null(dbdir)) dbdir <- getwd()
 if(is.null(version)) stop("Please provide database version")

 if(task=="download") {
  GAlist <- createDB(Glist=NULL, version=version, dbdir=dbdir)

  options(download.file.method="libcurl", url.method="libcurl", timeout=600)

  # Step 2: Download data from database:
  GAlist <- downloadDB(GAlist=GAlist, what="marker")
  GAlist <- downloadDB(GAlist=GAlist, what="gsets")
  #GAlist <- downloadDB(GAlist=GAlist, what="gsea")
  GAlist <- downloadDB(GAlist=GAlist, what="gstat")
  GAlist <- downloadDB(GAlist=GAlist, what="ensembl")
  GAlist <- downloadDB(GAlist=GAlist, what="reactome")
  GAlist <- downloadDB(GAlist=GAlist, what="string")
  GAlist <- downloadDB(GAlist=GAlist, what="stitch")
  #GAlist <- downloadDB(GAlist=GAlist, what="pubmed")
  #GAlist <- downloadDB(GAlist=GAlist, what="pubchem")
  GAlist <- downloadDB(GAlist=GAlist, what="dgi")
  #GAlist <- downloadDB(GAlist=GAlist, what="pharmgkb")
  #GAlist <- downloadDB(GAlist=GAlist, what="opentargets")
  GAlist <- downloadDB(GAlist=GAlist, what="gwascatalog")
  GAlist <- downloadDB(GAlist=GAlist, what="atc")
  #GAlist$study <- as.data.frame(GAlist$study)
 }
 # Step 3: Create marker sets from database:
 if(what=="createSets") {
  message("Creating full marker sets - this may take some time")
  GAlist <- createSetsDB(GAlist=GAlist)
 }

 return(GAlist)
}


#' createDB - Creates a database for genetic association studies
#'
#' @param Glist A list of genetic data. Default is NULL.
#' @param version Character string indicating the name of the database.
#'               It is a required argument.
#' @param dbdir Character string indicating the directory where the database should be created.
#'             Default is NULL.
#' @param what Character string indicating the type of database.
#'             Default is "lite".
#'
#' @return A list containing the following items:
#'         - version: Name of the database
#'         - traits: NULL
#'         - dirs: List of subdirectories within the database directory
#'         - features: List of features in the database
#'         - markers: Data frame of markers in the database
#'         - rsids: List of rsids in the database
#'         - cpra: List of cpra in the database
#'
#' @details
#' The createDB function creates a database for genetic association studies.
#' If `dbdir` is not specified, the database will be created in the current working directory.
#' The function creates subdirectories within the database directory, including "glist", "gstat",
#' "gsets", "gsea", "ldsc", "gbayes", "marker", "download", and "drugdb".
#' The list `Glist` should contain the genetic data and will be filtered to only keep
#' markers that are present in `Glist$rsidsLD`. The filtered markers are stored in a data frame
#' and saved in the "marker" subdirectory as "markers.txt.gz".
#'
#' @examples
#' \dontrun{
#' createDB(Glist = genetic_data, version = "db1", dbdir = "/data/db1")
#' createDB(Glist = genetic_data, version = "db2", dbdir = "/data/db2", what = "full")
#' }
#'
#' @export
#'
createDB <- function(Glist=NULL, version=NULL, dbdir=NULL, what="lite") {
 if (is.null(version)) stop("Please include a database name using the version argument")

 dbdir <- file.path(dbdir, version)
 if (dir.exists(dbdir)) stop(paste("Directory:",dbdir,"already exists"))
 dir.create(dbdir)

 dirs <- c(glist = "glist",
           gstat = "gstat",
           gsets = "gsets",
           gsea = "gsea",
           gbayes = "gbayes",
           ldsc = "ldsc",
           marker = "marker",
           drugdb = "drugdb",
           download = "download")

 lapply(names(dirs), function(x) {
  dir.create(file.path(dbdir, dirs[x]))
 })

 GAlist <- list(version = version,
                traits = NULL,
                dirs = file.path(dbdir, dirs),
                features = c("Markers", "Genes", "Proteins", "GO", "Pathways", "ProteinComplexes", "ChemicalComplexes"))
 GAlist$dbdir <- dbdir

 names(GAlist$dirs) <- c("glist","gstat","gsets","gsea", "gbayes", "ldsc", "marker", "drugdb", "download")

 if(!is.null(Glist)) {
  keep <- Glist$rsids %in% Glist$rsidsLD
  GAlist$markers <- data.frame(rsids = Glist$rsids[keep],
                               chr = Glist$chr[keep],
                               pos = Glist$pos[keep],
                               ea = Glist$a1[keep],
                               nea = Glist$a2[keep],
                               eaf = Glist$af[keep],
                               stringsAsFactors = FALSE)

  GAlist$rsids <- Glist$rsids[keep]
  GAlist$cpra <- Glist$cpra[keep]

  file_markers <- file.path(GAlist$dirs["marker"], "markers.txt.gz")
  fwrite(GAlist$markers, file = file_markers)
 }

 return(GAlist)
}


#' @export
#'
createSetsDB <- function(GAlist = NULL) {
 ensg2rsids <- GAlist$gsets[["ensg2rsids_10kb"]]

 fset_go <- getSetsDB(GAlist = GAlist, feature = "GO")
 sets_go <- lapply(fset_go, function(x){unique(unlist(ensg2rsids[x]))})
 sets_go <- sets_go[!sapply(sets_go, is.null)]
 setsfile_go <- file.path(GAlist$dirs["gsets"], "go2rsids.rds")
 saveRDS(sets_go, file = setsfile_go)

 fset_reactome <- getSetsDB(GAlist = GAlist, feature = "Pathways2Genes")
 sets_reactome <- lapply(fset_reactome, function(x){unique(unlist(ensg2rsids[x]))})
 sets_reactome <- sets_reactome[!sapply(sets_reactome, is.null)]
 setsfile_reactome <- file.path(GAlist$dirs["gsets"], "reactome2rsids.rds")
 saveRDS(sets_reactome, file = setsfile_reactome)

 fset_string <- getSetsDB(GAlist = GAlist, feature = "ProteinComplexes2Genes")
 sets_string <- lapply(fset_string, function(x){unique(unlist(ensg2rsids[x]))})
 sets_string <- sets_string[!sapply(sets_string, is.null)]
 setsfile_string <- file.path(GAlist$dirs["gsets"], "string2rsids.rds")
 saveRDS(sets_string, file = setsfile_string)

 fset_stitch <- getSetsDB(GAlist = GAlist, feature = "ChemicalComplexes2Genes")
 sets_stitch <- lapply(fset_stitch, function(x){unique(unlist(ensg2rsids[x]))})
 sets_stitch <- sets_stitch[!sapply(sets_stitch, is.null)]
 setsfile_stitch <- file.path(GAlist$dirs["gsets"], "stitch2rsids.rds")
 saveRDS(sets_stitch, file = setsfile_stitch)

 GAlist$gsetsfiles[12] <- setsfile_go
 GAlist$gsetsfiles[13] <- setsfile_reactome
 GAlist$gsetsfiles[14] <- setsfile_string
 GAlist$gsetsfiles[15] <- setsfile_stitch

 names(GAlist$gsetsfiles[12]) <- "go2rsids"
 names(GAlist$gsetsfiles[13]) <- "reactome2rsids"
 names(GAlist$gsetsfiles[14]) <- "string2rsids"
 names(GAlist$gsetsfiles[15]) <- "stitch2rsids"
 return(GAlist)
}

#' @export
#'
downloadDB <- function(GAlist=NULL, what=NULL, min_combined_score=900,  min_interactions=5) {


 if(is.null(what)) stop("Please specify what to download e.g. what=gsets")

 if(what=="gsets") {
  # download gsets files in the database
  message("Downloading annotation and marker sets")
  urls <- c("https://www.dropbox.com/s/ijtc7l6hgpaieo1/eg2rsids_10kb.rds?dl=1",
            "https://www.dropbox.com/s/0aqbqa7ihrg6i2e/ensg2rsids_10kb.rds?dl=1",
            "https://www.dropbox.com/s/p3ut5dwfx0zw4v1/ensp2rsids_10kb.rds?dl=1",
            "https://www.dropbox.com/s/1py37zd92ttsvnp/ensg2sym.rds?dl=1",
            "https://www.dropbox.com/s/2ggu4u5hp406cif/go.rds?dl=1",
            "https://www.dropbox.com/s/uryyxnjyhxa9azf/reactome.rds?dl=1",
            "https://www.dropbox.com/s/wnci7lldztnb93k/string2ensp.rds?dl=1",
            "https://www.dropbox.com/s/ny94ibdbqhtg62h/stitch.rds?dl=1",
            "https://www.dropbox.com/s/q83q3mnvos8wdxk/reactome2ensg.rds?dl=1",
            "https://www.dropbox.com/s/9ah6aw0fborrp0z/string2ensg.rds?dl=1",
            "https://www.dropbox.com/s/7gj36rdec6spk9u/stitch2ensg.rds?dl=1",
            "https://www.dropbox.com/s/7f7370ae4k7b6hc/ensg2eg.rds?dl=1",
            "https://www.dropbox.com/s/3y9z8liv54jnhai/eg2ensg.rds?dl=1")

  names(urls) <- c("eg2rsids_10kb.rds",
                   "ensg2rsids_10kb.rds",
                   "ensp2rsids_10kb.rds",
                   "ensg2sym.rds",
                   "go.rds",
                   "reactome.rds",
                   "string.rds",
                   "stitch.rds",
                   "reactome2ensg.rds",
                   "string2ensg.rds",
                   "stitch2ensg.rds",
                   "ensg2eg.rds",
                   "eg2ensg.rds")
  for (feature in names(urls)) {
   message(paste("Downloading file:",feature))
   #destfile <- paste0(GAlist$dirs["gsets"],feature)
   destfile <- file.path(GAlist$dirs["gsets"],feature)
   download.file(url=urls[feature], mode = "wb", dest=destfile)
  }
  GAlist$gsetsfiles <- file.path(GAlist$dirs["gsets"],names(urls))
  #GAlist$gsetsfiles <- paste0(GAlist$dirs["gsets"],names(urls))
  names(GAlist$gsetsfiles) <- gsub(".rds","",names(urls))

  GAlist$gsets <- vector(mode = "list", length = length(GAlist$gsetsfiles))
  for(i in 1:length(GAlist$gsetsfiles)) {
   GAlist$gsets[[i]] <- readRDS(GAlist$gsetsfiles[i])
  }
  GAlist$gsets[[1]] <- qgg:::mapSets(sets=GAlist$gsets[[1]], rsids=GAlist$rsids, index=FALSE)
  GAlist$gsets[[2]] <- qgg:::mapSets(sets=GAlist$gsets[[2]], rsids=GAlist$rsids, index=FALSE)
  GAlist$gsets[[3]] <- qgg:::mapSets(sets=GAlist$gsets[[3]], rsids=GAlist$rsids, index=FALSE)
  names(GAlist$gsets) <- gsub(".rds","",names(urls))

 }

 if(what=="gsea") {
  # download gsea files in the database
  message("Downloading gsea results")
  urls <- c("https://www.dropbox.com/s/bm8iq8jw0rtlge7/gseaGenes.rds?dl=1",
            "https://www.dropbox.com/s/spe569pwkct1w61/gseaGO.rds?dl=1",
            "https://www.dropbox.com/s/t9bd3ulc8wlo53t/gseaPathways.rds?dl=1",
            "https://www.dropbox.com/s/2qyh12h0kaph2z3/gseaProteinComplexes.rds?dl=1",
            "https://www.dropbox.com/s/lhtf2mk41nhwxa0/gseaChemicalComplexes.rds?dl=1")


  urls <- c("https://www.dropbox.com/s/uxm4wz9l5jls3pf/ct_gseaChromosomes_gdtdb.rds?dl=1",
  "https://www.dropbox.com/s/bocpdcb3whqp60e/ct_gseaGenes_gdtdb.rds?dl=1",
  "https://www.dropbox.com/s/b35i6h4k9vrx8rf/ct_gseaGO_gdtdb.rds?dl=1",
  "https://www.dropbox.com/s/vtcs5xufpxdache/ct_gseaPathways_gdtdb.rds?dl=1",
  "https://www.dropbox.com/s/fa73q57wcfldk5o/ct_gseaProteinComplexes_gdtdb.rds?dl=1",
  "https://www.dropbox.com/s/wohvb1hfiuxzlp3/ct_gseaChemicalComplexes_gdtdb.rds?dl=1")

  urls <- c("https://www.dropbox.com/s/ia5wmwrwiatwvqa/ct_gseaChromosomes_gdtdb.rds?dl=1",
  "https://www.dropbox.com/s/bocpdcb3whqp60e/ct_gseaGenes_gdtdb.rds?dl=1",
  "https://www.dropbox.com/s/b35i6h4k9vrx8rf/ct_gseaGO_gdtdb.rds?dl=1",
  "https://www.dropbox.com/s/vtcs5xufpxdache/ct_gseaPathways_gdtdb.rds?dl=1",
  "https://www.dropbox.com/s/fa73q57wcfldk5o/ct_gseaProteinComplexes_gdtdb.rds?dl=1",
  "https://www.dropbox.com/s/3puru40dph8zu33/ct_gseaChemicalComplexes_gdtdb.rds?dl=1")

  names(urls) <- c("ct_gseaChromosomes_gdtdb.rds",
                   "ct_gseaGenes_gdtdb.rds",
                   "ct_gseaGO_gdtdb.rds",
                   "ct_gseaPathways_gdtdb.rds",
                   "ct_gseaProteinComplexes_gdtdb.rds",
                   "ct_gseaChemicalComplexes_gdtdb.rds")

  for (feature in names(urls)) {
   message(paste("Downloading file:",feature))
   destfile <- file.path(GAlist$dirs["gsea"],feature)
   download.file(url=urls[feature], mode = "wb", dest=destfile)
  }
  GAlist$gseafiles <- file.path(GAlist$dirs["gsea"],names(urls))
  names(GAlist$gseafiles) <- gsub(".rds","",names(urls))
  #names(GAlist$gseafiles) <- gsub("gsea","",names(GAlist$gseafiles))

 }
 if(what=="gstat") {
  # download gstat files in the database
  message("Downloading GWAS summary statistics")
  url_stat <- "https://www.dropbox.com/s/qrcivih31iuuril/stat.rds?dl=1"
  destfile <- file.path(GAlist$dirs["gstat"],"gstat.rds")
  download.file(url=url_stat, mode = "wb", dest=destfile)
  GAlist$gstatfiles <- file.path(GAlist$dirs["gstat"],"gstat.rds")

  url_stat <- "https://www.dropbox.com/s/0sizkeuw0sl51tn/GWAS_information.csv?dl=1"
  destfile <- file.path(GAlist$dirs["gstat"],"GWAS_information.csv")
  download.file(url=url_stat, mode = "wb", dest=destfile)
  GAlist$study <- as.list(read.csv2(destfile))

  url_stat <- c("https://www.dropbox.com/s/iqd7c4bbds03nrg/GWAS1.txt.gz?dl=1",
                "https://www.dropbox.com/s/pdv1fg280n86dwg/GWAS2.txt.gz?dl=1",
                "https://www.dropbox.com/s/t3y05ex1uouo4qg/GWAS3.txt.gz?dl=1",
                "https://www.dropbox.com/s/34yd6sxltqiv6e8/GWAS4.txt.gz?dl=1",
                "https://www.dropbox.com/s/xyfadehraaajkol/GWAS5.txt.gz?dl=1",
                "https://www.dropbox.com/s/2b5u6m60l6a5e82/GWAS6.txt.gz?dl=1",
                "https://www.dropbox.com/s/zilng15j7c0kl8n/GWAS7.txt.gz?dl=1",
                "https://www.dropbox.com/s/oo5o7suu2bx04bj/GWAS8.txt.gz?dl=1")

  for(study in 1:length(url_stat)) {
   destfile <- file.path(GAlist$dirs["gstat"],substring(url_stat[study],43,nchar(url_stat[study])-5))
   download.file(url=url_stat[study], mode = "wb", dest=destfile)
  }
  GAlist$studyfiles <- file.path(GAlist$dirs["gstat"],paste0(GAlist$study$file,".gz"))
  names(GAlist$studyfiles) <- GAlist$study$id

 }
 if(what=="marker") {
  # download marker files in the database
  message("Downloading marker information")
  url_stat <- "https://www.dropbox.com/s/4k54owkby3uf2hf/markers.txt.gz?dl=1"
  destfile <- file.path(GAlist$dirs["marker"],"markers.txt.gz")
  download.file(url=url_stat, mode = "wb", dest=destfile)
  GAlist$markerfiles <-file.path(GAlist$dirs["marker"],"markers.txt.gz")
  GAlist$markers <- fread(GAlist$markerfiles, data.table=FALSE)
  GAlist$rsids <- GAlist$markers$rsids
  GAlist$cpra <- paste(GAlist$markers$chr,
                       GAlist$markers$pos,
                       GAlist$markers$ea,
                       GAlist$markers$nea,sep="_")
 }
 if(what=="ensembl") {
  url_ensembl <- "https://ftp.ensembl.org/pub/release-109/tsv/homo_sapiens/Homo_sapiens.GRCh38.109.entrez.tsv.gz"
  destfile <- file.path(GAlist$dirs["gsets"],"Homo_sapiens.GRCh38.109.entrez.tsv.gz")
  download.file(url=url_ensembl, mode = "wb", dest=destfile)
  file_ensembl <-file.path(GAlist$dirs["gsets"],"Homo_sapiens.GRCh38.109.entrez.tsv.gz")
  ensembl <- fread(file_ensembl, data.table=FALSE)
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

 }

 if(what=="pubmed") {
  url_pubmed <- "https://ftp.ncbi.nih.gov/gene/DATA/gene2pubmed.gz"
  destfile <- file.path(GAlist$dirs["gsets"],"gene2pubmed.gz")
  download.file(url=url_pubmed, mode = "wb", dest=destfile)
  pubmed <- fread(destfile, data.table = FALSE)
  pubmed <- pubmed[pubmed[,1]%in%9606,-1]
  eg2pmid <- split(pubmed$PubMed_ID,f=as.factor(pubmed$GeneID))
  pmid2eg <- split(pubmed$GeneID,f=as.factor(pubmed$PubMed_ID))
  saveRDS(eg2pmid,file=file.path(GAlist$dirs["gsets"],"eg2pmid.rds"))
  saveRDS(pmid2eg,file=file.path(GAlist$dirs["gsets"],"pmid2eg.rds"))
 }
 if(what=="pubchem") {
  #https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/xrefs/PubMedID/TXT
  #https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/
 }

 if(what=="string") {
  url_string <- "https://stringdb-static.org/download/protein.links.v11.5/9606.protein.links.v11.5.txt.gz"
  destfile <- file.path(GAlist$dirs["gsets"],"9606.protein.links.v11.5.txt.gz")
  download.file(url=url_string, mode = "wb", dest=destfile)
  string <- fread(destfile, data.table=FALSE)
  string$protein1 <- gsub("9606.","",string$protein1)
  string$protein2 <- gsub("9606.","",string$protein2)
  fwrite(string, file=destfile)

  #Create sets
  string  <- string[string$combined_score>=min_combined_score,]
  string <- string[string$protein2%in%names(GAlist$gsets$ensp2ensg),]
  string <- split( string$protein2,f=as.factor(string$protein1))
  sets <- string[sapply(string ,length)>=min_interactions]
  sets <- lapply(sets,function(x){na.omit(unlist(GAlist$gsets$ensp2ensg[x]))})
  GAlist$gsets[[7]] <- sets
 }

 if(what=="stitch") {
  file_stitch <- "http://stitch.embl.de/download/protein_chemical.links.v5.0/9606.protein_chemical.links.v5.0.tsv.gz"
  stitch <- fread(file_stitch, data.table=FALSE)
  stitch$protein <- gsub("9606.","",stitch$protein)
  stitch <- stitch[stitch$protein%in%names(GAlist$gsets$ensp2ensg),]
  stitch  <- stitch[stitch$combined_score>=min_combined_score,]
  stitch <- split( stitch$protein,f=as.factor(stitch$chemical))
  sets  <- stitch[sapply(stitch ,length)>=min_interactions]
  sets <- lapply(sets,function(x){na.omit(unlist(GAlist$gsets$ensp2ensg[x]))})
  GAlist$gsets[[8]] <- sets
 }

 if(what=="1000G") {
  url <- "https://www.dropbox.com/s/jk3p47jf8ser6se/1000G_EUR_Phase3_plink.zip?dl=1"
  dest <- file.path(GAlist$dirs["marker"],"1000G_EUR_Phase3_plink.zip")
  download.file(url=url, mode = "wb", dest=dest)
  unzip(dest, exdir=GAlist$dirs["marker"])
 }

 if(what=="reactome") {
  url_db <- "https://reactome.org/download/current/ReactomePathways.txt"
  destfile <- file.path(GAlist$dirs["gsets"],"ReactomePathways.txt")
  download.file(url=url_db, mode = "wb", dest=destfile)

  url_db <- "https://reactome.org/download/current/Ensembl2Reactome.txt"
  destfile <- file.path(GAlist$dirs["gsets"],"Ensembl2Reactome.txt")
  download.file(url=url_db, mode = "wb", dest=destfile)

  pathway <- fread(file.path(GAlist$dirs["gsets"],"ReactomePathways.txt"), data.table=FALSE, header=FALSE)
  isHSA <- grep("R-HSA",pathway[,1])
  pathway <- pathway[isHSA,]
  path2names <- pathway[,2]
  names(path2names) <- pathway[,1]

  GAlist$gsets$path2names <- path2names


  reactome <- fread(file.path(GAlist$dirs["gsets"],"Ensembl2Reactome.txt"),
                    data.table=FALSE, header=FALSE)
  isHSA <- grep("R-HSA",reactome[,2])
  reactome <- reactome[isHSA,]
  head(reactome)

  reac2ens <- split(reactome[,1],f=reactome[,2])
  ens2reac <- split(reactome[,2],f=reactome[,1])

  GAlist$gsets$reac2ens <- reac2ens
  GAlist$gsets$ens2reac <- ens2reac
 }


 if(what=="dgi") {
  # download dgidb files in the database
  message("Downloading Drug Gene Interaction database")
  url_db <- "https://www.dgidb.org/data/monthly_tsvs/2022-Feb/interactions.tsv"
  destfile <- file.path(GAlist$dirs["drugdb"],"interactions.tsv")
  download.file(url=url_db, mode = "wb", dest=destfile)
  url_db <- "https://www.dgidb.org/data/monthly_tsvs/2022-Feb/genes.tsv"
  destfile <- file.path(GAlist$dirs["drugdb"],"genes.tsv")
  download.file(url=url_db, mode = "wb", dest=destfile)
  url_db <- "https://www.dgidb.org/data/monthly_tsvs/2022-Feb/drugs.tsv"
  destfile <- file.path(GAlist$dirs["drugdb"],"drugs.tsv")
  download.file(url=url_db, mode = "wb", dest=destfile)
  url_db <- "https://www.dgidb.org/data/monthly_tsvs/2022-Feb/categories.tsv"
  destfile <- file.path(GAlist$dirs["drugdb"],"categories.tsv")
  download.file(url=url_db, mode = "wb", dest=destfile)

  drugdb <- fread(file.path(GAlist$dirs["drugdb"], "interactions.tsv"),
                  quote = "", data.table = FALSE)

  drug2eg <- split( drugdb$entrez_id, f=as.factor(drugdb$drug_name) )
  drug2eg <- lapply(drug2eg,function(x){as.character(x)})
  drug2eg <- drug2eg[!names(drug2eg)==""]

  drugGenes <- lapply(drug2eg,function(x){na.omit(unlist(GAlist$gsets$eg2ensg[x]))})
  drugGenes <- lapply(drugGenes, function(x){unique(x)})
  drugGenes <- drugGenes[sapply(drugGenes, function(x){ !any(is.na(x)) } )]
  drugGenes <- drugGenes[ sapply(drugGenes, length)>0]
  saveRDS(drugGenes,file=file.path(GAlist$dirs["gsets"],"drugGenes.rds"))

  string <- getSetsDB(GAlist=GAlist, feature="String")
  drug2ensp <- lapply(drugGenes,function(x){na.omit(unlist(GAlist$gsets$ensg2ensp[x]))})
  drugComplex <- lapply(drug2ensp,function(x){na.omit(unlist(string[x]))})
  drugComplex <- lapply(drugComplex,function(x){na.omit(unlist(GAlist$gsets$ensp2ensg[x]))})
  drugComplex <- lapply(drugComplex, function(x){unique(x)})

  for(i in 1:length(drugComplex)) {
   drugComplex[[i]] <- unique(c(drugGenes[[i]], drugComplex[[i]]))
  }
  saveRDS(drugComplex,file=file.path(GAlist$dirs["gsets"],"drugComplex.rds"))

  drugGenesSets <- lapply(drugGenes,function(x){unique(unlist(GAlist$gsets$ensg2rsids_10kb[x]))})
  drugGenesSets <- drugGenesSets[!sapply(drugGenesSets,is.null)]
  saveRDS(drugGenesSets,file=file.path(GAlist$dirs["gsets"],"drugGenesSets.rds"))

  drugComplexSets <- lapply(drugComplex,function(x){unique(unlist(GAlist$gsets$ensg2rsids_10kb[x]))})
  drugComplexSets <- drugComplexSets[!sapply(drugComplexSets,is.null)]
  saveRDS(drugComplexSets,file=file.path(GAlist$dirs["gsets"],"drugComplexSets.rds"))

 }

 if(what=="pharmgkb") {
  cwd <- getwd()
  setwd(GAlist$dirs["drugdb"])
  url <- "https://api.pharmgkb.org/v1/download/file/data/drugLabels.zip"
  output_file <-  file.path(GAlist$dirs["drugdb"],"drugLabels.zip")
  download.file(url, destfile = output_file, mode = "wb")
  unzip(output_file)
  file.remove(output_file)

  url <- "https://api.pharmgkb.org/v1/download/file/data/relationships.zip"
  output_file <-  file.path(GAlist$dirs["drugdb"],"relationships.zip")
  download.file(url, destfile = output_file, mode = "wb")
  unzip(output_file)
  file.remove(output_file)

  url <- "https://api.pharmgkb.org/v1/download/file/data/clinicalVariants.zip"
  output_file <-  file.path(GAlist$dirs["drugdb"],"clinicalVariants.zip")
  download.file(url, destfile = output_file, mode = "wb")
  unzip(output_file)
  file.remove(output_file)

  url <- "https://api.pharmgkb.org/v1/download/file/data/automated_annotations.zip"
  output_file <-  file.path(GAlist$dirs["drugdb"],"automated_annotations.zip")
  download.file(url, destfile = output_file, mode = "wb")
  unzip(output_file)
  file.remove(output_file)
  setwd(cwd)

 }

 if(what=="opentargets") {
  dir.create(file.path(GAlist$dirs["drugdb"], "targets"))
  dir.create(file.path(GAlist$dirs["drugdb"], "diseases"))
  dir.create(file.path(GAlist$dirs["drugdb"], "diseaseToPhenotype"))
  dir.create(file.path(GAlist$dirs["drugdb"], "significantAdverseDrugReactions"))
  dir.create(file.path(GAlist$dirs["drugdb"], "associationByDatasourceDirect"))
  dir.create(file.path(GAlist$dirs["drugdb"], "associationByDatatypeDirect"))
  dir.create(file.path(GAlist$dirs["drugdb"], "associationByOverallDirect"))

  # Download opentarget
  urls <- c("ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.11/output/etl/json/targets/",
            "ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.11/output/etl/json/diseases/",
            "ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.11/output/etl/json/diseaseToPhenotype/",
            "ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.11/output/etl/json/fda/significantAdverseDrugReactions/",
            "ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.11/output/etl/json/associationByDatasourceDirect/",
            "ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.11/output/etl/json/associationByDatatypeDirect/",
            "ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.11/output/etl/json/associationByOverallDirect/")

  for (i in 1:length(urls)) {
   urldir <- unlist(strsplit(urls[i],split="/"))
   urldir <- urldir[length(urldir)]
   files = getURL(urls[i], ftp.use.epsv = FALSE, dirlistonly = TRUE)
   files <- stringr::str_c(urls[i], stringr::str_split(files, "\n")[[1]])
   files <- stringr::str_trim(files)
   is.json <- grep(".json",files, fixed=TRUE)
   files <- files[is.json]
   for (j in 1:length(files)) {
    destfile <- unlist(strsplit(files[j],split="/"))
    destfile <- destfile[length(destfile)]
    destfile <- file.path(file.path(GAlist$dirs["drugdb"], urldir),destfile)
    download.file(url=files[j], mode = "wb", dest=destfile)
   }
  }
  files <- dir(file.path(GAlist$dirs["drugdb"], "associationByOverallDirect"), full.names = TRUE)
  df <- NULL
  for (i in 1:length(files)) {
   con <- file(files[i],"r")
   df <- rbind(df,jsonlite::stream_in(con))
   close(con)
  }
  filename <- file.path(GAlist$dirs["drugdb"], "associationByOverallDirect.tsv")
  fwrite(df, file=filename)

  files <- dir(file.path(GAlist$dirs["drugdb"], "associationByDatatypeDirect"), full.names = TRUE)
  df <- NULL
  for (i in 1:length(files)) {
   con <- file(files[i],"r")
   df <- rbind(df,jsonlite::stream_in(con))
   close(con)
  }
  filename <- file.path(GAlist$dirs["drugdb"], "associationByDatatypeDirect.tsv")
  fwrite(df, file=filename)

  files <- dir(file.path(GAlist$dirs["drugdb"], "associationByDatasourceDirect"), full.names = TRUE)
  df <- NULL
  for (i in 1:length(files)) {
   con <- file(files[i],"r")
   df <- rbind(df,jsonlite::stream_in(con))
   close(con)
  }
  filename <- file.path(GAlist$dirs["drugdb"], "associationByDatasourceDirect.tsv")
  fwrite(df, file=filename)
 }

 if(what=="gwascatalog") {
  # http://ftp.ebi.ac.uk/pub/databases/gwas/releases/
  # timeout=600 required for large file
  options(download.file.method="libcurl", url.method="libcurl", timeout=600)
  dbdir <- file.path(GAlist$dbdir, "gwas")
  if(!dir.exists(dbdir)) dir.create(dbdir)

  file_studies <- "http://ftp.ebi.ac.uk/pub/databases/gwas/releases/2023/04/25/gwas-catalog-studies_ontology-annotated.tsv"
  destfile <- file.path(dbdir, "gwas-catalog-studies_ontology-annotated.tsv")
  download.file(file_studies, destfile = destfile, mode = "wb")

  file_gwas <- "http://ftp.ebi.ac.uk/pub/databases/gwas/releases/2023/04/25/gwas-catalog-associations_ontology-annotated.tsv"
  destfile <- file.path(dbdir, "gwas-catalog-associations_ontology-annotated.tsv")
  download.file(file_gwas, destfile = destfile, mode = "wb")

  gwas <- fread(destfile, data.table=FALSE, quote="")
 }

 if(what=="drugbank") {
  setwd(GAlist$dirs["drugdb"])
  url <- "https://www.dropbox.com/s/2rf7plojpkmy18z/drugid2drugname.rds?dl=1"
  file_drugbank <- file.path(GAlist$dirs["drugdb"],"drugid2drugname.rds")
  download.file(url, destfile = file_drugbank, mode = "wb")
  GAlist$drugbank <- readRDS(file=file_drugbank)
 }

 if(what=="atc") {
  # Add targets to GAlist
 df <- fread("https://www.dropbox.com/s/n5ehglmhhs0kcue/WHO%20ATC-DDD%202023-03-28.csv?dl=1",data.table = FALSE)
 GAlist$atc <- NULL
 GAlist$atc$code <- df$atc_code
 GAlist$atc$name <- df$atc_name

 # # Add drug target data frame with ATC information to GAlist
 drugGenes <- getSetsDB(GAlist = GAlist, feature = "DrugGenes")
 nreps <- sapply(drugGenes,length)
 drugs <- rep(names(drugGenes), times=nreps)
 ensg <- unlist(drugGenes)
 ensg2drug <- split(drugs, f=as.factor(ensg))
 df <- data.frame(Drug=drugs, Target=ensg)
 df$ATC <- rep("Unknown",nrow(df))
 has_atc <- match(tolower(df$Drug),tolower(GAlist$atc$name))
 df$ATC[!is.na(has_atc)] <- as.character(GAlist$atc$code[has_atc[!is.na(has_atc)]])
 GAlist$targets <- df
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

 return(GAlist)
}



#' @export
#'
annotationDB <- function(GAlist=NULL,
                         kb=10,
                         string=NULL,
                         stitch=NULL,
                         string_min_interaction=5,
                         stitch_min_interaction=5) {
 file_string <- string
 file_stitch <- stitch

 # http://bioconductor.org/packages/release/bioc/html/AnnotationHub.html
 #install.packages("BiocManager")
 #BiocManager::install("org.Hs.eg.db")
 #BiocManager::install("reactome.db")
 #BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

 #library(data.table)
 #library(org.Hs.eg.db)
 #library(reactome.db)

 if(is.null(GAlist)) {
  message("Please provide GAlist used for creating marker sets")
 }

 # Input annotation files from string and stitch (more can be added)
 if(is.null(file_string)) {
  file_string <- "https://stringdb-static.org/download/protein.links.v11.5/9606.protein.links.v11.5.txt.gz"
 }
 if(is.null(file_stitch)) {
  file_stitch <- "http://stitch.embl.de/download/protein_chemical.links.v5.0/9606.protein_chemical.links.v5.0.tsv.gz"
 }
 message(paste("STRING file used:",file_string))
 message(paste("STITCH file used:",file_string))

 ###########################################
 # Map rsids to chromosomal location
 ###########################################

 #rsids2chrpos <- data.frame(rsids=unlist(Glist$rsids),
 #                           chr=unlist(Glist$chr),
 #                           pos=unlist(Glist$pos),
 #                           stringsAsFactors=FALSE)
 rsids2chrpos <- GAlist$markers[,c("rsids","chr","pos")]

 ###########################################
 # Map ENSEMBL GENE IDS to ENTREZ IDs
 ###########################################
 ensg2eg <- as.list(org.Hs.egENSEMBL2EG)
 saveRDS(ensg2eg,file="ensg2eg.rds")

 ###########################################
 # Map ENSEMBL PROTEIN IDS to ENTREZ IDs
 ###########################################
 ensp2eg <- as.list(org.Hs.egENSEMBLPROT2EG)

 ###########################################
 # Map ENTREZ IDs to ENSEMBL GENE IDS
 ###########################################
 eg2ensg <- org.Hs.egENSEMBL
 mapped_genes <- mappedkeys(eg2ensg)
 eg2ensg <- as.list(eg2ensg[mapped_genes])
 saveRDS(eg2ensg,file="eg2ensg.rds")

 ###########################################
 # Map ENTREZ IDs to ENSEMBL PROTEIN IDS
 ###########################################
 eg2ensp <- org.Hs.egENSEMBLPROT
 mapped_genes <- mappedkeys(eg2ensp)
 eg2ensp <- as.list(eg2ensp[mapped_genes])

 ###########################################
 # Map ENTREZ IDs to CHROMOSOME
 ###########################################
 eg2chr <- org.Hs.egCHR
 mapped_genes <- mappedkeys(eg2chr)
 eg2chr <- as.list(eg2chr[mapped_genes])
 eg2chr <- sapply(eg2chr,function(x){x[[1]]})

 ###########################################
 # Map ENTREZ IDs to START POSITION
 ###########################################
 eg2chrloc <- org.Hs.egCHRLOC
 mapped_genes <- mappedkeys(eg2chrloc)
 eg2chrloc <- as.list(eg2chrloc[mapped_genes])
 eg2chrloc <- sapply(eg2chrloc,function(x){x[[1]]})
 eg2chrloc <- abs(eg2chrloc)

 ###########################################
 # Map ENTREZ IDs to END POSITION
 ###########################################
 eg2chrlocend <- org.Hs.egCHRLOCEND
 mapped_genes <- mappedkeys(eg2chrlocend)
 eg2chrlocend <- as.list(eg2chrlocend[mapped_genes])
 eg2chrlocend <- sapply(eg2chrlocend,function(x){x[[1]]})
 eg2chrlocend <- abs(eg2chrlocend)

 ###########################################
 # Map ENTREZ IDs to chromosomal position
 ###########################################
 eg2chrpos <- data.frame(chr=eg2chr[names(eg2chrloc)],
                         start=eg2chrloc,
                         stop=eg2chrlocend[names(eg2chrloc)],
                         stringsAsFactors=FALSE)
 eg2chrpos <- eg2chrpos[!eg2chrpos$chr=="X",]
 eg2chrpos <- eg2chrpos[!eg2chrpos$chr=="Y",]
 eg2chrpos <- eg2chrpos[!eg2chrpos$chr=="MT",]
 eg2chrpos <- eg2chrpos[!eg2chrpos$chr=="Un",]


 ###########################################
 # Map ENSEMBL GENE IDS to symbol
 ###########################################
 ensg2eg <- as.list(org.Hs.egENSEMBL2EG)
 eg2ensg <- org.Hs.egENSEMBL
 mapped_genes <- mappedkeys(eg2ensg)
 eg2ensg <- as.list(eg2ensg[mapped_genes])
 eg2sym <- org.Hs.egSYMBOL
 mapped_genes <- mappedkeys(eg2sym)
 eg2sym <- as.list(eg2sym[mapped_genes])

 ensg2sym <- sapply(ensg2eg, function(x){
  paste(unlist(eg2sym[x], use.names=FALSE),collapse=" ")})
 saveRDS(ensg2sym,file="ensg2sym.rds")

 ensp2ensg <- sapply(ensp2eg, function(x){
  unlist(eg2ensg[x],use.names=FALSE)})
 saveRDS(ensp2ensg,file="ensp2ensg.rds")

 ###########################################
 # Map ENTREZ IDs to rsids/cpra
 ###########################################
 eg2rsids <- vector(mode="list", length=nrow(eg2chrpos))
 eg2cpra <- vector(mode="list", length=nrow(eg2chrpos))
 names(eg2rsids) <- names(eg2cpra) <- rownames(eg2chrpos)

 for (i in 1:nrow(eg2chrpos)) {
  region <- rsids2chrpos$chr==eg2chrpos$chr[i] &  rsids2chrpos$pos>(eg2chrpos$start[i]-kb*1000) & rsids2chrpos$pos<(eg2chrpos$stop[i]+kb*1000)
  if(sum(region)>0) eg2rsids[[i]] <- rsids2chrpos$rsids[region]
  if(sum(region)>0) eg2cpra[[i]] <- rownames(rsids2chrpos)[region]
  message(paste("Processing gene:",i,"out of:",nrow(eg2chrpos)))
 }
 eg2rsids <- eg2rsids[!sapply(eg2rsids,is.null)]
 eg2cpra <- eg2cpra[!sapply(eg2cpra,is.null)]
 length(eg2rsids)
 length(eg2cpra)

 saveRDS(eg2rsids,file="eg2rsids.rds")
 saveRDS(eg2cpra,file="eg2cpra.rds")

 eg2rsids <- readRDS(file="eg2rsids.rds")
 eg2cpra <- readRDS(file="eg2cpra.rds")

 ###########################################
 # Map ENSEMBL GENE IDs to rsids/cpra
 ###########################################
 egs <- names(eg2rsids)
 ensg2rsids <- ensg2cpra <- vector(mode="list", length=length(ensg2eg))
 names(ensg2rsids) <- names(ensg2cpra) <- names(ensg2eg)
 for (i in 1:length(ensg2rsids)) {
  eg <- ensg2eg[[i]]
  eg <- eg[eg%in%egs]
  if(length(eg)>0) ensg2rsids[[i]] <-  unique(unlist(eg2rsids[eg]))
  if(length(eg)>0) ensg2cpra[[i]] <-  unique(unlist(eg2cpra[eg]))
  message(paste("Processing gene:",i,"out of:",length(ensg2rsids)))
 }
 ensg2rsids <- ensg2rsids[!sapply(ensg2rsids,is.null)]
 ensg2cpra <- ensg2cpra[!sapply(ensg2cpra,is.null)]
 length(ensg2rsids)
 length(ensg2cpra)

 saveRDS(ensg2rsids,file="ensg2rsids.rds")
 saveRDS(ensg2cpra,file="ensg2cpra.rds")

 ###########################################
 # Map ENSEMBL PROTEIN IDs to rsids/cpra
 ###########################################
 ensp2eg <- as.list(org.Hs.egENSEMBLPROT2EG)
 egs <- names(eg2rsids)
 ensp2rsids <- ensp2cpra <- vector(mode="list", length=length(ensp2eg))
 names(ensp2rsids) <- names(ensp2cpra) <- names(ensp2eg)
 for (i in 1:length(ensp2rsids)) {
  eg <- ensp2eg[[i]]
  eg <- eg[eg%in%egs]
  if(length(eg)>0) ensp2rsids[[i]] <-  unique(unlist(eg2rsids[eg]))
  if(length(eg)>0) ensp2cpra[[i]] <-  unique(unlist(eg2cpra[eg]))
  message(paste("Processing gene:",i,"out of:",length(ensp2rsids)))
 }
 ensp2rsids <- ensp2rsids[!sapply(ensp2rsids,is.null)]
 ensp2cpra <- ensp2cpra[!sapply(ensp2cpra,is.null)]
 length(ensp2rsids)
 length(ensp2cpra)

 saveRDS(ensp2rsids,file="ensp2rsids.rds")
 saveRDS(ensp2cpra,file="ensp2cpra.rds")


 ###########################################
 # Map GO IDs to rsids/cpra
 ###########################################
 go2eg <- as.list(org.Hs.egGO2EG)
 go2ensg <- lapply(go2eg,function(x){
  ensg <- na.omit(unlist(eg2ensg[x]))
  ensg[!duplicated(ensg)]
 })
 length(go2ensg)
 go2ensg <- go2ensg[!sapply(go2ensg,is.null)]
 length(go2ensg)
 ngenes <- sapply(go2ensg,length)
 go2ensg <- go2ensg[ngenes>=10]
 length(go2ensg)

 go2rsids <- go2cpra <- vector(mode="list", length=length(go2ensg))
 names(go2rsids) <- names(go2cpra) <- names(go2ensg)
 for (i in 1:length(go2ensg)) {
  ensg <- na.omit(unlist(go2ensg[i]))
  rsids <- unlist(ensg2rsids[ensg])
  go2rsids[[i]] <- rsids[!duplicated(rsids)]
  cpra <- unlist(ensg2cpra[ensg])
  go2cpra[[i]] <- cpra[!duplicated(cpra)]
 }
 go2rsids <- go2rsids[!sapply(go2rsids,is.null)]
 go2cpra <- go2cpra[!sapply(go2cpra,is.null)]
 length(go2rsids)
 length(go2cpra)

 saveRDS(go2ensg,file="go.rds")
 saveRDS(go2rsids,file="go2rsids.rds")
 saveRDS(go2cpra,file="go2cpra.rds")


 ######################################################
 # Map rsids/cpra to protein complexes from string
 ######################################################

 ensp <- unique(names(ensp2rsids))

 string <- fread(file_string, data.table=FALSE)
 string$protein1 <- gsub("9606.","",string$protein1)
 string$protein2 <- gsub("9606.","",string$protein2)
 string <- string[string$protein1%in%ensp,]
 string <- string[string$protein2%in%ensp,]
 string <- split( string$protein2,f=as.factor(string$protein1))
 string  <- string[sapply(string ,length)>=stitch_min_interactions]
 length(string)
 string2ensg <- sapply(string, function(x){unlist(ensp2ensg[x],use.names=FALSE)})
 string2ensg <- lapply(string2ensg,unique)


 string2rsids <- lapply(string,function(x){unique(unlist(ensp2rsids[x]))})
 string2cpra <- lapply(string,function(x){unique(unlist(ensp2cpra[x]))})
 string2rsids <- string2rsids[!sapply(string2rsids,is.null)]
 string2cpra <- string2cpra[!sapply(string2cpra,is.null)]
 length(string2rsids)
 length(string2cpra)

 saveRDS(string,file="string2ensp.rds")
 saveRDS(string2ensg,file="string2ensg.rds")
 saveRDS(string2rsids,file="string2rsids.rds")
 saveRDS(string2cpra,file="string2cpra.rds")

 ######################################################
 # Map rsids/cpra to chemical complexes from stitch
 ######################################################
 ensp <- unique(names(ensp2rsids))

 stitch <- fread(file_stitch, data.table=FALSE)
 stitch$protein <- gsub("9606.","",stitch$protein)
 stitch <- stitch[stitch$protein%in%ensp,]
 stitch <- split( stitch$protein,f=as.factor(stitch$chemical))
 length(stitch)
 stitch  <- stitch[sapply(stitch ,length)>=stitch_min_interaction]
 length(stitch)

 stitch2ensg <- sapply(stitch, function(x){unlist(ensp2ensg[x],use.names=FALSE)})
 stitch2ensg <- lapply(stitch2ensg,unique)

 stitch2rsids <- lapply(stitch,function(x){unique(unlist(ensp2rsids[x]))})
 stitch2cpra <- lapply(stitch,function(x){unique(unlist(ensp2cpra[x]))})
 stitch2rsids <- stitch2rsids[!sapply(stitch2rsids,is.null)]
 stitch2cpra <- stitch2cpra[!sapply(stitch2cpra,is.null)]
 length(stitch2rsids)
 length(stitch2cpra)

 saveRDS(stitch,file="stitch.rds")
 saveRDS(stitch2ensg,file="stitch2ensg.rds")
 saveRDS(stitch2rsids,file="stitch2rsids.rds")
 saveRDS(stitch2cpra,file="stitch2cpra.rds")


 ######################################################
 # Map rsids/cpra to biological pathways from reactome
 ######################################################

 reactome <- as.list(reactomePATHID2EXTID)
 reactome2ensg <- sapply(reactome, function(x){unlist(eg2ensg[x],use.names=FALSE)})
 length(reactome2ensg)
 reactome2ensg <- reactome2ensg[!sapply(reactome2engs,is.null)]
 length(reactome2ensg)
 reactome2ensg <- lapply(reactome2ensg,unique)

 reactome2rsids <- lapply(reactome,function(x){unique(unlist(eg2rsids[x]))})
 reactome2cpra <- lapply(reactome,function(x){unique(unlist(eg2cpra[x]))})

 reactome2rsids <- reactome2rsids[!sapply(reactome2rsids,is.null)]
 reactome2cpra <- reactome2cpra[!sapply(reactome2cpra,is.null)]
 length(reactome2rsids)
 length(reactome2cpra)
 str(reactome2rsids)

 saveRDS(reactome,file="reactome.rds")
 saveRDS(reactome2ensg,file="reactome2ensg.rds")
 saveRDS(reactome2rsids,file="reactome2rsids.rds")
 saveRDS(reactome2cpra,file="reactome2cpra.rds")

 eg2rsids <- readRDS(file="eg2rsids.rds")
 eg2rsids2indx <- qgg:::mapSets(sets=eg2rsids,rsids=rsids)
 saveRDS(eg2rsids2indx, file="indx_eg2rsids.rds")

 ensg2rsids <- readRDS(file="ensg2rsids.rds")
 ensg2rsids2indx <- qgg:::mapSets(sets=ensg2rsids,rsids=rsids)
 saveRDS(ensg2rsids2indx, file="indx_ensg2rsids.rds")

 ensp2rsids <- readRDS(file="ensp2rsids.rds")
 ensp2rsids2indx <- qgg:::mapSets(sets=ensp2rsids,rsids=rsids)
 saveRDS(ensp2rsids2indx, file="indx_ensp2rsids.rds")

 go2rsids <- readRDS(file="go2rsids.rds")
 go2rsids2indx <- qgg:::mapSets(sets=go2rsids,rsids=rsids)
 saveRDS(go2rsids2indx, file="indx_go2rsids.rds")

 reactome2rsids <- readRDS(file="reactome2rsids.rds")
 reactome2rsids2indx <- qgg:::mapSets(sets=reactome2rsids,rsids=rsids)
 saveRDS(reactome2rsids2indx, file="indx_reactome2rsids.rds")

 string2rsids <- readRDS(file="string2rsids.rds")
 string2rsids2indx <- qgg:::mapSets(sets=string2rsids,rsids=rsids)
 saveRDS(string2rsids2indx, file="indx_string2rsids.rds")

 stitch2rsids <- readRDS(file="stitch2rsids.rds")
 stitch2rsids2indx <- qgg:::mapSets(sets=stitch2rsids,rsids=rsids)
 saveRDS(stitch2rsids2indx, file="indx_stitch2rsids.rds")

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

 # internal summary statistic column format
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

 marker <- GAlist$markers
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

 if(any(isInteger)) {
  nlevs <- sapply(stat[,isInteger],function(x){length(unique(x))})
  isChr <- col_names[isInteger][nlevs<26]
  isPos <- col_names[isInteger][nlevs>23]
 }
 if(any(isCharacter)) {
  nlevs <- sapply(stat[,isCharacter],function(x){length(unique(x))})
  isAllele <- col_names[isCharacter][nlevs<5]
  if(length(isAllele)==1) isEA <- isNEA <- isAllele
  if(length(isAllele)==2) {
   isEA <- isAllele[1]
   isNEA <- isAllele[2]
  }
  isMarker <- col_names[isCharacter][nlevs>5]
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
  isEAF <- col_names[isNumeric][0<=minVal & maxVal<1 & medianVal>0.1 & medianVal<0.7]
  if(length(isEAF)==0) isEAF <- "Unknown"
  isInfo <- col_names[isNumeric][0.3<minVal & maxVal<=1 & medianVal>0.7]
  if(length(isInfo)==0) isInfo <- "Unknown"
 }
 best_guess <- c(isMarker,isChr,isPos, isEA, isNEA, isEAF,
                 isB, isSEB, isP, isN, isNca, isNco, isInfo)
 best_guess
}
