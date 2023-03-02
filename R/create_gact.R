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

  # Step 2: Download data from database:
  GAlist <- downloadDB(GAlist=GAlist, what="marker")
  GAlist <- downloadDB(GAlist=GAlist, what="gsets")
  GAlist <- downloadDB(GAlist=GAlist, what="gsea")
  GAlist <- downloadDB(GAlist=GAlist, what="gstat")
  GAlist <- downloadDB(GAlist=GAlist, what="dgidb")

  # Step 3: Create marker sets from database:
  if(what=="full") {
  message("Creating full marker sets - this may take some time")
  GAlist <- mapSetsDB(GAlist=GAlist)
  }
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
#' "gsets", "gsea", "ldsc", "gbayes", "marker", "raw", and "dgidb".
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
 GAlist$dbdir <- dbdir

 dirs <- c(glist = "glist",
           gstat = "gstat",
           gsets = "gsets",
           gsea = "gsea",
           ldsc = "ldsc",
           gbayes = "gbayes",
           marker = "marker",
           raw = "raw",
           dgidb = "dgidb",
           pharmgkb = "pharmgkb")

 lapply(names(dirs), function(x) {
  dir.create(file.path(dbdir, dirs[x]))
 })

 GAlist <- list(version = version,
                traits = NULL,
                dirs = file.path(dbdir, dirs),
                features = c("Markers", "Genes", "Proteins", "GO", "Pathways", "ProteinComplexes", "ChemicalComplexes"))

 names(GAlist$dirs) <- c("glist","gstat","gsets","gsea", "ldsc", "gbayes", "marker", "raw", "dgidb","pharmgkb")

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


# createDB <- function(Glist=NULL, version=NULL, dbdir=NULL, what="lite") {
#
#  if(is.null(version)) stop(paste("Please include a database name using the version argument"))
#  #if(is.null(Glist)) stop(paste("Please include a Glist using the Glist argument"))
#  dbdir <- paste0(dbdir,"/",version)
#  gstatdir <- paste0(dbdir,"/gstat/")
#  gsetsdir <- paste0(dbdir,"/gsets/")
#  gseadir <- paste0(dbdir,"/gsea/")
#  glistdir <- paste0(dbdir,"/glist/")
#  ldscdir <- paste0(dbdir,"/ldsc/")
#  gbayesdir <- paste0(dbdir,"/gbayes/")
#  markerdir <- paste0(dbdir,"/marker/")
#  rawdir <- paste0(dbdir,"/raw/")
#  dgidir <- paste0(dbdir,"/dgidb/")
#  if(dir.exists(dbdir)) stop(paste("Directory:",dbdir,"allready exists"))
#  if(!dir.exists(dbdir)) {
#   dir.create(dbdir)
#   dir.create(glistdir)
#   dir.create(gstatdir)
#   dir.create(gsetsdir)
#   dir.create(gseadir)
#   dir.create(ldscdir)
#   dir.create(gbayesdir)
#   dir.create(markerdir)
#   dir.create(rawdir)
#   dir.create(dgidir)
#  }
#
#  GAlist <- NULL
#  GAlist$version <- version
#
#  GAlist$traits <- NULL
#
#  GAlist$dirs <- c(glistdir,gstatdir,gsetsdir,gseadir, ldscdir, gbayesdir, markerdir, rawdir, dgidir)
#  names(GAlist$dirs) <- c("glist","gstat","gsets","gsea", "ldsc", "gbayes", "marker", "raw", "dgidb")
#
#  # features in the database
#  GAlist$features <- c("Markers","Genes","Proteins","GO","Pathways",
#                       "ProteinComplexes","ChemicalComplexes")
#
#  if(!is.null(Glist)) {
#   keep <- unlist(Glist$rsids)%in%unlist(Glist$rsidsLD)
#   GAlist$markers <- data.frame(rsids=unlist(Glist$rsids),
#                                chr=unlist(Glist$chr),
#                                pos=unlist(Glist$pos),
#                                ea=unlist(Glist$a1),
#                                nea=unlist(Glist$a2),
#                                eaf=unlist(Glist$af),
#                                stringsAsFactors=FALSE)[keep,]
#   GAlist$rsids <- unlist(Glist$rsids)[keep]
#   GAlist$cpra <- unlist(Glist$cpra)[keep]
#   file_markers <- paste0(GAlist$dirs["marker"],"markers.txt.gz")
#   fwrite(GAlist$markers, file=file_markers)
#  }
#
#  # update Glist$ldfiles
#  #ldfiles <- list.files(path=paste0(dbdir,"/glist/ldfiles"),pattern=".ld")
#  #rws <- sapply(ldfiles,function(x){grep(x,Glist$ldfiles)})
#  #rws <- order(rws)
#  #Glist$ldfiles <- paste0(dbdir,"/glist/",ldfiles[rws])
#
#  return(GAlist)
# }

#' @export
#'
mapSetsDB <- function(GAlist = NULL) {
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

# mapSetsDB <- function(GAlist=NULL) {
#  ensg2rsids <- GAlist$gsets[["ensg2rsids_10kb"]]
#
#  fset <- getSetsDB(GAlist=GAlist,feature="GO")
#  sets <- lapply(fset,function(x){unique(unlist(ensg2rsids[x]))})
#  sets <- sets[!sapply(sets,is.null)]
#  setsfile <- paste0(GAlist$dirs["gsets"],"go2rsids.rds")
#  saveRDS(sets,file=setsfile)
#
#  fset <- getSetsDB(GAlist=GAlist,feature="Pathways2Genes")
#  sets <- lapply(fset,function(x){unique(unlist(ensg2rsids[x]))})
#  sets <- sets[!sapply(sets,is.null)]
#  setsfile <- paste0(GAlist$dirs["gsets"],"reactome2rsids.rds")
#  saveRDS(sets,file=setsfile)
#
#  fset <- getSetsDB(GAlist=GAlist,feature="ProteinComplexes2Genes")
#  sets <- lapply(fset,function(x){unique(unlist(ensg2rsids[x]))})
#  sets <- sets[!sapply(sets,is.null)]
#  setsfile <- paste0(GAlist$dirs["gsets"],"string2rsids.rds")
#  saveRDS(sets,file=setsfile)
#
#  fset <- getSetsDB(GAlist=GAlist,feature="ChemicalComplexes2Genes")
#  sets <- lapply(fset,function(x){unique(unlist(ensg2rsids[x]))})
#  sets <- sets[!sapply(sets,is.null)]
#  setsfile <- paste0(GAlist$dirs["gsets"],"stitch2rsids.rds")
#  saveRDS(sets,file=setsfile)
#
#  GAlist$gsetsfiles[12] <- paste0(GAlist$dirs["gsets"],"go2rsids.rds")
#  GAlist$gsetsfiles[13] <- paste0(GAlist$dirs["gsets"],"reactome2rsids.rds")
#  GAlist$gsetsfiles[14] <- paste0(GAlist$dirs["gsets"],"string2rsids.rds")
#  GAlist$gsetsfiles[15] <- paste0(GAlist$dirs["gsets"],"stitch2rsids.rds")
#
#  names(GAlist$gsetsfiles[12]) <- "go2rsids"
#  names(GAlist$gsetsfiles[13]) <- "reactome2rsids"
#  names(GAlist$gsetsfiles[14]) <- "string2rsids"
#  names(GAlist$gsetsfiles[15]) <- "stitch2rsids"
#  return(GAlist)
# }

#' @export
#'
downloadDB <- function(GAlist=NULL, what=NULL) {

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

 if(what=="dgidb") {
  # download dgidb files in the database
  message("Downloading Drug Gene Interaction database")
  url_db <- "https://www.dgidb.org/data/monthly_tsvs/2022-Feb/interactions.tsv"
  destfile <- file.path(GAlist$dirs["dgidb"],"interactions.tsv")
  download.file(url=url_db, mode = "wb", dest=destfile)
  url_db <- "https://www.dgidb.org/data/monthly_tsvs/2022-Feb/genes.tsv"
  destfile <- file.path(GAlist$dirs["dgidb"],"genes.tsv")
  download.file(url=url_db, mode = "wb", dest=destfile)
  url_db <- "https://www.dgidb.org/data/monthly_tsvs/2022-Feb/drugs.tsv"
  destfile <- file.path(GAlist$dirs["dgidb"],"drugs.tsv")
  download.file(url=url_db, mode = "wb", dest=destfile)
  url_db <- "https://www.dgidb.org/data/monthly_tsvs/2022-Feb/categories.tsv"
  destfile <- file.path(GAlist$dirs["dgidb"],"categories.tsv")
  download.file(url=url_db, mode = "wb", dest=destfile)
 }
 return(GAlist)
}





#' Get summary statistics from the database
#'
#' This function retrieves summary statistics from the gact database of GSEA results,
#' such as marker statistics, gene statistics, protein statistics, etc.
#'
#' @param GAlist A list object providing information and infrastructure of the gact database.
#' @param feature The type of statistics to retrieve from the database, such as "Markers", "Genes", "Proteins", etc.
#' @param featureID A character vector of specific features to extract from the database. If not specified, all features will be returned.
#' @param file The file name of the database of GSEA results.
#' @param studyID The study ID to retrieve the results for. If not specified, results will be returned for all studies.
#' @param trait The trait to retrieve the results for. If not specified, results will be returned for all traits.
#' @param threshold A numeric p-value threshold used to filter results.
#' @param format The format to return the results in, either "list" (default) or "data.frame".
#' @return A list or data.frame of summary statistics for the specified feature, study ID, and trait.
#' @examples
#'
#' \dontrun{
#' # load example data
#' GAlist <- readRDS(file="GAlist.rds")
#'
#' # get marker statistics
#' stat <- getStatDB(GAlist = GAlist, feature = "Markers", studyID="GWAS1")
#'
#' # get gene statistics
#' stat <- getStatDB(GAlist = GAlist, feature = "Genes")
#'
#' # get gene statistics for a specific feature
#' stat <- getStatDB(GAlist = GAlist, feature = "Genes", featureID = "TP53")
#'
#' # get gene statistics filtered by a threshold
#' stat <- getStatDB(GAlist = GAlist, feature = "Genes", threshold = 0.05)
#' }


#' @export
#'
getStatDB <- function(GAlist=NULL, feature=NULL, featureID=NULL,file=NULL,
                    studyID=NULL, trait=NULL, threshold=0.95,
                    format="list") {
 features <- c("Markers","Genes","Proteins","GO","Pathways",
               "ProteinComplexes","ChemicalComplexes","Chromosomes")
 header <- c("Marker ID","Gene ID","Protein ID","GO ID","Pathway ID",
             "Protein ID","Chemical ID", "Chromosome ID")
 names(header) <- features
 if(!feature%in%features) stop(paste("feature:",feature,"not in database"))

 #gseafile <- paste0(GAlist$dirs["gsea"],"ct_gsea",feature,"_gdtdb.rds")
 gseafile <- GAlist$gseafiles[paste0("ct_gsea",feature,"_gdtdb")]
 res <- readRDS(gseafile)
 message(paste("Extract statistics based p-value threshold:",threshold))
 cls5 <- grep("_0.05", colnames(res$stat))
 cls95 <- grep("_0.95", colnames(res$stat))
 cls <- 1:ncol(res$stat)
 if(threshold==0.05) cls <- cls5
 if(threshold==0.95) cls <- cls95

 colnames(res$stat) <- gsub("z_","",colnames(res$stat))
 colnames(res$p) <- gsub("z_","",colnames(res$p))
 colnames(res$stat) <- gsub("_0.05","",colnames(res$stat))
 colnames(res$p) <- gsub("_0.05","",colnames(res$p))
 colnames(res$stat) <- gsub("_0.95","",colnames(res$stat))
 colnames(res$p) <- gsub("_0.95","",colnames(res$p))
 res$p[res$stat==0] <- 1
 rws <- rep(TRUE,lenth=nrow(res$stat))
 if(!is.null(featureID)) rws <- rownames(res$stat)%in%featureID
 if(sum(rws)==0) stop("None of featureIDs found in database")
 res$stat <- res$stat[rws,cls]
 res$p <- res$p[rws,cls]
 if(!is.null(studyID)) {
  res$stat <- res$stat[,colnames(res$stat)%in%studyID]
  res$p <- res$p[,colnames(res$p)%in%studyID]
 }

 if(format=="data.frame") {
  res <- as.data.frame(res)
  if(feature=="Genes") {
   res <- cbind(rownames(res),GAlist$gsets[["ensg2sym"]][rownames(res)], res)
   colnames(res)[1:2] <- c("Ensembl Gene ID","Symbol")
   ensg2sym_list <- lapply(res[,"Symbol"], function(x){
    unlist(strsplit(x,split=" "))})
   gsym <- unlist(ensg2sym_list)
   rws <- rep(1:nrow(res),times=sapply(ensg2sym_list,length))
   res <- res[rws,]
   res[,"Symbol"] <- gsym
   if(!is.null(featureID)) {
    select <- tolower(res[,"Ensembl Gene ID"])%in%tolower(featureID) | tolower(res[,"Symbol"])%in%tolower(featureID)
    if(any(select)) res <- res[select,]
    if(!any(select)) stop("None of featureIDs found in the database")
   }
  }
  if(!feature=="Genes") {
   res <- cbind(rownames(res), res)
   colnames(res)[1] <- header[feature]
  }

 }
 return(res)
}


#' @export
#'
writeStatDB <- function(GAlist=NULL, feature=NULL, featureID=NULL,
                      studyID=NULL, threshold=0.95,
                      format="data.frame", file.csv=NULL, hyperlink=TRUE) {
 res <- getStatDB(GAlist=GAlist, feature=feature, featureID=featureID,
                 studyID=studyID, threshold=threshold,
                 format=format)
 features <- c("Markers","Genes","Proteins","GO","Pathways",
               "ProteinComplexes","ChemicalComplexes","Chromosomes")
 if(!feature%in%features) stop(paste("feature:",feature,"not in database"))
 header <- c("Marker ID","Gene ID","Protein ID","GO ID","Pathway ID",
             "Protein ID","Chemical ID", "Chromosome ID")
 names(header) <- features
 if(feature=="Genes") {
  if(hyperlink) {
   res <- cbind(rownames(res), res)
   res2hyperlink_ensembl <- paste0("http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",res[,1])
   res2hyperlink_opentarget <- paste0("https://platform.opentargets.org/target/",res[,1])
   colnames(res)[1:3] <- c("Ensembl Gene ID","Open Target","Symbol")

   res2hyperlink_ensembl <- paste0("=Hyperlink(",'"',res2hyperlink_ensembl,'"',";",'"',res[,1],'"',")")
   res2hyperlink_opentarget <- paste0("=Hyperlink(",'"',res2hyperlink_opentarget,'"',";",'"',res[,1],'"',")")
   res[,1] <- res2hyperlink_ensembl
   res[,2] <- res2hyperlink_opentarget
  }
 }
 res2html <- c("https://www.ncbi.nlm.nih.gov/snp/",
               "http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",
               "http://www.ensembl.org/Homo_sapiens/Protein/Summary?p=",
               "https://www.ebi.ac.uk/QuickGO/term/",
               "https://reactome.org/content/detail/",
               "https://string-db.org/network/9606.",
               "https://pubchem.ncbi.nlm.nih.gov/compound/")

 names(res2html) <- c("Marker","Genes","Proteins","GO","Pathways",
                      "ProteinComplexes","ChemicalComplexes")

 if(!feature=="Genes") {
  res2hyperlink <- paste0(res2html[feature],res[,1])
  if(feature=="ChemicalComplexes") res2hyperlink <- paste0(res2html[feature],substring(res[,1],5,nchar(as.character(res[,1]))))
  res2hyperlink <- paste0("=Hyperlink(",'"',res2hyperlink,'"',";",'"',res[,1],'"',")")
  if(hyperlink) res[,1] <- res2hyperlink
  colnames(res)[1] <- header[feature]
 }

 write.csv2(res,file=file.csv,row.names=FALSE)
}

#' @export
#'
getStat <- function(GAlist=NULL, feature=NULL, featureID=NULL,
                    studyID=NULL, trait="t2d", threshold=1,
                    format="data.frame", hyperlink=FALSE, cls=NULL) {

 features <- c("Markers","Genes","Proteins","GO","Pathways",
               "ProteinComplexes","ChemicalComplexes")
 header <- c("Marker ID","Gene ID","Protein ID","GO ID","Pathway ID",
             "Protein ID","Chemical ID")
 names(header) <- features

 if(feature=="Markers") {
  res <- readRDS(GAlist$gstatfiles)
  res <- as.data.frame(res)
  rownames(res) <- res$rsids
  return(res)
 }

 if(!feature%in%GAlist$features) stop(paste("feature:",feature,"not in GACT database"))

 res <- readRDS(GAlist$gseafiles[feature])
 #cls <- c("z_0.001","z_0.05","z_0.95")
 if(!is.null(cls)) {
  #add check
 }
 if(is.null(cls)) cls <- c("p=0.01")
 colnames(res$stat) <- gsub("z_","p=",colnames(res$stat))
 colnames(res$p) <- gsub("z_","p=",colnames(res$p))
 res$p[res$stat==0] <- 1
 res$stat <- res$stat[,cls]
 res$p <- res$p[,cls]

 if(format=="data.frame") {
  res <- as.data.frame(res)

  if(feature=="Genes") {
   if(!hyperlink) {
    res <- cbind(rownames(res),GAlist$gsets[["ensg2sym"]][rownames(res)], res)
    colnames(res)[1:2] <- c("Ensembl Gene ID","Symbol")
   }
   if(hyperlink) {
    res <- cbind(rownames(res),rownames(res),GAlist$gsets[["ensg2sym"]][rownames(res)], res)
    res2hyperlink_ensembl <- paste0("http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",res[,1])
    res2hyperlink_opentarget <- paste0("https://platform.opentargets.org/target/",res[,1])
    colnames(res)[1:3] <- c("Ensembl Gene ID","Open Target","Symbol")

    res2hyperlink_ensembl <- paste0("=Hyperlink(",'"',res2hyperlink_ensembl,'"',";",'"',res[,1],'"',")")
    res2hyperlink_opentarget <- paste0("=Hyperlink(",'"',res2hyperlink_opentarget,'"',";",'"',res[,1],'"',")")
    res[,1] <- res2hyperlink_ensembl
    res[,2] <- res2hyperlink_opentarget
   }
  }

  res2html <- c("https://www.ncbi.nlm.nih.gov/snp/",
                "http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",
                "http://www.ensembl.org/Homo_sapiens/Protein/Summary?p=",
                "https://www.ebi.ac.uk/QuickGO/term/",
                "https://reactome.org/content/detail/",
                "https://string-db.org/network/9606.",
                "https://pubchem.ncbi.nlm.nih.gov/compound/")

  names(res2html) <- c("Marker","Genes","Proteins","GO","Pathways",
                       "ProteinComplexes","ChemicalComplexes")



  if(!feature=="Genes") {
   res <- cbind(rownames(res), res)
   res2hyperlink <- paste0(res2html[feature],res[,1])
   if(feature=="ChemicalComplexes") res2hyperlink <- paste0(res2html[feature],substring(res[,1],5,nchar(as.character(res[,1]))))
   res2hyperlink <- paste0("=Hyperlink(",'"',res2hyperlink,'"',";",'"',res[,1],'"',")")
   if(hyperlink) res[,1] <- res2hyperlink
   colnames(res)[1] <- header[feature]
  }
 }

 return(res)
}

#' #' @export
#' #'
#' writeStat <- function(GAlist=NULL, feature=NULL, featureID=NULL,
#'                       studyID=NULL, trait="T2D", threshold=1,
#'                       format="data.frame", file.csv=NULL, hyperlink=TRUE) {
#'  stat <- getStat(GAlist=GAlist, feature=feature, featureID=featureID,
#'                  studyID=studyID, trait=trait, threshold=threshold,
#'                  format=format, hyperlink=hyperlink)
#'  write.csv2(stat,file=file.csv,row.names=FALSE)
#' }

#' @export
#'
getStudiesDB <- function(GAlist=NULL) {
 return(as.data.frame(GAlist$study))
}

#' @export
#'
designSetsDB <- function(GAlist=NULL, feature=NULL, featureID=NULL, rowFeatureID=NULL) {
 if(is.null(GAlist)) stop ("Please provide GAlist")
 if(is.null(feature)) stop ("Please provide feature")
 sets <- getSetsDB(GAlist=GAlist, feature=feature)
 if(!is.null(featureID)) {
  select <- names(sets)%in%featureID
  if(sum(select)==0) stop("None of the fetureIDs found in sets")
  sets <- sets[select]
 }
 if(is.null(rowFeatureID)) rowFeatureID <- unique(unlist(sets))
 sets <- qgg:::mapSets(sets=sets,rsids=rowFeatureID, index=TRUE)
 W <- matrix(0,nrow=length(rowFeatureID), ncol=length(sets))
 colnames(W) <- names(sets)
 rownames(W) <- rowFeatureID
 for(i in 1:length(sets)) {
  W[sets[[i]],i] <- 1
 }
 return(W)
}

#' Retrieve sets of features from a GAlist object
#'
#' The `getSetsDB` function retrieves sets of features from a GAlist object. The
#' feature sets are specified by the `feature` argument, and the specific IDs of
#' the features can be filtered using the `featureID` argument.
#'
#' @param GAlist A list object providing information and infrastructure of the gact database.
#' @param feature A character string specifying the type of feature set to be
#'   retrieved. Possible values are "Entres Genes", "Genes", "Proteins", "Gene
#'   Symbol", "GO", "Pathways", "ProteinComplexes", "ChemicalComplexes",
#'   "ProteinComplexes2Genes", and "ChemicalComplexes2Genes".
#' @param featureID A character vector of IDs specifying the specific features to
#'   be retrieved.
#' @return A list of the specified feature sets. If `featureID` is not `NULL`,
#'   returns only the sets with IDs in `featureID`.
#'
#' @examples
#' \dontrun{
#' sets <- getSetsDB(GAlist, feature="Genes")
#' getSetsDB(GAlist, feature="Genes", featureID=c("gene1", "gene2"))
#' }
#'
#' @export
#'
getSetsDB <- function(GAlist=NULL, feature=NULL, featureID=NULL) {
 sets <- NULL
 if(feature=="Entres Genes") sets <- GAlist$gsets[[1]]
 if(feature=="Genes") sets <- GAlist$gsets[[2]]
 if(feature=="Proteins") sets <- GAlist$gsets[[3]]
 if(feature=="Gene Symbol") sets <- GAlist$gsets[[4]]
 if(feature=="GO") sets <- GAlist$gsets[[5]]
 if(feature=="Pathways") sets <- GAlist$gsets[[9]]
 if(feature=="Pathways2Genes") sets <- GAlist$gsets[[9]]
 if(feature=="ProteinComplexes") sets <- GAlist$gsets[[7]]
 if(feature=="ChemicalComplexes") sets <- GAlist$gsets[[8]]
 if(feature=="ProteinComplexes2Genes") sets <- GAlist$gsets[[10]]
 if(feature=="ChemicalComplexes2Genes") sets <- GAlist$gsets[[11]]
 if(!is.null(featureID)) {
  select <- names(sets)%in%featureID
  if(sum(select)==0) stop("None of the fetureIDs found in sets")
  #if(any(!select)) message(paste("Some IDs not in data base:",featureID[!select]))
  sets <- sets[select]
 }
 return(sets)
 }


#' Get Marker Sets from database
#'
#' This function retrieves marker sets based on a given feature and feature ID from a file. The feature can be one of the following: Genes, Inter Genes, Gene Symbol, Proteins, Pathways, GO, ProteinComplexes, ChemicalComplexes. The feature ID refers to the ID of the feature for which the marker sets are to be retrieved. The `rsids` argument is optional and refers to a list of SNP IDs for which the corresponding marker sets are to be retrieved.
#'
#' @param GAlist A list object providing information and infrastructure of the gact database.
#' @param feature A string specifying the type of biological feature for which the marker sets are to be retrieved. The options are: "Genes", "Entrez Genes", "Gene Symbol", "Proteins", "Pathways", "GO", "Protein Complexes", and "Chemical Complexes".
#' @param featureID The ID of the feature for which the marker sets are to be retrieved.
#' @param rsids A character vector of rsids to subset the marker sets.
#'
#' @return A list of marker sets for the specified feature and feature ID, where each set is a character vector of rsids.
#'
#' @examples
#' \dontrun{
#' sets <- getMarkerSetsDB(GAlist=GAlist, feature="Genes", featureID=c("ENSG00000165879", "ENSG00000012048"))
#' }
#' @export
#'
getMarkerSetsDB <- function(GAlist=NULL, feature=NULL, featureID=NULL, rsids=NULL) {

 if(feature=="Genes") setsfile <- file.path(GAlist$dirs["gsets"],"ensg2rsids_10kb.rds")
 if(feature=="Entres Genes") setsfile <- file.path(GAlist$dirs["gsets"],"eg2rsids_10kb.rds")
 if(feature=="Gene Symbol") setsfile <- file.path(GAlist$dirs["gsets"],"sym2rsids_10kb.rds")
 if(feature=="Proteins") setsfile <- file.path(GAlist$dirs["gsets"],"ensp2rsids_10kb.rds")
 if(feature=="Pathways") setsfile <- file.path(GAlist$dirs["gsets"],"reactome2rsids.rds")
 if(feature=="GO") setsfile <- file.path(GAlist$dirs["gsets"],"go2rsids.rds")
 if(feature=="ProteinComplexes") setsfile <- file.path(GAlist$dirs["gsets"],"string2rsids.rds")
 if(feature=="ChemicalComplexes") setsfile <- file.path(GAlist$dirs["gsets"],"stitch2rsids.rds")
 sets <- readRDS(file=setsfile)
 if(!is.null(featureID)) {
  inSet <- featureID%in%names(sets)
  if(any(!inSet)) warning(paste("Some IDs not in data base:",featureID[!inSet]))
  featureID <- featureID[featureID%in%names(sets)]
  sets <- sets[featureID]
 }
 if(!is.null(rsids)) sets <- qgg:::mapSets(sets=sets, rsids=rsids, index=FALSE)
 return(sets)
}


#' @export
#'
getFeatureDB <- function(GAlist=GAlist, feature=NULL, featureID=NULL, format="list") {
 if(feature=="Drug Gene Interactions") {
  eg2ensg <- readRDS(file.path(GAlist$dirs["gsets"],"eg2ensg.rds"))
  egs <- rep(names(eg2ensg),times=sapply(eg2ensg,length))
  ensg <- unlist(eg2ensg, use.names=FALSE)
  df1 <- data.frame(entrez_id=egs,ENSG=ensg)

  df2 <- fread(file.path(GAlist$dirs["dgidb"],"interactions.tsv"), data.table=FALSE)
  df2$entrez_id <- as.character(df2$entrez_id)

  df <- merge(x=df1, y=df2, by.x="entrez_id", by.y="entrez_id",  all.y=TRUE)
  colnames(df)[1:3] <- c("Entrez Gene","Ensembl Gene ID","Symbol")
 }
 return(df)
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
 GAlist$studyfiles <- file.path(GAlist$dirs["gstat"],GAlist$study$file,".gz")
 names(GAlist$studyfiles) <- GAlist$study$id


 # processing of summary statistics
 if(writeStatDB) {
  message("Perform quality control of external summary statistics")
  stat <- qcStatDB(GAlist=GAlist,stat=stat, excludeMAFDIFF=excludeMAFDIFF)
  message(paste("Writing processed summary statistics to internal file:",
                GAlist$study$file[study_number]))
  file_stat <- file.path(GAlist$dirs["gstat"],GAlist$study$file[study_number],".gz")
  if(file.exists(file_stat)) stop(paste("GWAS summary statistics file allready exists:",
                                        file_stat))
  fwrite(stat, file_stat)
 }

 return(GAlist)
}


#' @export
#'
getMarkerStatDB <- function(GAlist=NULL, studyID=NULL, what="all", format="list", rm.na=TRUE, rsids=NULL, cpra=NULL) {

 if(!is.null(cpra)) {
  cpra1 <- paste(GAlist$markers[,"chr"],
                 GAlist$markers[,"pos"],
                 toupper(GAlist$markers[,"ea"]),
                 toupper(GAlist$markers[,"nea"]), sep="_")
  cpra2 <- paste(GAlist$markers[,"chr"],
                 GAlist$markers[,"pos"],
                 toupper(GAlist$markers[,"nea"]),
                 toupper(GAlist$markers[,"ea"]),sep="_")

  mapped <- cpra1%in%cpra | cpra2%in%cpra
  message("Map markers based on cpra")
  message(paste("Number of markers in cpra mapped to marker ids in GAlist:",sum(mapped)))
  rsids <- GAlist$rsids[mapped]
 }

 if(is.null(studyID)) studyID <- GAlist$study$id
 names(GAlist$study$neff) <-GAlist$study$id

 if (length(studyID)==1) format <- "data.frame"
 if (length(studyID)>1) format <- "list"

 if(format=="list" && what=="all") {
  b <- seb <- z <- p <- n <- matrix(NA,ncol=length(studyID),nrow=length(GAlist$rsids))
  colnames(b) <- colnames(seb) <- colnames(z) <- colnames(p) <- colnames(n) <- studyID
  rownames(b) <- rownames(seb) <- rownames(z) <- rownames(p) <- rownames(n) <- GAlist$rsids
  for (study in studyID) {
   message(paste("Extracting data from study:",study))
   stat <- fread(GAlist$studyfiles[study], data.table=FALSE)
   if(is.null(stat[["n"]])) stat$n <- rep(GAlist$study$neff[study],length(stat$b))
   b[stat$rsids,study] <- stat$b
   seb[stat$rsids,study] <- stat$seb
   z[stat$rsids,study] <- stat$b/stat$seb
   p[stat$rsids,study] <- stat$p
   n[stat$rsids,study] <- stat$n
  }
  if(!is.null(rsids)) rsids <- rsids[rsids%in%GAlist$rsids]
  if(is.null(rsids)) rsids <- stat$rsids[stat$rsids%in%GAlist$rsids]
  rsids <- match(rsids,GAlist$rsids)
  if(rm.na) return(list(b=na.omit(b[rsids,]),seb=na.omit(seb[rsids,]),z=na.omit(z[rsids,]),
                        p=na.omit(p[rsids,]), n=na.omit(n[rsids,]) ))
  if(!rm.na) return(list(b=b[rsids,],seb=seb[rsids,],z=z[rsids,],p=p[rsids,],n=n[rsids,] ))
 }

 if(format=="data.frame" && what=="all") {
  if(length(studyID)>1) stop("Only one study allowed")
  study <- studyID
  message(paste("Extracting data from study:",study))
  stat <- fread(GAlist$studyfiles[study], data.table=FALSE)
  if(is.null(stat[["n"]])) stat$n <- rep(GAlist$study$neff[study],nrow(stat))
  #if(is.null(stat[["z"]])) stat$z <- stat$b/stat$seb

  if(!is.null(rsids)) rsids <- rsids[rsids%in%stat$rsids]
  if(is.null(rsids)) rsids <- stat$rsids
  rsids <- match(rsids,stat$rsids)
  if(rm.na) return(na.omit(stat[rsids,]))
  if(!rm.na) return(na.omit(stat[rsids,]))
 }

 if(what%in%c("rsids","b","seb","eaf","ea","nea","z","p")) {
  if(length(what)>1) stop("Only one feature allowed")
  res <- matrix(NA,ncol=length(studyID),nrow=length(GAlist$rsids))
  colnames(res) <- studyID
  rownames(res) <- GAlist$rsids
  for (study in studyID) {
   message(paste("Extracting data from study:",study))
   stat <- fread(GAlist$studyfiles[study], data.table=FALSE)
   if(what=="z") res[stat$rsids,study] <- stat$b/stat$seb
   if(!what=="z") res[stat$rsids,study] <- stat[,what]
  }
  if(!is.null(rsids)) rsids <- rsids[rsids%in%stat$rsids]
  if(is.null(rsids)) rsids <- stat$rsids
  rsids <- match(rsids,stat$rsids)
  res <- res[rsids,]

  if(rm.na) res <- na.omit(res)
  if(what=="rsids") res <- res[,1]
  return(res)
 }

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



#' Generate gene sets for pathway analysis
#'
#' Generates gene sets for pathway analysis based on a given set of SNPs, using a
#' pathway database contained in a GAlist object. The function
#' samples causal SNPs and constructs gene sets of a given size and number of
#' causal genes.
#'
#' @param GAlist a list object containing SNP and pathway
#'   data
#' @param mcausal an integer indicating the number of causal SNPs to sample.
#'   Default is 500.
#' @param causal a vector of SNP IDs representing the causal SNPs to be used in
#'   the gene set construction. Default is \code{NULL} and SNPs will be randomly
#'   sampled from the available SNPs in \code{GAlist}.
#' @param ngenes an integer or vector indicating the number of genes to include
#'   in each gene set. If a vector is supplied, a gene set will be generated for
#'   each value in the vector. Default is 10.
#' @param ncgenes an integer or vector indicating the number of causal genes to
#'   include in each gene set. If a vector is supplied, a gene set will be
#'   generated for each combination of values in \code{ngenes} and \code{ncgenes}.
#'   Default is 2.
#' @param ngsets an integer indicating the number of gene sets to generate.
#'   Default is 10.
#'
#' @return a list of gene sets, where each gene set is represented as a vector of
#'   SNP IDs.
#'
#'
#' @examples
#' \dontrun{
#' # Generate gene sets
#' sets <- gpath(GAlist = GAlist, ngenes = 10, ncgenes = 2, ngsets = 5)
#' }
#'
#' @export
#'
gpath <- function(GAlist = NULL,
                  mcausal=500, causal=NULL, overlap=FALSE,
                  ngenes=10, ncgenes=2, ngsets=10, format="markers"){
 # mcausal is the number of causal markers
 # causal is a vector of causal markers
 # ngenes is the number (or vector) of genes in pathway
 # ncgenes is the number (or vector) of causal genes in pathway
 # ngsets is the number of pathways

 rsids <- unlist(GAlist$rsids)

 if(min(ngenes)-max(ncgenes)<0) stop("Number of causal genes (ncgenes) larger than number of genes in pathway (ngenes)")

 # Sample causal markers
 if(is.null(causal)) causal <- sample(rsids, size = mcausal)
 mcausal <- length(causal)

 # Marker sets for genes in database
 sets <- getMarkerSetsDB(GAlist = GAlist, feature="Genes")

 # Identify genes with causal markers
 causal_sets <- sapply(sets, function(x){any(x %in% causal)})

 # These are causal genes
 causal_genes <- names(sets)[causal_sets]

 # These are non causal markers
 non_causal_genes <- names(sets)[!causal_sets]

 gsets <- NULL
 setnames <- NULL
 ijk <- 0
 for (i in 1:ngsets) {
  for (j in 1:length(ngenes)) {
   for (k in 1:length(ncgenes)) {
    ijk <- ijk + 1
    csets <- sample(causal_genes, ncgenes[k])
    ncsets <- sample(non_causal_genes, ngenes[j]-ncgenes[k])
    gsets[[ijk]] <- c(csets, ncsets)
    setnames <- c(setnames, paste(c(i,ngenes[j],ncgenes[k]),collapse="_"))
   }
  }
 }

 if(format=="markers") gsets <- lapply(gsets, function(x){unlist(sets[x])})
 names(gsets) <- setnames
 return(gsets)
}

#' @export
designSets <- function(sets=NULL, rsids=NULL) {
 sets <- qgg:::mapSets(sets=sets,rsids=rsids, index=TRUE)
 W <- matrix(0,nrow=length(rsids), ncol=length(sets))
 colnames(W) <- names(sets)
 rownames(W) <- rsids
 for(i in 1:length(sets)) {
  W[sets[[i]],i] <- 1
 }
 return(W)
}


