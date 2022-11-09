#######################################################################################
# Create GACT database functions
#######################################################################################
#'
#' Create or download GACT data
#'
#' @description
#' The gact function is used to create or download a GACT database with information
#' on genomic assocations for complex traits. The genomic associations
#' include single marker associations mapped to genetic markers in a Glist.
#' This includes stringent quality control. In addition genomic associations
#' linked to different genomic features (e.g. genes, proteins, chemical-complexes,
#' protein-complexes, biological pathways) are provided.
#'
#' @param GAlist list object providing information about genomic associations for different traits
#' @param feature name of feature including "Marker","Genes","Proteins","GO","Pathways","ProteinComplexes","ChemicalComplexes"
#' @param featureID is the feature specific IDs such as "GO:0000002"
#' @param studyID is the study specific ID (not used currently)
#' @param trait is the trait or disease name such as "t2d"
#' @param threshold is the p-value (or posterior inclusion probability) for which information is extracted
#' @param format output format which currently is a data.frame
#' @param hyperlink logical if TRUE then a feature specific hyperlink is provided in the data frame
#' @param file.csv is the name of the csvs file used for writing the data
'
#'
#' @return Returns a data frame with genomic associations for a specific feature

#' @author Peter Soerensen

#' @export
#'
createDB <- function(Glist=NULL, version=NULL, dbdir=NULL, what="lite") {

 if(is.null(version)) stop(paste("Please include a database name using the version argument"))
 #if(is.null(Glist)) stop(paste("Please include a Glist using the Glist argument"))
 dbdir <- paste0(dbdir,"/",version)
 gstatdir <- paste0(dbdir,"/gstat/")
 gsetsdir <- paste0(dbdir,"/gsets/")
 gseadir <- paste0(dbdir,"/gsea/")
 glistdir <- paste0(dbdir,"/glist/")
 ldscdir <- paste0(dbdir,"/ldsc/")
 gbayesdir <- paste0(dbdir,"/gbayes/")
 markerdir <- paste0(dbdir,"/marker/")
 rawdir <- paste0(dbdir,"/raw/")
 if(dir.exists(dbdir)) stop(paste("Directory:",dbdir,"allready exists"))
 if(!dir.exists(dbdir)) {
  dir.create(dbdir)
  dir.create(glistdir)
  dir.create(gstatdir)
  dir.create(gsetsdir)
  dir.create(gseadir)
  dir.create(ldscdir)
  dir.create(gbayesdir)
  dir.create(markerdir)
  dir.create(rawdir)
 }

 GAlist <- NULL
 GAlist$version <- version

 GAlist$traits <- NULL

 GAlist$dirs <- c(glistdir,gstatdir,gsetsdir,gseadir, ldscdir, gbayesdir, markerdir, rawdir)
 names(GAlist$dirs) <- c("glist","gstat","gsets","gsea", "ldsc", "gbayes", "marker", "raw")

 # features in the database
 GAlist$features <- c("Markers","Genes","Proteins","GO","Pathways",
                      "ProteinComplexes","ChemicalComplexes")

 if(!is.null(Glist)) {
  keep <- unlist(Glist$rsids)%in%unlist(Glist$rsidsLD)
  GAlist$markers <- data.frame(rsids=unlist(Glist$rsids),
                               chr=unlist(Glist$chr),
                               pos=unlist(Glist$pos),
                               ea=unlist(Glist$a1),
                               nea=unlist(Glist$a2),
                               eaf=unlist(Glist$af),
                               stringsAsFactors=FALSE)[keep,]
  GAlist$rsids <- unlist(Glist$rsids)[keep]
  GAlist$cpra <- unlist(Glist$cpra)[keep]
  file_markers <- paste0(GAlist$dirs["marker"],"markers.txt.gz")
  fwrite(GAlist$markers, file=file_markers)
 }

 # update Glist$ldfiles
 #ldfiles <- list.files(path=paste0(dbdir,"/glist/ldfiles"),pattern=".ld")
 #rws <- sapply(ldfiles,function(x){grep(x,Glist$ldfiles)})
 #rws <- order(rws)
 #Glist$ldfiles <- paste0(dbdir,"/glist/",ldfiles[rws])

 return(GAlist)
}

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
            "https://www.dropbox.com/s/9ah6aw0fborrp0z/string2ensg.rds?dl=1",
            "https://www.dropbox.com/s/ny94ibdbqhtg62h/stitch.rds?dl=1",
            "https://www.dropbox.com/s/q83q3mnvos8wdxk/reactome2ensg.rds?dl=1",
            "https://www.dropbox.com/s/9ah6aw0fborrp0z/string2ensg.rds?dl=1",
            "https://www.dropbox.com/s/7gj36rdec6spk9u/stitch2ensg.rds?dl=1")

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
                   "stitch2ensg.rds")
  for (feature in names(urls)) {
   message(paste("Downloading file:",feature))
   destfile <- paste0(GAlist$dirs["gsets"],feature)
   download.file(url=urls[feature], mode = "wb", dest=destfile)
  }
  GAlist$gsetsfiles <- paste0(GAlist$dirs["gsets"],names(urls))
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

  names(urls) <- c("gseaGenes.rds",
                   "gseaGO.rds",
                   "gseaPathways.rds",
                   "gseaProteinComplexes.rds",
                   "gseaChemicalComplexes.rds")

  for (feature in names(urls)) {
   message(paste("Downloading file:",feature))
   destfile <- paste0(GAlist$dirs["gsea"],feature)
   download.file(url=urls[feature], mode = "wb", dest=destfile)
  }
  GAlist$gseafiles <- paste0(GAlist$dirs["gsea"],names(urls))
  names(GAlist$gseafiles) <- gsub(".rds","",names(urls))
  names(GAlist$gseafiles) <- gsub("gsea","",names(GAlist$gseafiles))

 }
 if(what=="gstat") {
  # download gstat files in the database
  message("Downloading GWAS summary statistics")
  url_stat <- "https://www.dropbox.com/s/qrcivih31iuuril/stat.rds?dl=1"
  destfile <- paste0(GAlist$dirs["gstat"],"gstat.rds")
  download.file(url=url_stat, mode = "wb", dest=destfile)
  GAlist$gstatfiles <- paste0(GAlist$dirs["gstat"],"gstat.rds")

  url_stat <- "https://www.dropbox.com/s/jgdragns3a77wr4/GWAS_information.csv?dl=1"
  destfile <- paste0(GAlist$dirs["gstat"],"GWAS_information.csv")
  download.file(url=url_stat, mode = "wb", dest=destfile)
  GAlist$study <- as.list(fread(destfile,data.table=FALSE))

  url_stat <- c("https://www.dropbox.com/s/iqd7c4bbds03nrg/GWAS1.txt.gz?dl=1",
                "https://www.dropbox.com/s/pdv1fg280n86dwg/GWAS2.txt.gz?dl=1",
                "https://www.dropbox.com/s/t3y05ex1uouo4qg/GWAS3.txt.gz?dl=1",
                "https://www.dropbox.com/s/34yd6sxltqiv6e8/GWAS4.txt.gz?dl=1",
                "https://www.dropbox.com/s/xyfadehraaajkol/GWAS5.txt.gz?dl=1",
                "https://www.dropbox.com/s/2b5u6m60l6a5e82/GWAS6.txt.gz?dl=1")

  for(study in 1:length(url_stat)) {
   destfile <- paste0(GAlist$dirs["gstat"],substring(url_stat[study],43,nchar(url_stat[study])-5))
   download.file(url=url_stat[study], mode = "wb", dest=destfile)
  }


 }
 if(what=="marker") {
  # download marker files in the database
  message("Downloading marker information")
  url_stat <- "https://www.dropbox.com/s/4k54owkby3uf2hf/markers.txt.gz?dl=1"
  destfile <- paste0(GAlist$dirs["marker"],"markers.txt.gz")
  download.file(url=url_stat, mode = "wb", dest=destfile)
  GAlist$markerfiles <-paste0(GAlist$dirs["marker"],"markers.txt.gz")
  GAlist$markers <- fread(GAlist$markerfiles, data.table=FALSE)
  GAlist$rsids <- GAlist$markers$rsids
  GAlist$cpra <- paste(GAlist$markers$chr,
                       GAlist$markers$pos,
                       GAlist$markers$ea,
                       GAlist$markers$nea,sep="_")
 }



 return(GAlist)
}



#' @export
#'
gact <- function(GAlist=NULL, version="t2dm-gact-0.0.1", task="download",
                 dbdir=NULL, what="lite") {

 if(is.null(dbdir)) dbdir <- getwd()

 if(task=="download") {
  GAlist <- createDB(Glist=NULL, version=version, dbdir=dbdir)

  # Step 2: Download data from database:
  GAlist <- downloadDB(GAlist=GAlist, what="marker")
  GAlist <- downloadDB(GAlist=GAlist, what="gsets")
  GAlist <- downloadDB(GAlist=GAlist, what="gsea")
  GAlist <- downloadDB(GAlist=GAlist, what="gstat")

 }
 # if(task=="download") {
 #
 #  dbdir <- paste0(dbdir,"/",version)
 #  gseadir <- paste0(dbdir,"/gsea/")
 #  gstatdir <- paste0(dbdir,"/gstat/")
 #  gsetsdir <- paste0(dbdir,"/gsets/")
 #  if(dir.exists(dbdir)) stop(paste("Directory:",dbdir,"allready exists"))
 #  if(!dir.exists(dbdir)) {
 #   dir.create(dbdir)
 #   dir.create(gseadir)
 #   dir.create(gstatdir)
 #   dir.create(gsetsdir)
 #  }
 #
 #  # download feature files in the database
 #  features <- c("Genes","GO","Pathways",
 #                "ProteinComplexes","ChemicalComplexes")
 #
 #  for( feature in features) {
 #   url <- paste0("https://github.com/psoerensen/gdtdb/raw/main/",version,"/gsea/gsea",feature,".rds")
 #   destfile <- paste0(gseadir,"gsea",feature,".rds")
 #   download.file( url=url, mode = "wb",  destfile=destfile)
 #  }
 #
 #  message("Downloading summary statistics")
 #
 #  url_stat <- "https://www.dropbox.com/s/qrcivih31iuuril/stat.rds?dl=1"
 #  destfile <- paste0(gstatdir,"gstat.rds")
 #  download.file(url=url_stat, mode = "wb", dest=destfile)
 #
 #
 #  urls <- c("https://www.dropbox.com/s/ijtc7l6hgpaieo1/eg2rsids_10kb.rds?dl=1",
 #            "https://www.dropbox.com/s/0aqbqa7ihrg6i2e/ensg2rsids_10kb.rds?dl=1",
 #            "https://www.dropbox.com/s/p3ut5dwfx0zw4v1/ensp2rsids_10kb.rds?dl=1",
 #            "https://www.dropbox.com/s/1py37zd92ttsvnp/ensg2sym.rds?dl=1",
 #            "https://www.dropbox.com/s/2ggu4u5hp406cif/go.rds?dl=1",
 #            "https://www.dropbox.com/s/uryyxnjyhxa9azf/reactome.rds?dl=1",
 #            "https://www.dropbox.com/s/9ah6aw0fborrp0z/string2ensg.rds?dl=1",
 #            "https://www.dropbox.com/s/ny94ibdbqhtg62h/stitch.rds?dl=1",
 #            "https://www.dropbox.com/s/q83q3mnvos8wdxk/reactome2ensg.rds?dl=1",
 #            "https://www.dropbox.com/s/9ah6aw0fborrp0z/string2ensg.rds?dl=1",
 #            "https://www.dropbox.com/s/7gj36rdec6spk9u/stitch2ensg.rds?dl=1")
 #
 #  names(urls) <- c("eg2rsids_10kb.rds",
 #                   "ensg2rsids_10kb.rds",
 #                   "ensp2rsids_10kb.rds",
 #                   "ensg2sym.rds",
 #                   "go.rds",
 #                   "reactome.rds",
 #                   "string.rds",
 #                   "stitch.rds",
 #                   "reactome2ensg.rds",
 #                   "string2ensg.rds",
 #                   "stitch2ensg.rds")
 #
 #  for (feature in names(urls)) {
 #   message(paste("Downloading file:",feature))
 #   destfile <- paste0(gsetsdir,feature)
 #   download.file(url=urls[feature], mode = "wb", dest=destfile)
 #  }
 #
 #
 #  GAlist <- NULL
 #  GAlist$version <- version
 #
 #  GAlist$traits <- c("t2d")
 #  GAlist$dirs <- list.dirs(dbdir)
 #
 #  GAlist$features <- features
 #
 #  featurefiles <- paste0(gseadir,"gsea",features,".rds")
 #  names(featurefiles) <- features
 #  GAlist$featurefiles <- featurefiles
 #  GAlist$gsetsfiles <- paste0(gsetsdir,names(urls))
 #  GAlist$gstatfiles <- paste0(gstatdir,"gstat.rds")
 #
 #
 #
 #  # load Glist
 #  glistfile <- paste0(dbdir,"/glist/",list.files(path="./glist",pattern="Glist"))
 #  if(file.exists(glistfile))  {
 #   Glist <- readRDS(glistfile)
 #   GAlist$glistfile <- glistfile
 #
 #   # update Glist$ldfiles
 #   ldfiles <- list.files(path=paste0(dbdir,"/glist/ldfiles"),pattern=".ld")
 #   rws <- sapply(ldfiles,function(x){grep(x,Glist$ldfiles)})
 #   rws <- order(rws)
 #   Glist$ldfiles <- paste0(dbdir,"/glist/",ldfiles[rws])
 #  }
 #
 #  GAlist$ensg2sym <- readRDS(paste0(gsetsdir,"ensg2sym.rds"))
 #
 #  GAlist$gsets <- vector(mode = "list", length = length(GAlist$gsetsfiles))
 #  for(i in 1:length(GAlist$gsetsfiles)) {
 #   GAlist$gsets[[i]] <- readRDS(GAlist$gsetsfile[i])
 #  }
 #
 #  stat <- readRDS(GAlist$gstatfiles)
 #
 #  GAlist$gsets[[1]] <- qgg:::mapSets(sets=GAlist$gsets[[1]], rsids=stat$rsids, index=FALSE)
 #  GAlist$gsets[[2]] <- qgg:::mapSets(sets=GAlist$gsets[[2]], rsids=stat$rsids, index=FALSE)
 #  GAlist$gsets[[3]] <- qgg:::mapSets(sets=GAlist$gsets[[3]], rsids=stat$rsids, index=FALSE)
 #
 #  names(GAlist$gsets) <- c("eg2rsids",
 #                   "ensg2rsids",
 #                   "ensp2rsids",
 #                   "ensg2sym",
 #                   "go",
 #                   "reactome",
 #                   "string",
 #                   "stitch",
 #                   "reactome2ensg",
 #                   "string2ensg",
 #                   "stitch2ensg")
 #
 # }

 return(GAlist)
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

#' @export
#'
writeStat <- function(GAlist=NULL, feature=NULL, featureID=NULL,
                      studyID=NULL, trait="T2D", threshold=1,
                      format="data.frame", file.csv=NULL, hyperlink=TRUE) {
 stat <- getStat(GAlist=GAlist, feature=feature, featureID=featureID,
                 studyID=studyID, trait=trait, threshold=threshold,
                 format=format, hyperlink=hyperlink)
 write.csv2(stat,file=file.csv,row.names=FALSE)
}

#' @export
#'
getStudies <- function(GAlist=NULL) {
 return(as.data.frame(GAlist$study))
}


#' @export
#'
getSets <- function(GAlist=NULL, feature=NULL, featureID=NULL) {
 sets <- NULL
 if(feature=="Entres Genes") sets <- GAlist$gsets[[1]]
 if(feature=="Genes") sets <- GAlist$gsets[[2]]
 if(feature=="Proteins") sets <- GAlist$gsets[[3]]
 if(feature=="Gene Symbol") sets <- GAlist$gsets[[4]]
 if(feature=="GO") sets <- GAlist$gsets[[5]]
 if(feature=="Pathways") sets <- GAlist$gsets[[6]]
 if(feature=="ProteinComplexes") sets <- GAlist$gsets[[7]]
 if(feature=="ChemicalComplexes") sets <- GAlist$gsets[[8]]
 if(feature=="Pathways2Genes") sets <- GAlist$gsets[[9]]
 if(feature=="ProteinComplexes2Genes") sets <- GAlist$gsets[[10]]
 if(feature=="ChemicalComplexes2Genes") sets <- GAlist$gsets[[11]]
 if(!is.null(featureID)) {
  inSet <- featureID%in%names(sets)
  if(any(!inSet)) warning(paste("Some IDs not in data base:",featureID[!inSet]))
  featureID <- featureID[featureID%in%names(sets)]
  sets <- unlist(sets[featureID])
 }
 return(sets)
 }


#' @export
#'
updateStatDB <- function(GAlist=NULL,
                         stat=NULL,
                         source=NULL,
                         trait="unknown",
                         type="unknown",
                         gender="unknown",
                         n=NULL,
                         ncase=NULL,
                         ncontrol=NULL,
                         reference="unknown") {

 message("Perform quality control of external summary statistics")
 stat <- qcStatDB(GAlist=GAlist,stat=stat)

 message("Collecting information on external summary statistics")
 study_number <- length(GAlist$study$id)+1
 studyID <- paste0("GWAS",study_number)
 GAlist$study$id[study_number] <- studyID
 GAlist$study$file[study_number] <- paste0(studyID,".txt")
 GAlist$study$trait[study_number] <- trait
 GAlist$study$type[study_number] <- type
 GAlist$study$gender[study_number] <- gender
 GAlist$study$n[study_number] <- n
 GAlist$study$ncase[study_number] <- ncase
 GAlist$study$ncontrol[study_number] <- ncontrol
 GAlist$study$reference[study_number] <- reference
 GAlist$study$source[study_number] <- source
 message(paste("Writing processed summary statistics til internal file:",
               GAlist$study$file[study_number]))
 file_stat <- paste0(GAlist$dirs["gstat"],GAlist$study$file[study_number],".gz")
 if(file.exists(file_stat)) stop(paste("GWAS summary statistics file allready exists:",
                                       file_stat))
 fwrite(stat, file_stat)
 file_stat_information <- paste0(GAlist$dirs["gstat"],"GWAS_information.csv")
 #fwrite(as.data.frame(GAlist$study), file_stat_information)
 write.csv2(as.data.frame(GAlist$study),file=file_stat_information,row.names=FALSE)
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
 if(all(fm_internal[1:9]%in%colnames(stat))) format <- "internal"

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
 if(!is.null(stat$eaf))effect_allele_freq <- stat[,"eaf"]

 # flip is not aligned
 stat[!aligned,"b"] <- -effect[!aligned]
 stat[!aligned,"ea"] <- non_effect_allele[!aligned]
 stat[!aligned,"nea"] <- effect_allele[!aligned]
 if(!is.null(stat$eaf)) {
  stat[!aligned,"eaf"] <- 1-effect_allele_freq[!aligned]
  excludeMAFDIFF <- abs(marker$eaf-stat$eaf) > excludeMAFDIFF
 }

 message(paste("Number of markers excluded by large difference between MAF difference:", sum(excludeMAFDIFF)))
 message("")

 if(!is.null(stat$eaf)) stat <- stat[!excludeMAFDIFF,]
 if(!is.null(stat$eaf)) marker <- marker[!excludeMAFDIFF,]
 if(is.null(stat$n)) stat$n <- neff(seb=stat$seb,af=stat$eaf)
 colnames(stat)[1] <- "rsids"
 return(stat)
}
