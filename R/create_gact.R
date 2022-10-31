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
gact <- function(GAlist=NULL, version="t2dm-gact-0.0.1", task="download", wkdir=NULL, what="lite") {

 if(is.null(wkdir)) wkdir <- getwd()

 if(task=="download") {

  #version="t2dm-gact-0.0.1"
  #wkdir <- "C:/Users/au223366/Dropbox/Projects/balder"
  if(is.null(wkdir)) wkdir <- getwd()
  dbdir <- paste0(wkdir,"/",version)
  gtestdir <- paste0(wkdir,"/",version,"/gtest/")
  gstatdir <- paste0(wkdir,"/",version,"/gstat/")
  gsetsdir <- paste0(wkdir,"/",version,"/gsets/")
  if(dir.exists(dbdir)) stop(paste("Directory:",dbdir,"allready exists"))
  if(!dir.exists(dbdir)) {
   dir.create(dbdir)
   dir.create(gtestdir)
   dir.create(gstatdir)
   dir.create(gsetsdir)
  }

  # download feature files in the database
  features <- c("Genes","GO","Pathways",
                "ProteinComplexes","ChemicalComplexes")

  for( feature in features) {
   url <- paste0("https://github.com/psoerensen/gdtdb/raw/main/",version,"/gtest/gsea",feature,".rds")
   destfile <- paste0(gtestdir,"gsea",feature,".rds")
   download.file( url=url, mode = "wb",  destfile=destfile)
  }

  message("Downloading summary statistics")

  url_stat <- "https://www.dropbox.com/s/qrcivih31iuuril/stat.rds?dl=1"
  destfile <- paste0(gstatdir,"gstat.rds")
  download.file(url=url_stat, mode = "wb", dest=destfile)


  urls <- c("https://www.dropbox.com/s/ijtc7l6hgpaieo1/eg2rsids_10kb.rds?dl=1",
            "https://www.dropbox.com/s/0aqbqa7ihrg6i2e/ensg2rsids_10kb.rds?dl=1",
            "https://www.dropbox.com/s/p3ut5dwfx0zw4v1/ensp2rsids_10kb.rds?dl=1",
            "https://www.dropbox.com/s/1py37zd92ttsvnp/ensg2sym.rds?dl=1",
            "https://www.dropbox.com/s/2ggu4u5hp406cif/go.rds?dl=1",
            "https://www.dropbox.com/s/uryyxnjyhxa9azf/reactome.rds?dl=1",
            "https://www.dropbox.com/s/v5d4pvdwnerfq7o/string.rds?dl=1",
            "https://www.dropbox.com/s/ny94ibdbqhtg62h/stitch.rds?dl=1")

  names(urls) <- c("eg2rsids_10kb.rds",
                   "ensg2rsids_10kb.rds",
                   "ensp2rsids_10kb.rds",
                   "ensg2sym.rds",
                   "go.rds",
                   "reactome.rds",
                   "string.rds",
                   "stitch.rds")

  for (feature in names(urls)) {
   message(paste("Downloading file:",feature))
   destfile <- paste0(gsetsdir,feature)
   download.file(url=urls[feature], mode = "wb", dest=destfile)
  }


  GAlist <- NULL
  GAlist$version <- version

  GAlist$traits <- c("t2d")
  GAlist$dirs <- list.dirs(dbdir)

  GAlist$features <- features

  featurefiles <- paste0(gtestdir,"gsea",features,".rds")
  names(featurefiles) <- features
  GAlist$featurefiles <- featurefiles
  GAlist$gsetsfiles <- paste0(gsetsdir,names(urls))
  GAlist$gstatfiles <- paste0(gstatdir,"gstat.rds")



  # load Glist
  glistfile <- paste0(dbdir,"/glist/",list.files(path="./glist",pattern="Glist"))
  if(file.exists(glistfile))  {
   Glist <- readRDS(glistfile)
   GAlist$glistfile <- glistfile

   # update Glist$ldfiles
   ldfiles <- list.files(path=paste0(wkdir,"/glist/ldfiles"),pattern=".ld")
   rws <- sapply(ldfiles,function(x){grep(x,Glist$ldfiles)})
   rws <- order(rws)
   Glist$ldfiles <- paste0(wkdir,"/glist/",ldfiles[rws])
  }

  GAlist$ensg2sym <- readRDS(paste0(gsetsdir,"ensg2sym.rds"))

 }

 if(task=="prepare") {
  GAlist <- NULL
  GAlist$version <- version

  GAlist$traits <- c("t2d","cad")
  GAlist$dirs <- list.dirs()

  # load Glist
  glistfile <- paste0(wkdir,"/glist/",list.files(path="./glist",pattern="Glist"))
  if(file.exists(glistfile))  Glist <- readRDS(glistfile)
  GAlist$glistfile <- glistfile

  # update Glist$ldfiles
  ldfiles <- list.files(path=paste0(wkdir,"/glist/ldfiles"),pattern=".ld")
  rws <- sapply(ldfiles,function(x){grep(x,Glist$ldfiles)})
  rws <- order(rws)
  Glist$ldfiles <- paste0(wkdir,"/glist/",ldfiles[rws])

  # features in the database
  features <- c("Marker","Genes","Proteins","GO","Pathways",
                "ProteinComplexes","ChemicalComplexes")

  GAlist$features <- features

  featurefiles <- paste0(wkdir,"/statistics/gsea",features,".rds")
  names(featurefiles) <- features
  GAlist$featurefiles <- featurefiles

 }

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
  rownames(res) <- res$rsids
  return(res)
 }

 if(!feature%in%GAlist$features) stop(paste("feature:",feature,"not in GACT database"))



 res <- readRDS(GAlist$featurefiles[feature])
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
    res <- cbind(rownames(res),GAlist$ensg2sym[rownames(res)], res)
    colnames(res)[1:2] <- c("Ensembl Gene ID","Symbol")
   }
   if(hyperlink) {
    res <- cbind(rownames(res),rownames(res),GAlist$ensg2sym[rownames(res)], res)
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
