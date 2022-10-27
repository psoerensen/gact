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
#' @param GACTdb list object providing information about genomic associations for different traits
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
gact <- function(version="gact-0.01", task="download", wkdir=NULL, what="lite") {
 if(is.null(wkdir)) wkdir <- getwd()
 #if(dir.exists(wkdir))
 #dir.create(wkdir)


 library(org.Hs.eg.db)
 library(reactome.db)

 # load Glist
 Glist_full_name <- paste0(wkdir,"/glist/",list.files(path="./glist",pattern="Glist"))
 Glist <- readRDS(Glist_full_name)

 # update Glist$ldfiles
 ldfiles <- list.files(path=paste0(wkdir,"/glist"),pattern=".ld")
 rws <- sapply(ldfiles,function(x){grep(x,Glist$ldfiles)})
 rws <- order(rws)
 Glist$ldfiles <- paste0(wkdir,"/glist/",ldfiles[rws])

 gactdb <- NULL
 gactdb$version <- version

 gactdb$traits <- c("t2d","cad")
 gactdb$dirs <- list.dirs()


 glistfile <- paste0(wkdir,"/glist/",list.files(path="./glist",pattern="Glist"))
 gactdb$glistfile <- glistfile

 features <- c("Marker","Genes","Proteins","GO","Pathways",
               "ProteinComplexes","ChemicalComplexes")
 gactdb$features <- features

 featurefiles <- paste0(wkdir,"/statistics/gsea",features,".rds")
 names(featurefiles) <- features
 gactdb$featurefiles <- featurefiles

 ensg2eg <- as.list(org.Hs.egENSEMBL2EG)
 eg2ensg <- org.Hs.egENSEMBL
 mapped_genes <- mappedkeys(eg2ensg)
 eg2ensg <- as.list(eg2ensg[mapped_genes])
 eg2sym <- org.Hs.egSYMBOL
 mapped_genes <- mappedkeys(eg2sym)
 eg2sym <- as.list(eg2sym[mapped_genes])
 ensg2sym <- sapply(ensg2eg, function(x){paste(unlist(eg2sym[x], use.names=FALSE),collapse=" ")})

 gactdb$ensg2eg <- ensg2eg
 gactdb$eg2ensg <- eg2ensg
 gactdb$eg2sym <- eg2sym
 gactdb$ensg2sym <- ensg2sym

 gactdb$pathways <- as.list(reactomePATHID2EXTID)

 return(gactdb)
}

#' @export
#'
getStat <- function(GACTdb=NULL, feature=NULL, featureID=NULL,
                    studyID=NULL, trait="T2D", threshold=1,
                    format="data.frame", hyperlink=FALSE, cls=NULL) {

 features <- c("Marker","Genes","Proteins","GO","Pathways",
               "ProteinComplexes","ChemicalComplexes")
 header <- c("Marker ID","Gene ID","Protein ID","GO ID","Pathway ID",
             "Protein ID","Chemical ID")
 names(header) <- features

 if(!feature%in%GACTdb$features) stop(paste("feature:",feature,"not in GACT database"))
 res <- readRDS(GACTdb$featurefiles[feature])
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
 res <- as.data.frame(res)

 if(feature=="Genes") {
  if(!hyperlink) {
   res <- cbind(rownames(res),GACTdb$ensg2sym[rownames(res)], res)
   colnames(res)[1:2] <- c("Ensembl Gene ID","Symbol")
  }
  if(hyperlink) {
   res <- cbind(rownames(res),rownames(res),GACTdb$ensg2sym[rownames(res)], res)
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
 return(res)
}

#' @export
#'
writeStat <- function(GACTdb=NULL, feature=NULL, featureID=NULL,
                      studyID=NULL, trait="T2D", threshold=1,
                      format="data.frame", file.csv=NULL, hyperlink=TRUE) {
 stat <- getStat(GACTdb=GACTdb, feature=feature, featureID=featureID,
                 studyID=studyID, trait=trait, threshold=threshold,
                 format=format, hyperlink=hyperlink)
 write.csv2(stat,file=file.csv,row.names=FALSE)
}
