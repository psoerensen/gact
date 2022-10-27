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
   res <- cbind(rownames(res),gactdb$ensg2sym[rownames(res)], res)
   colnames(res)[1:2] <- c("Ensembl Gene ID","Symbol")
  }
  if(hyperlink) {
   res <- cbind(rownames(res),rownames(res),gactdb$ensg2sym[rownames(res)], res)
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
