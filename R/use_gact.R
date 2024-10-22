#######################################################################################
# GACT database functions
#######################################################################################
#'
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
#' stat <- getFeatureStat(GAlist = GAlist, feature = "Markers", studyID="GWAS1")
#'
#' # get gene statistics
#' stat <- getFeatureStat(GAlist = GAlist, feature = "Genes")
#'
#' # get gene statistics for a specific feature
#' stat <- getFeatureStat(GAlist = GAlist, feature = "Genes", featureID = "TP53")
#'
#' # get gene statistics filtered by a threshold
#' stat <- getFeatureStat(GAlist = GAlist, feature = "Genes", threshold = 0.05)
#' }


#' @export
#'
getFeatureStat <- function(GAlist=NULL, feature=NULL, featureID=NULL,file=NULL,
                    studyID=NULL, trait=NULL, threshold=0.95,
                    format="list") {
 features <- c("Markers","Genes","Proteins","GO","Pathways",
               "ProteinComplexes","ChemicalComplexes","Chromosomes", "VEGAS")
 header <- c("Marker ID","Gene ID","Protein ID","GO ID","Pathway ID",
             "Protein ID","Chemical ID", "Chromosome ID", "VEGAS")
 names(header) <- features
 if(!feature%in%features) stop(paste("feature:",feature,"not in database"))

 # #gseafile <- paste0(GAlist$dirs["gsea"],"ct_gsea",feature,"_gdtdb.rds")
 # gseafile <- GAlist$gseafiles[paste0("ct_gsea",feature,"_gdtdb")]
 # res <- readRDS(gseafile)
 # message(paste("Extract statistics based p-value threshold:",threshold))
 # cls5 <- grep("_0.05", colnames(res$stat))
 # cls95 <- grep("_0.95", colnames(res$stat))
 # cls <- 1:ncol(res$stat)
 # if(threshold==0.05) cls <- cls5
 # if(threshold==0.95) cls <- cls95
 #
 # colnames(res$stat) <- gsub("z_","",colnames(res$stat))
 # colnames(res$p) <- gsub("z_","",colnames(res$p))
 # colnames(res$stat) <- gsub("_0.05","",colnames(res$stat))
 # colnames(res$p) <- gsub("_0.05","",colnames(res$p))
 # colnames(res$stat) <- gsub("_0.95","",colnames(res$stat))
 # colnames(res$p) <- gsub("_0.95","",colnames(res$p))
 # res$p[res$stat==0] <- 1
 # rws <- rep(TRUE,lenth=nrow(res$stat))
 # if(!is.null(featureID)) rws <- rownames(res$stat)%in%featureID
 # if(sum(rws)==0) stop("None of featureIDs found in database")
 # res$stat <- res$stat[rws,cls]
 # res$p <- res$p[rws,cls]
 # if(!is.null(studyID)) {
 #  res$stat <- res$stat[,colnames(res$stat)%in%studyID]
 #  res$p <- res$p[,colnames(res$p)%in%studyID]
 # }

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
 res <- getFeatureStat(GAlist=GAlist, feature=feature, featureID=featureID,
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
getFeatureStat <- function(GAlist=NULL, feature=NULL, featureID=NULL,
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
designMatrixDB <- function(GAlist=NULL, feature=NULL, featureID=NULL, rowFeatureID=NULL, scale=FALSE) {
 if(is.null(GAlist)) stop ("Please provide GAlist")
 if(is.null(feature)) stop ("Please provide feature")
 sets <- getFeatureSets(GAlist=GAlist, feature=feature)
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
  if(scale) W[,i] <- scale(W[,i])
 }
 return(W)
}

#' Retrieve Feature Sets from Database
#'
#' The `getFeatureSets` function retrieves sets of features from a `GAlist` object based on the specified `feature` argument.
#' It allows filtering of these feature sets using the `featureID` argument. This function is designed to handle a variety of
#' feature types, including those related to genomic, proteomic, or pharmacological information. The successful operation
#' of this function depends on the structure and content of the `GAlist` parameter.
#'
#' @param GAlist A list object providing information and infrastructure of the gact database.
#' @param feature A character string specifying the type of data set to retrieve.
#'        Valid features include:
#'        \itemize{
#'          \item \code{"GO"}: Gene Ontology sets.
#'          \item \code{"Pathways"}: Pathway sets.
#'          \item \code{"ProteinComplexes"}: Protein complex sets.
#'          \item \code{"ChemicalComplexes"}: Chemical complex sets.
#'          \item \code{"DrugGenes"}: Drug-gene interaction sets.
#'          \item \code{"DrugATCGenes"}: Drug ATC gene sets.
#'          \item \code{"DrugComplexes"}: Drug complex sets.
#'          \item \code{"DiseaseGenes"}: Disease-gene interaction sets.
#'          \item \code{"DiseaseGenesEXP"}: Experimentally validated disease-gene sets.
#'          \item \code{"DiseaseGenesKB"}: Knowledge-based disease-gene sets.
#'          \item \code{"DiseaseGenesTM"}: Text-mined disease-gene sets.
#'          \item \code{"DiseaseGenesEXPplus"}: Comprehensive experimentally validated disease-gene sets.
#'          \item \code{"DiseaseGenesKBplus"}: Comprehensive knowledge-based disease-gene sets.
#'          \item \code{"DiseaseGenesTMplus"}: Comprehensive text-mined disease-gene sets.
#'          \item \code{"ATC1Genes"}, \code{"ATC2Genes"}, \code{"ATC3Genes"}, \code{"ATC4Genes"}: ATC classification gene sets.
#'          \item \code{"GTEx"}, \code{"GTExV7"}, \code{"GTExV8"}: GTEx project eQTL sets.
#'          \item \code{"GWAScatalog"}, \code{"GWAScatalogPlus"}: GWAS catalog sets.
#'          \item \code{"String"}: STRING database protein interaction sets.
#'          \item \code{"Stitch"}: STITCH database protein-chemical interaction sets.
#'        }
#' @param featureID Optional; a character vector of specific feature IDs to select.
#' @param minsets Optional; a numeric value specifying the minimum number of sets to include.
#' @param upstream Optional; a logical value indicating whether to include upstream data.
#' @param downstream Optional; a logical value indicating whether to include downstream data.
#' @param min_combined_score Optional; a numeric value specifying the minimum combined score for inclusion.
#' @param min_interactions Optional; a numeric value specifying the minimum number of interactions for inclusion.

#' @return A list of the specified feature sets. If `featureID` is not `NULL`,
#'   returns only the sets with IDs in `featureID`.
#'
#' @examples
#' \dontrun{
#' sets <- getFeatureSets(GAlist, feature="Genes")
#' getFeatureSets(GAlist, feature="Genes", featureID=c("gene1", "gene2"))
#' }
#'
#' @export
#'
getFeatureSets <- function(GAlist=NULL, feature=NULL, featureID=NULL, minsets=NULL,
                      upstream=FALSE,downstream=FALSE,
                      min_combined_score=700, min_interactions=5) {
 sets <- NULL
 if(feature=="GO") sets <- readRDS(file=file.path(GAlist$dirs["gsets"],"go.rds"))
 if(feature=="Pathways") sets <- readRDS(file=file.path(GAlist$dirs["gsets"],"reactome2ensg.rds"))
 if(feature=="ProteinComplexes") sets <- readRDS(file=file.path(GAlist$dirs["gsets"],"string2ensg.rds"))
 if(feature=="ChemicalComplexes") sets <- readRDS(file=file.path(GAlist$dirs["gsets"],"stitch2ensg.rds"))
 if(feature=="DrugGenes") sets <- readRDS(file.path(GAlist$dirs["gsets"],"drugGenes.rds"))
 if(feature=="DrugATCGenes") {
  hasATC <- !GAlist$target$ATC=="Unknown"
  sets <- split(GAlist$target$Target[hasATC],GAlist$target$Drug[hasATC])
  sets <- lapply(sets,unique)
 }
 if(feature=="DrugComplexes") sets <- readRDS(file.path(GAlist$dirs["gsets"],"drugComplex.rds"))
 if(feature=="DiseaseGenes") sets <- readRDS(file = file.path(GAlist$dirs["gsets"],"disease2ensg_human_disease_integrated_full.rds"))
 if(feature=="DiseaseGenesEXP") sets <- readRDS(file = file.path(GAlist$dirs["gsets"],"disease2ensg_human_disease_experiments_filtered.rds"))
 if(feature=="DiseaseGenesKB") sets <- readRDS(file = file.path(GAlist$dirs["gsets"],"disease2ensg_human_disease_knowledge_filtered.rds"))
 if(feature=="DiseaseGenesTM") sets <- readRDS(file = file.path(GAlist$dirs["gsets"],"disease2ensg_human_disease_textmining_filtered.rds"))
 if(feature=="DiseaseGenesEXPplus") sets <- readRDS(file = file.path(GAlist$dirs["gsets"],"disease2ensg_human_disease_experiments_full.rds"))
 if(feature=="DiseaseGenesKBplus") sets <- readRDS(file = file.path(GAlist$dirs["gsets"],"disease2ensg_human_disease_knowledge_full.rds"))
 if(feature=="DiseaseGenesTMplus") sets <- readRDS(file = file.path(GAlist$dirs["gsets"],"disease2ensg_human_disease_textmining_full.rds"))

 if(feature%in%c("ATC1Genes","ATC2Genes","ATC3Genes","ATC4Genes")) {
  if(feature=="ATC1Genes") atcSets <- readRDS(file = file.path(GAlist$dirs["gsets"], "atcSets1.rds"))
  if(feature=="ATC2Genes") atcSets <- readRDS(file = file.path(GAlist$dirs["gsets"], "atcSets2.rds"))
  if(feature=="ATC3Genes") atcSets <- readRDS(file = file.path(GAlist$dirs["gsets"], "atcSets3.rds"))
  if(feature=="ATC4Genes") atcSets <- readRDS(file = file.path(GAlist$dirs["gsets"], "atcSets4.rds"))
  drugSets <- readRDS(file=file.path(GAlist$dirs["gsets"],"drug2ensg.rds"))
  atcSets <- mapSetsDB(sets=atcSets, featureID=names(drugSets), index=FALSE)
  sets <- lapply(atcSets, function(x){
   unlist(drugSets[x])
  })
 }
 if (feature %in% c("GTEx", "GTExV7", "GTExV8")) {
  sets <- list()
  dbdir_suffix <- if (feature == "GTEx" || feature == "GTExV8") {
   "gtex/GTEx_Analysis_v8_eQTL"
  } else if (feature == "GTExV7") {
   "gtex/GTEx_Analysis_v7_eQTL"
  }
  dbdir <- file.path(GAlist$dbdir, dbdir_suffix)
  files <- list.files(dbdir, pattern = "egenes", full.names = TRUE)
  tissue <- gsub("\\.v[78]\\.egenes\\.txt\\.gz$", "", basename(files))
  for (i in seq_along(files)) {
   df <- fread(files[i], data.table = FALSE)
   df$gene_id <- substr(df$gene_id, 1, 15)
   df$variant_id <- gsub("chr|_b37|_b38", "", df$variant_id)
   rsid_column <- if (feature == "GTExV7") "rs_id_dbSNP147_GRCh37p13" else "rs_id_dbSNP151_GRCh38p7"
   sets[[i]] <- data.frame(ensg = df$gene_id,
                               rsids = df[,rsid_column],
                               cpra = df$variant_id,
                               p = df$qval)
   message(paste("Processing file:", basename(files[i])))
  }
  names(sets) <- tissue
 }
 if (feature %in% c("GWAScatalog", "GWAScatalogPlus")) {
  dbdir <- file.path(GAlist$dbdir, "gwas")
  gwasfile <- file.path(dbdir, "gwas-catalog-associations_ontology-annotated.tsv")
  gwas <- fread(gwasfile, data.table = FALSE, quote = "")
  processGenes <- function(genes) {
   genes <- lapply(genes, function(x) {
    set <- unlist(strsplit(x, split = ","))
    gsub(" ", "", set)
   })
   genes
  }
  gwasGenes <- split(gwas$SNP_GENE_IDS, f = gwas$MAPPED_TRAIT)
  gwasGenesUp <- split(gwas$UPSTREAM_GENE_ID, f = gwas$MAPPED_TRAIT)
  gwasGenesDown <- split(gwas$DOWNSTREAM_GENE_ID, f = gwas$MAPPED_TRAIT)
  gwasGenes <- processGenes(gwasGenes)
  gwasGenesUp <- processGenes(gwasGenesUp)
  gwasGenesDown <- processGenes(gwasGenesDown)
  # Check for mismatches in trait names
  if (any(!names(gwasGenes) == names(gwasGenesUp))) {
   stop("Mismatch detected in trait names (upstream)")
  }
  if (any(!names(gwasGenes) == names(gwasGenesDown))) {
   stop("Mismatch detected in trait names (downstream)")
  }
  # Combine genes if GWAScatalogPlus is selected
  if (feature == "GWAScatalogPlus") {
   gwasGenes <- mapply(c, gwasGenes, gwasGenesUp, gwasGenesDown, SIMPLIFY = FALSE)
  }
  sets <- lapply(gwasGenes, unique)
 }

 # if(feature%in%c("GWAScatalog","GWAScatalogPlus")) {
 #  dbdir <- file.path(GAlist$dbdir, "gwas")
 #  gwasfile <- file.path(dbdir, "gwas-catalog-associations_ontology-annotated.tsv")
 #  gwas <- fread(gwasfile, data.table=FALSE, quote="")
 #  gwasGenes <- split(gwas$SNP_GENE_IDS,f=gwas$MAPPED_TRAIT)
 #  gwasGenesUp <- split(gwas$UPSTREAM_GENE_ID,f=gwas$MAPPED_TRAIT)
 #  gwasGenesDown <- split(gwas$DOWNSTREAM_GENE_ID,f=gwas$MAPPED_TRAIT)
 #  gwasGenes <- lapply(gwasGenes, function(x) {
 #   set <- unlist(strsplit(x,split=","))
 #   set <- gsub(" ", "",set)
 #  })
 #  gwasGenesUp <- lapply(gwasGenesUp, function(x) {
 #   set <- unlist(strsplit(x,split=","))
 #   set <- gsub(" ", "",set)
 #  })
 #  gwasGenesDown <- lapply(gwasGenesDown, function(x) {
 #   set <- unlist(strsplit(x,split=","))
 #   set <- gsub(" ", "",set)
 #  })
 #  if(any(!names(gwasGenes)==names(gwasGenesUp))) stop("Mismatch detected in trait names")
 #  if(any(!names(gwasGenes)==names(gwasGenesDown))) stop("Mismatch detected in trait names")
 #  if(feature=="GWAScatalogPlus") {
 #   for(i in 1:length(gwasGenes)) {
 #    gwasGenes[[i]] <- c(gwasGenes[[i]],gwasGenesUp[[i]],gwasGenesDown[[i]])
 #   }
 #  }
 #  sets <- lapply(gwasGenes,unique)
 # }

 if(feature=="String") {
  file_string <- file.path(GAlist$dirs["gsets"],"9606.protein.links.v11.5.txt.gz")
  string <- fread(file_string, data.table=FALSE)
  string  <- string[string$combined_score>=min_combined_score,]
  string <- split( string$protein2,f=as.factor(string$protein1))
  sets <- string[sapply(string ,length)>=min_interactions]
 }
 if(feature=="Stitch") {
  file_stitch <- "http://stitch.embl.de/download/protein_chemical.links.v5.0/9606.protein_chemical.links.v5.0.tsv.gz"
  stitch <- fread(file_stitch, data.table=FALSE)
  stitch$protein <- gsub("9606.","",stitch$protein)
  stitch <- stitch[stitch$protein%in%names(GAlist$gsets$ensp2ensg),]
  stitch  <- stitch[stitch$combined_score>=min_combined_score,]
  stitch <- split( stitch$protein,f=as.factor(stitch$chemical))
  sets  <- stitch[sapply(stitch ,length)>=min_interactions]
 }

 if(!is.null(featureID)) {
  select <- names(sets)%in%featureID
  if(sum(select)==0) stop("None of the featureIDs found in sets")
  sets <- sets[select]
 }
 mset <- sapply(sets,length)
 if(!is.null(minsets)) sets <- sets[mset>minsets]
 return(sets)
}



#' @export
#'
getDrugComplexesDB <- function(GAlist=NULL, min_interactions=1, min_combined_score=900) {
 drugGenes <- readRDS(file=file.path(GAlist$dirs["gsets"],"drug2ensg.rds"))
 #file_string <- file.path(GAlist$dirs["gsets"],"9606.protein.links.v11.5.txt.gz")
 file_string <- file.path(GAlist$dirs["gsets"],"9606.protein.links.v12.0.txt.gz")
 string <- fread(file_string, data.table=FALSE)
 string  <- string[string$combined_score>=min_combined_score,]
 string <- split( string$protein2,f=as.factor(string$protein1))
 string <- string[sapply(string ,length)>=min_interactions]

 drug2ensp <- lapply(drugGenes,function(x){na.omit(unlist(GAlist$gsets$ensg2ensp[x]))})
 drugComplex <- lapply(drug2ensp,function(x){na.omit(unlist(string[x]))})
 drugComplex <- lapply(drugComplex,function(x){na.omit(unlist(GAlist$gsets$ensp2ensg[x]))})
 drugComplex <- lapply(drugComplex, function(x){unique(x)})

 for(i in 1:length(drugComplex)) {
  drugComplex[[i]] <- unique(c(drugGenes[[i]], drugComplex[[i]]))
 }
 return(drugComplex)
}

#' Retrieve Marker Sets from Database
#'
#' Retrieves marker sets based on a given feature and feature ID. Supports various features
#' such as Genes, Entrez Genes, Gene Symbol, Proteins, Pathways, GO, Protein Complexes,
#' Chemical Complexes, and others. Optionally, can filter by a list of SNP IDs.
#'
#' @param GAlist A list object with the gact database infrastructure.
#' @param feature A string specifying the biological feature type for the marker sets.
#' @param featureID ID of the feature to retrieve marker sets for.
#' @param rsids Optional vector of rsids to filter the marker sets.
#' @param threshold Threshold value for filtering (default is 0.01).
#'
#' @return A list of marker sets for the specified feature and feature ID, each set
#'         being a character vector of rsids.
#'
#' @examples
#' \dontrun{
#'   sets <- getMarkerSets(GAlist = GAlist, feature = "Genes",
#'                           featureID = c("ENSG00000165879", "ENSG00000012048"))
#' }
#' @importFrom msigdbr msigdbr
#' @export
getMarkerSets <- function(GAlist = NULL, feature = NULL, featureID = NULL,
                            rsids = NULL, threshold = 0.01) {
 if(is.null(feature)) stop("Feature is required")

 setsfile <- NULL
 if(feature=="Genesplus") setsfile <- file.path(GAlist$dirs["gsets"],"ensg2rsids.rds")
 if(feature=="Genes") setsfile <- file.path(GAlist$dirs["gsets"],"ensg2rsids_10kb.rds")
 if(feature=="Entrez Genes") setsfile <- file.path(GAlist$dirs["gsets"],"eg2rsids_10kb.rds")
 if(feature=="Gene Symbol") setsfile <- file.path(GAlist$dirs["gsets"],"sym2rsids_10kb.rds")
 if(feature=="Proteins") setsfile <- file.path(GAlist$dirs["gsets"],"ensp2rsids_10kb.rds")
 if(feature=="Pathways") setsfile <- file.path(GAlist$dirs["gsets"],"reactome2rsids.rds")
 if(feature=="GO") setsfile <- file.path(GAlist$dirs["gsets"],"go2rsids.rds")
 if(feature=="ProteinComplexes") setsfile <- file.path(GAlist$dirs["gsets"],"string2rsids.rds")
 if(feature=="ChemicalComplexes") setsfile <- file.path(GAlist$dirs["gsets"],"stitch2rsids.rds")
 if(feature=="DrugGenes") setsfile <- file.path(GAlist$dirs["gsets"],"drug2ensg.rds")
 if(feature=="DrugComplexes") setsfile <- file.path(GAlist$dirs["gsets"],"drug2string2ensg.rds")
 if(!is.null(setsfile)) sets <- readRDS(file=setsfile)

 if(feature=="Chromosomes") {
  markers <- fread(file.path(GAlist$dirs["marker"],"markers.txt.gz"),
                   data.table=FALSE)
  sets <- split( markers$rsids, f=as.factor(markers$chr) )
 }
 if(feature=="KEGG") {
  msets <- readRDS(file.path(GAlist$dirs["gsets"],"ensg2rsids.rds"))
  msigdb <- msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG")
  sets <- split(msigdb$ensembl_gene, f=msigdb$gs_name)
  sets <- mapSets(sets,names(msets), index=FALSE)
  sets <- lapply(sets, function(x) {unlist(msets[x])})
 }
 if(feature%in%c("GTEx","GTExV7","GTExV8")) {
  gtexsets <- getFeatureSets(GAlist=GAlist, feature=feature)
  sets <- sapply(gtexsets, function(x){ x$rsids[x$p<threshold]})
 }
 if(feature%in%c("Regulatory Categories","Regulatory Regions","Promoter", "Enhancer","OCR",
                 "TF", "CTCF")) {
  if(feature=="Regulatory Categories") setsfile <- file.path(GAlist$dirs["gsets"], "reg2rsids.rds")
  if(feature=="Regulatory Regions") setsfile <- file.path(GAlist$dirs["gsets"], "ensr2rsids.rds")

  sets <- readRDS(file=setsfile)
 }

 # Filter by featureID if provided
 if(!is.null(featureID)) {
  inSet <- featureID %in% names(sets)
  if(any(!inSet)) warning("Some IDs not in database:", paste(featureID[!inSet], collapse = ", "))
  sets <- sets[featureID]
 }

 # Filter by rsids if provided
 if(!is.null(rsids)) {
  sets <- qgg:::mapSets(sets = sets, rsids = rsids, index = FALSE)
 }

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


#' Retrieve Marker Statistics from Database
#'
#' This function fetches statistical data for markers from a genomic association database.
#' It allows filtering based on study IDs and various parameters. Options include fetching
#' specific statistics like effect sizes, standard errors, p-values, etc., for one or more studies.
#'
#' @param GAlist A list object providing information and infrastructure of the genomic association database.
#' @param studyID A vector of study IDs for which statistics are to be retrieved.
#' @param what Specifies the type of statistics to retrieve; "all" for all statistics or specific
#'        types like "rsids", "b", "seb", "eaf", "ea", "nea", "z", "p" (default: "all").
#' @param format Specifies the format of the output; either "list" or "data.frame" (default: "list").
#' @param rm.na Boolean indicating whether to remove NA values (default: TRUE).
#' @param rsids An optional vector of SNP IDs to subset the marker statistics.
#' @param cpra An optional vector of chromosome-position-ref-alt values to map to marker IDs.
#'
#' @return Depending on the 'format' and 'what' parameters, the function returns a list or
#'         data.frame of the requested marker statistics.
#'
#' @examples
#' \dontrun{
#'   stats <- getMarkerStat(GAlist = GAlist, studyID = "GWAS1")
#' }
#' @export
#'
getMarkerStat <- function(GAlist=NULL, studyID=NULL, what="all", format="list", rm.na=TRUE, rsids=NULL, cpra=NULL) {

 if(!is.null(cpra)) {
  markers <- fread(GAlist$markerfiles, data.table=FALSE)
  cpra1 <- paste(markers[,"chr"],
                 markers[,"pos"],
                 toupper(markers[,"ea"]),
                 toupper(markers[,"nea"]), sep="_")
  cpra2 <- paste(markers[,"chr"],
                 markers[,"pos"],
                 toupper(markers[,"nea"]),
                 toupper(markers[,"ea"]),sep="_")

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
   stat$p <- as.numeric(stat$p)
   stat$p[stat$p<.Machine$double.xmin] <- .Machine$double.xmin
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
  stat$p <- as.numeric(stat$p)
  stat$p[stat$p<.Machine$double.xmin] <- .Machine$double.xmin
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
   stat$p <- as.numeric(stat$p)
   stat$p[stat$p<.Machine$double.xmin] <- .Machine$double.xmin
   if(what=="z") res[stat$rsids,study] <- stat$b/stat$seb
   if(!what=="z") res[stat$rsids,study] <- stat[,what]
  }
  if(!is.null(rsids)) rsids <- rsids[rsids%in%GAlist$rsids]
  if(is.null(rsids)) rsids <- stat$rsids[stat$rsids%in%GAlist$rsids]
  rsids <- match(rsids,GAlist$rsids)

  if(rm.na) res <- na.omit(res[rsids,])
  if(!rm.na) res <- res[rsids,]
  if(what=="rsids") res <- res[,1]
  return(res)
 }

}

#' Retrieve GWAS Results
#'
#' This function is a wrapper around `getMarkerStat()` and inherits its documentation. It retrieves GWAS results for a specified study based on the given parameters. Additional parameters can be passed to filter the results.
#'
#' @inherit getMarkerStat params return
#'
#' @export
getGWAS <- function(GAlist=NULL, studyID=NULL, what="all", format="list", rm.na=TRUE, rsids=NULL, cpra=NULL) {

 getMarkerStat(GAlist=GAlist, studyID=studyID, what=what, format=format, rm.na=rm.na, rsids=rsids, cpra=cpra)
}

#' Retrieve VEP (Variant Effect Predictor) information for specified genes or rsids.
#'
#' This function fetches VEP annotations for either a list of genes (via `ensg`) or specific rsids.
#' It searches for the corresponding rsids if genes are provided and then extracts VEP data
#' from the Ensembl VEP files.
#'
#' @param GAlist A list object containing necessary data, including the directory (`dbdir`) where VEP files are located.
#' @param ensg A character vector of gene Ensembl IDs (e.g., ENSG00000139618) for which rsids should be retrieved and analyzed.
#' @param enst A character vector of transcript Ensembl IDs (optional, not used in the current function).
#' @param ensp A character vector of protein Ensembl IDs (optional, not used in the current function).
#' @param rsids A character vector of rsids (e.g., rs123) for which VEP annotations should be retrieved.
#'
#' @return A data frame containing the VEP annotations for the requested rsids, with columns:
#' \itemize{
#'   \item \code{rsids}: The rsid (variant ID).
#'   \item \code{Location}: The chromosomal location of the variant.
#'   \item \code{REF_ALLELE}: The reference allele.
#'   \item \code{Allele}: The alternate allele.
#'   \item \code{Consequence}: The predicted consequence of the variant.
#'   \item \code{IMPACT}: The predicted impact of the variant.
#'   \item \code{SYMBOL}: The gene symbol.
#'   \item \code{Gene}: The Ensembl gene ID.
#'   \item \code{Feature_type}: The type of feature (e.g., transcript).
#'   \item \code{Feature}: The Ensembl transcript ID.
#'   \item \code{BIOTYPE}: The biotype of the transcript (e.g., protein_coding).
#'   \item Additional columns for SIFT, PolyPhen, CLIN_SIG, and others.
#' }
#'
#' @details
#' The function first checks whether rsids are provided directly or need to be retrieved
#' based on the provided Ensembl gene IDs (`ensg`). It then reads chromosome-specific
#' VEP files from the provided `GAlist` directory and extracts relevant variant annotations
#' for the specified rsids. Results are filtered and returned as a data frame.
#'
#' @note The VEP files must be located in the `dbdir` of the `GAlist` object, and the files should be named
#' as `vep01.txt.gz`, `vep02.txt.gz`, ..., corresponding to each chromosome.
#'
#' @examples
#' # Example with Ensembl gene IDs
#' GAlist <- list(dbdir = "path_to_vep_files")
#' vep_data <- getVEP(GAlist, ensg = c("ENSG00000139618"))
#'
#' # Example with rsids
#' vep_data <- getVEP(GAlist, rsids = c("rs123", "rs456"))
#'
#' @export
getVEP <- function(GAlist=NULL, ensg=NULL, enst=NULL, ensp=NULL, rsids=NULL) {
 rs2vep_list <- list()  # Store data frames temporarily

 # Check if ensg or rsids are provided
 if(is.null(ensg) && is.null(rsids)) stop("Please provide either 'ensg' or 'rsids' to retrieve VEP data.")

 # Retrieve rsids for genes if ensg is provided
 if(!is.null(ensg)) {
  ensgSets <- getMarkerSets(GAlist=GAlist, feature="Genesplus")
  print(paste("Found information for", sum(ensg %in% names(ensgSets)), "genes in database"))
  rsids <- unlist(ensgSets[names(ensgSets) %in% ensg])
  rsids <- unique(rsids)
  if(length(rsids) == 0) stop("No rsids were found for the provided genes")
 }

 # Read chromosome-specific VEP files
 files <- paste0(file.path(GAlist$dbdir, "vep"), "/vep", sprintf("%02d", 1:22), ".txt.gz")
 for (i in seq_along(files)) {
  df <- fread(files[i], data.table = FALSE)
  rws <- df[,"#Uploaded_variation"] %in% rsids
  rs2vep_list[[i]] <- df[rws, ]
  print(paste("Finished processing chromosome", i))
 }

 # Combine all chromosome data into one data frame
 rs2vep <- do.call(rbind, rs2vep_list)

 # Define columns to retain
 cls <- c("#Uploaded_variation","Location","REF_ALLELE","Allele","Consequence","IMPACT",
          "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON",
          "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
          "DISTANCE", "SYMBOL_SOURCE", "CANONICAL", "SIFT", "PolyPhen", "CLIN_SIG",
          "LOEUF", "am_class", "am_pathogenicity", "MPC", "CADD_PHRED")

 # Check if the required columns are available and subset
 available_columns <- colnames(rs2vep)
 rs2vep <- rs2vep[, cls[cls %in% available_columns], drop = FALSE]

 # Rename the first column to "rsids"
 colnames(rs2vep)[1] <- "rsids"

 return(rs2vep)
}
# Some processing notes on how we obtained the VEP information from Ensembl
# path <- getwd()
# gact:::download_zenodo(doi = "10.5281/zenodo.10467174", path=path)
# markers <- fread(file.path(path,"markers.txt.gz"),data.table=FALSE)
#
# for(CHR in 1:22){
#  tmp <- markers[markers$chr==CHR,]
#  vcf <- tmp[,c("chr","pos","rsids","ea","nea")]
#  write.table(vcf, paste0(path,"/markers_chr",CHR,".vcf"), quote=FALSE, col.names=FALSE, row.names=FALSE)
# }
#
# # Upload filerne til VEP, GRCh37
# #https://grch37.ensembl.org/Homo_sapiens/Tools/VEP
# # Restrict output to MANE, i.e., MANE_SELECT IS NOT -
#
# # Specify the directory containing the .gz files
# gz_dir <- "path_to_your_gz_files"
#
# # List all the .gz files in the directory
# gz_files <- list.files(file.path(GAlist$dbdir, "vep"), pattern = "\\.gz$", full.names = TRUE)
#
#
# # Specify the output tarball file path
# tarfile <- file.path(GAlist$dbdir, "vep", "vep.tar.gz")
#
# # Combine the .gz files into a tarball
# tar(tarfile = tarfile, files = gz_files, compression = "gzip")



#' Get GTEx data
#'
#' This function retrieves gene expression data from the GTEx project,
#' based on the specified GTEx version, tissue type, and RSIDs (if provided).
#' It supports both GTEx V7 and V8 versions. GTEx downloaded from here:
#' https://www.gtexportal.org/home/downloads/adult-gtex/qtl
#'
#' @param GAlist A list containing the directory path with GTEx data.
#' @param version A character string specifying the GTEx data version,
#'        either "V8" or "V7".
#' @param rsids A vector of RSIDs for which data is to be retrieved.
#'        If NULL, data for all available RSIDs is retrieved.
#' @param format A character string specifying the output format,
#'        either "data.frame" or "list".
#' @param tissue A vector of tissue names for which data is to be retrieved.
#' @param tissue A vector of tissue names for which data is to be retrieved.
#'        Possible tissues include: "Adipose_Subcutaneous", "Adipose_Visceral_Omentum",
#'        "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial",
#'        "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia",
#'        "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex",
#'        "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus",
#'        "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia",
#'        "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", "Breast_Mammary_Tissue",
#'        "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", "Colon_Sigmoid",
#'        "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa",
#'        "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle",
#'        "Kidney_Cortex", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal",
#'        "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate",
#'        "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg",
#'        "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis",
#'        "Thyroid", "Uterus", "Vagina", "Whole_Blood".
#'        If NULL, data for all available tissues is retrieved.
#' @return A list or a data frame (based on the `format` parameter)
#'         containing gene expression data.
#' @examples
#' GAlist <- list(dbdir = "path/to/GTEx/data")
#' gtx_data <- getGTX(GAlist, version = "V8", tissue = c("Liver", "Heart"))
#' @export
getGTX <- function(GAlist = NULL, version = "V8", rsids = NULL, ensg=NULL,
                   tissue = NULL, cls=NULL, format = "data.frame") {

 # Initialize an empty list to store results
 gtx <- list()

 if(!is.null(ensg)) {
  ensgSets <- getMarkerSets(GAlist=GAlist, feature="Genesplus")
  print(paste("Found information for",sum(ensg%in%names(ensgSets)),"genes in database"))
  rsids <- unlist(ensgSets[names(ensgSets) %in% ensg])
  rsids <- unique(rsids)
  if(length(rsids) == 0) stop("No rsids were found for the provided genes")
 }


 # Define the directory suffix based on the version
 dbdir_suffix <- if (version == "V8") {
  "gtex/GTEx_Analysis_v8_eQTL"
 } else if (version == "V7") {
  "gtex/GTEx_Analysis_v7_eQTL"
 } else {
  stop("Invalid GTEx version. Only 'V7' and 'V8' are supported.")
 }

 # Construct the directory path
 dbdir <- file.path(GAlist$dbdir, dbdir_suffix)

 # Get the list of files
 files <- list.files(dbdir, pattern = "egenes", full.names = TRUE)

 # Extract tissue names from filenames
 files_tissue <- gsub("\\.v[78]\\.egenes\\.txt\\.gz$", "", basename(files))

 # Filter files by tissue if specified
 if (!is.null(tissue)) {
  select <- files_tissue %in% tissue
  select <- files_tissue %in% tissue
  if (!any(select)) stop("Some of the selected tissues not in this GTEx version")
  files <- files[select]
  files_tissue <- files_tissue[select]
 }

 # Process each file
 for (i in seq_along(files)) {
  df <- fread(files[i], data.table = FALSE)
  if(i==1) dfcls <- colnames(df)
  if(i>1) {
   if(!any(colnames(df)%in%dfcls)) stop("Mismatch between column names in GTEx files")
  }
  df$gene_id <- substr(df$gene_id, 1, 15)
  df$chr <- gsub("chr","",df$chr)
  df$variant_id <- gsub("chr|_b37|_b38", "", df$variant_id)
  rsid_column <- if (version == "V7") "rs_id_dbSNP147_GRCh37p13" else "rs_id_dbSNP151_GRCh38p7"
  if (is.null(cls)) {
   gtx[[i]] <- data.frame(tissue = files_tissue[i],
                          ensg = df$gene_id,
                          rsids = df[, rsid_column],
                          p = df$pval_beta,
                          q = df$qval)
  }
  if (is.null(cls)) {
   rws <- 1:nrow(df)
   if(!is.null(rsids)) rws <- df[, rsid_column]%in%rsids
   gtx[[i]] <- cbind(tissue=files_tissue[i], df[rws,])
   colnames(gtx[[i]])[colnames(gtx[[i]])==rsid_column] <- "rsids"
   colnames(gtx[[i]])[colnames(gtx[[i]])=="gene_id"] <- "ensg"

  }

  # gtx[[i]] <- data.frame(tissue = tissue[i],
  #                        ensg = df$gene_id,
  #                        rsids = df[, rsid_column],
  #                        chr = df$chr,
  #                        ref = df$ref,
  #                        alt = df$alt,
  #                        maf = df$maf,
  #                        tss_distance= df$tss_distance,
  #                        qval = df$qval)
  message(paste("Processing file:", basename(files[i])))
 }

 # Combine all data frames if required
 if (format == "data.frame") {
  gtx <- do.call(rbind, gtx)
  if(!is.null(ensg)) gtx <- gtx[gtx$ensg%in%ensg,]
  if(!is.null(rsids)) gtx <- gtx[gtx$rsids%in%rsids,]
  gtx <- cbind(gtx[,"rsids"], gtx[, colnames(gtx) != "rsids"])
  colnames(gtx)[1] <- "rsids"
 }

 return(gtx)
}


#' Get VEGAS Data
#'
#' This function retrieves VEGAS (Versatile Gene-based Association Study) Z-score data from the GAlist object, allowing optional filtering by `studyID` and handling of missing values. The function returns a filtered matrix of Z-scores corresponding to the specified `studyID` or all data if no specific study is provided.
#'
#' @param GAlist A list object that contains directories and data required for the analysis, including the path to the Z-score data.
#' @param ensg Not currently used in the function, but can be reserved for future functionality.
#' @param studyID A character vector of study IDs to filter the Z-scores. Only Z-scores corresponding to these study IDs will be returned. If `NULL`, all studies are returned.
#' @param rm.na Logical; if `TRUE`, rows with `NA` values will be removed from the data. Defaults to `TRUE`.
#'
#' @return A matrix of Z-scores filtered by the provided `studyID` (if given) and cleaned of missing values if `rm.na` is `TRUE`.
#' @examples
#' # Assuming GAlist is a valid object with the necessary structure
#' Z_scores <- getVEGAS(GAlist, studyID = c("study1", "study2"))
#'
#' # Get all Z-scores and remove rows with NA values
#' Z_scores_all <- getVEGAS(GAlist, rm.na = TRUE)
#'
#' @export
getVEGAS <- function(GAlist = NULL, ensg = NULL, studyID = NULL, rm.na = TRUE) {
 # Read the Z-score matrix from the specified file path in the GAlist object
 Z <- readRDS(file = file.path(GAlist$dirs["gsea"], "Z_vegas.rds"))

 # If studyID is provided, filter the Z-scores matrix by matching column names
 if (!is.null(studyID)) {
  cls <- colnames(Z) %in% studyID
  if (any(!studyID %in% colnames(Z))) {
   warning("Some studyIDs not present in the database")
  }
  Z <- Z[, cls, drop = FALSE]  # Drop any unwanted columns
 }

 # If rm.na is TRUE, remove rows with NA values
 if (rm.na) Z <- na.omit(Z)

 # Return the filtered Z-scores matrix
 Z
}


#' Retrieve Disease Associations
#'
#' This function retrieves disease association data based on protein or gene information from the provided GAlist object. It supports filtering by `ENSP` (protein) or `ENSG` (gene), and can fetch data from multiple sources: Integrated, Experiments, Knowledge, or Textmining datasets.
#'
#' @param GAlist A list object that contains directories and data required for the analysis, including disease data files and gene-protein mappings.
#' @param ensg A character vector of gene IDs (ENSG) to filter the disease data. If not provided, the function will use `ensp` to filter.
#' @param ensp A character vector of protein IDs (ENSP) to filter the disease data.
#' @param what A character string specifying the source of disease data. Options are "Integrated" (default), "Experiments", "Knowledge", or "Textmining".
#'
#' @return A data frame of disease associations filtered by the provided gene or protein IDs, sorted by the disease score in decreasing order.
#' @examples
#' # Assuming GAlist is a valid object with the necessary structure
#' diseases <- getDISEASES(GAlist, ensg = c("ENSG00000139618"))
#'
#' # Retrieve data using protein IDs
#' diseases_protein <- getDISEASES(GAlist, ensp = c("ENSP00000269305"))
#'
#' @export
getDISEASES <- function(GAlist = NULL, ensg = NULL, ensp = NULL, what = "Integrated") {

 # Determine the file to load based on the 'what' parameter
 if (what == "Integrated") file <- "human_disease_integrated_full.tsv.gz"
 if (what == "Experiments") file <- "human_disease_experiments_filtered.tsv.gz"
 if (what == "Knowledge") file <- "human_disease_knowledge_filtered.tsv.gz"
 if (what == "Textmining") file <- "human_disease_textmining_filtered.tsv.gz"

 # Read the disease data file
 df <- fread(file.path(GAlist$dirs["gsets"], file), data.table = FALSE)
 colnames(df) <- c("ensp", "sym", "DOID", "Disease", "Score")  # Ensure proper column names

 # Filter by protein IDs (ENSP), if provided
 if (!is.null(ensp)) df <- df[df[,1] %in% ensp, ]

 # Filter by gene IDs (ENSG), if provided
 if (!is.null(ensg)) {
  # Get the ENSP entries for the given ENSG IDs
  ensp <- GAlist$gsets$ensg2ensp[names(GAlist$gsets$ensg2ensp) %in% ensg]
  ensg2ensp <- do.call(rbind, lapply(names(ensp), function(ensg) {
   data.frame(ensg = ensg, ensp = ensp[[ensg]], stringsAsFactors = FALSE)
  }))

  # Filter the disease data by matching ENSP values
  df <- df[df$ensp %in% unique(unlist(ensp)), ]

  # Merge the gene-protein mapping with the disease data
  df <- merge(x = ensg2ensp, y = df, by = "ensp", all = FALSE)
 }

 # Sort the data by the "Score" column in descending order
 o <- order(df$Score, decreasing = TRUE)
 df[o, ]  # Return the sorted data frame
}

#' @export
#'

# hyperg test
hgtestDB <- function(p = NULL, sets = NULL, threshold = 0.05) {
 population_size <- length(p)
 sample_size <- sapply(sets, length)
 n_successes_population <- sum(p < threshold)
 n_successes_sample <- sapply(sets, function(x) {
  sum(p[x] < threshold)
 })
 phyperg <- rep(1,length(sets))
 names(phyperg) <- names(sets)
 for (i in 1:length(sets)) {
  phyperg[i] <- 1.0-phyper(n_successes_sample[i]-1, n_successes_population,
                           population_size-n_successes_population,
                           sample_size[i])
 }
 # Calculate enrichment factor
 ef <- (n_successes_sample/sample_size)/
  (n_successes_population/population_size)

 # Create data frame for table
 df <- data.frame(feature = names(sets),
                  ng = sample_size,
                  nag = n_successes_sample,
                  ef=ef,
                  phgt = phyperg)
 colnames(df) <- c("Feature", "Number of Genes",
                   "Number of Associated Genes",
                   "Enrichment Factor",
                   "P-value")
 df
}




#' @export
#'
# Create a plot function
mhplotDB <- function(p=NULL, main=NULL) {
 pobs <- -log10(p)
 plot(pobs,
      pch = 20, ylim = c(0, max(pobs)), main=main,
      xlab = "Position", ylab = "-log10(p-value)", frame.plot=FALSE)
 abline(h = 8, lty = 2, col = "red")
}

#' @export
#'
# Create a plot function
qqplotDB <- function(p=NULL, main=NULL) {
 pobs <- -log10(p)
 pexp <- -log10((1:length(pobs))/length(pobs) )
 plot( y=sort(pobs), x=sort(pexp),
       xlab="Expected -log10(P)",
       ylab="Observed -log10(P)",
       frame.plot=FALSE, main=main)
 abline(a=0,b=1, col="red", lty=2, lwd=2)
}




#' @export
#'
getStudiesShinyDB <- function(GAlist=NULL) {
 df_studies <- as.data.frame(GAlist$study)[,-c(2,11)]
 df_studies$neff <- as.integer(df_studies$neff)
 df_studies$reference <- createURL(url="https://pubmed.ncbi.nlm.nih.gov/",
                                   urlid=gsub("PMID:","",df_studies$reference))
 df_studies
}


#' @export
#'
# Create a plot function
createURL <- function(url=NULL,urlid=NULL){
 url <- paste0(url, urlid)
 html_code <- paste0("<a href='", url, "' target='_blank'>", urlid, "</a>")
 html_code
}

#' @export
#'
addATC <- function(drugname=NULL) {
 atc <- rep("Unknown",length(drugname))
 has_atc <- match(tolower(drugname),tolower(GAlist$atc$name))
 atc[!is.na(has_atc)] <- as.character(GAlist$atc$code[has_atc[!is.na(has_atc)]])
 return(atc)
}

#' @export
#'
createDT <- function(df=NULL) {
 dt <- datatable(df, extensions = "Buttons",
                 escape = FALSE,
                 options = list(paging = TRUE,
                                scrollX = TRUE,
                                searching = TRUE,
                                ordering = TRUE,
                                dom = 'Bfrtip',
                                buttons = c('csv', 'excel'),
                                pageLength = 10,
                                lengthMenu = c(3, 5, 10)),
                 rownames = FALSE)
 return(dt)
}


# Extract gene information from Ensembl
#' @export
#'
getGeneDB <- function(symbol=NULL) {
 base_url <- "https://rest.ensembl.org"
 url <- paste0(base_url, "/lookup/symbol/homo_sapiens/", symbol)
 response <- GET(url = url, content_type("application/json"))
 gene <- fromJSON(content(response, "text"), flatten = TRUE)
 if(!is.null(gene$error)) gene_information <- c(symbol, rep(NA, 6))
 if(is.null(gene$error)) gene_information <- c(symbol, gene$id, gene$description, gene$biotype, gene$seq_region_name, gene$start, gene$end)
 return(gene_information)
}

# Define function to retrieve interaction partners
#' @export
#'
getInteractionsDB <- function(ids=NULL, species="9606", threshold=900) {
 ids <- paste0(ids, collapse = "%0d")
 url <- paste0("https://string-db.org/api/tsv/interaction_partners?identifiers=",
               ids,
               "&species=",
               species,
               "&required_score=",
               threshold)
 res <- GET(url)
 interactions <- read.table(text = content(res, "text"), sep = "\t",
                            stringsAsFactors = FALSE, header=TRUE)
 return(interactions)
}


#' @export
mapSetsDB <- function(sets = NULL, featureID = NULL, GAlist = NULL, index = TRUE) {
 nsets <- sapply(sets, length)
 if(is.null(names(sets))) names(sets) <- paste0("Set",1:length(sets))
 rs <- rep(names(sets), times = nsets)
 rsSets <- unlist(sets, use.names = FALSE)
 rsSets <- match(rsSets, featureID)
 inW <- !is.na(rsSets)
 rsSets <- rsSets[inW]
 if (!index) rsSets <- featureID[rsSets]
 rs <- rs[inW]
 rs <- factor(rs, levels = unique(rs))
 rsSets <- split(rsSets, f = rs)
 return(rsSets)
}

#' Retrieve Linkage Disequilibrium (LD) Scores for Specific Genetic Markers
#'
#' This function retrieves LD scores for specified genetic markers based on the population ancestry and version of the reference panel. The function reads in marker files and extracts LD scores, allowing for optional filtering by specific rsIDs.
#'
#' @param GAlist A list containing directory paths where marker data is stored. It is assumed that `GAlist$dirs["marker"]` points to the correct directory.
#' @param ancestry A character string specifying the ancestry population. The options are "EUR" (European), "EAS" (East Asian), or "SAS" (South Asian). Default is "EUR".
#' @param version A character string specifying the version of the reference panel. The options are "HapMap3", "1000G", or "Original". Default is "HapMap3".
#' @param rsids An optional character vector of rsIDs to filter the LD scores. If not specified, all LD scores will be returned.
#'
#' @return A named vector of LD scores, where the names are rsIDs.
#'
#' @details
#' The function reads in precomputed LD score files for the selected ancestry and version, extracts the LD scores, and returns them as a named vector. If the `rsids` argument is provided, the returned vector will only include LD scores for the specified rsIDs.
#'
#' @examples
#' # Example usage:
#' GAlist <- list(dirs = c(marker = "/path/to/marker/files"))
#' ldscores <- getLDscores(GAlist, ancestry = "EUR", version = "HapMap3", rsids = c("rs123", "rs456"))
#'
#' @export
#' @export
getLDscores <- function(GAlist=NULL, ancestry="EUR", version="HapMap3", rsids=NULL) {
 # Check if GAlist is provided
 if (is.null(GAlist) || !("marker" %in% names(GAlist$dirs))) {
  stop("Error: GAlist or GAlist$dirs['marker'] is missing.")
 }

 # Check if ancestry and version are valid
 valid_ancestries <- c("EUR", "EAS", "SAS")
 valid_versions <- c("HapMap3", "1000G", "Original")

 if (!(ancestry %in% valid_ancestries)) {
  stop("Error: Invalid ancestry. Choose from 'EUR', 'EAS', or 'SAS'.")
 }

 if (!(version %in% valid_versions)) {
  stop("Error: Invalid version. Choose from 'HapMap3', '1000G', or 'Original'.")
 }

 # Define the file path based on ancestry and version
 file_name <- switch(
  paste(ancestry, version, sep = "_"),
  "EUR_HapMap3" = "markers_1000G_eur_hm3.txt.gz",
  "EAS_HapMap3" = "markers_1000G_eas_hm3.txt.gz",
  "SAS_HapMap3" = "markers_1000G_sas_hm3.txt.gz",
  "EUR_1000G" = "markers_1000G_eur_filtered.txt.gz",
  "EAS_1000G" = "markers_1000G_eas_filtered.txt.gz",
  "SAS_1000G" = "markers_1000G_sas_filtered.txt.gz",
  "EUR_Original" = "markers_1000G_eur_w_ld.txt.txt.gz",
  stop("Error: No marker file found for the selected ancestry and version.")
 )

 # Read the marker file
 marker_file_path <- file.path(GAlist$dirs["marker"], file_name)

 if (!file.exists(marker_file_path)) {
  stop("Error: The marker file does not exist at the specified location.")
 }

 marker <- fread(marker_file_path, data.table = FALSE)

 # Check if 'ldscores' and 'rsids' columns exist
 if (!all(c("ldscores", "rsids") %in% colnames(marker))) {
  stop("Error: The marker file does not contain the required columns 'ldscores' and 'rsids'.")
 }

 # Extract ldscores and assign rsids as names
 ldscores <- marker[,"ldscores"]
 names(ldscores) <- marker$rsids

 # Filter by rsids if provided
 if (!is.null(rsids)) {
  ldscores <- ldscores[names(ldscores) %in% rsids]
 }

 return(ldscores)
}


# Define the function to translate gene symbols to Ensembl IDs
#' @export
getGeneInfo <- function(sym=NULL,ensg=NULL) {
 if(!is.null(sym)) {
  result <- as.data.frame(t(sapply(sym, function(x){sym2info(x)})))
 }
 if(!is.null(ensg)) {
  result <- as.data.frame(t(sapply(ensg, function(x){ensg2info(x)})))
 }
 colnames(result) <- c("Symbol", "Ensembl Gene Id", "Description", "Type", "Chr", "Start", "End")
 return(result)
}

sym2info <- function(symbol) {
 base_url <- "https://rest.ensembl.org"
 url <- paste0(base_url, "/lookup/symbol/homo_sapiens/", symbol)
 response <- GET(url = url, content_type("application/json"))
 gene <- fromJSON(content(response, "text"), flatten = TRUE)
 if(!is.null(gene$error)) gene_information <- c(symbol, rep(NA, 6))
 if(is.null(gene$error)) gene_information <- c(symbol, gene$id, gene$description, gene$biotype, gene$seq_region_name, gene$start, gene$end)
 return(gene_information)
}

ensg2info <- function(ensg) {
 base_url <- "https://rest.ensembl.org"
 url <- paste0(base_url, "/lookup/id/", ensg)
 response <- GET(url = url, content_type("application/json"))
 gene <- fromJSON(content(response, "text"), flatten = TRUE)
 if(!is.null(gene$error)) gene_information <- c(symbol, rep(NA, 6))
 if(is.null(gene$error)) gene_information <- c(gene$display_name, gene$id, gene$description, gene$biotype, gene$seq_region_name, gene$start, gene$end)
 return(gene_information)
}

#' Add Annotations to a Data Frame
#'
#' This function enriches a data frame with annotations based on Ensembl gene or regulatory feature IDs.
#' It supports the addition of gene or regulatory annotations and can create hyperlinks for Excel.
#'
#' @param df Data frame to which annotations will be added.
#' @param ensg Vector of Ensembl gene IDs; used if 'feature' is set to "Genes".
#' @param ensr Vector of Ensembl regulatory feature IDs; used if 'feature' is set to "Regulatory".
#' @param feature Type of feature for annotation, either "Genes" or "Regulatory" (default: "Genes").
#' @param hyperlinkEXCEL Boolean, if TRUE, adds hyperlinks to the Ensembl website for each feature ID.
#'
#' @return An annotated data frame with additional columns for Ensembl IDs, symbols, and genomic coordinates.
#'         If 'hyperlinkEXCEL' is TRUE, the data frame will include Excel hyperlinks.
#' @export
#'
addAnnotationDB <- function(df = NULL, ensg = NULL, ensr = NULL, feature="Genes",
                            hyperlinkEXCEL = FALSE) {

 # Load annotation data
 if(feature=="Genes") annotation <- readRDS(file.path(GAlist$dirs["gsets"], "genesplus_annotation.rds"))
 if(feature=="Regulatory") annotation <- readRDS(file.path(GAlist$dirs["gsets"], "regulatory_annotation.rds"))

 # Process data frame
 if (!is.null(df)) {
  if (is.matrix(df)) df <- as.data.frame(df)
  ensid <- rownames(df)
  if (is.null(ensid)) stop("Please provide rownames for df (e.g. Ensembl Gene or Regulatory IDs)")
 }

 if(feature=="Regulatory") {
  ensr <- ensid
  if(any(ensr%in%annotation$reg_id)) {
   df <- cbind(annotation[ensr, ], df)
  }
  # Add hyperlinks if required
  if (hyperlinkEXCEL) {
   warning("Not possible to add hyperlink for regulatory features")
  }
 }
 if(feature=="Genes") {
  ensg <- ensid
  # Merge data with annotations
  df <- cbind(ensg, GAlist$gsets[["ensg2sym"]][ensg], annotation[ensg, -1], df)
  colnames(df)[1:5] <- c("Ensembl Gene ID", "Symbol", "Chr", "Start", "Stop")

  # Split and reformat symbols
  ensg2sym_list <- lapply(df[,"Symbol"], function(x) unlist(strsplit(x, split = " ")))
  gsym <- unlist(ensg2sym_list)
  rws <- rep(1:nrow(df), times = sapply(ensg2sym_list, length))
  df <- df[rws, ]
  df[,"Symbol"] <- gsym

  # Add hyperlinks if required
  if (hyperlinkEXCEL) {
   df <- addHyperlinks(df)
  }
 }

 return(df)
}

# Function to add hyperlinks for Excel
addHyperlinks <- function(df) {
 df_with_links <- cbind(df[, 1], df[, 1], df)
 colnames(df_with_links)[1:3] <- c("Ensembl", "Open Target", "Ensembl Gene ID")
 df_with_links[, "Ensembl"] <- createHyperlink(df_with_links[, "Ensembl"], "http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=")
 df_with_links[, "Open Target"] <- createHyperlink(df_with_links[, "Open Target"], "https://platform.opentargets.org/target/")
 return(df_with_links)
}

# Function to create hyperlink string
createHyperlink <- function(ids, base_url) {
 paste0("=Hyperlink(",'"',paste0(base_url,ids),'"',";",'"',ids,'"',")")
 #paste0("=HYPERLINK(\"", base_url, ids, "\", \"", ids, "\")")
}


# QQ-plot
#' @export
qplot <- function(p=NULL, main = "") {
 mlogObs <- -log10(p)
 m <- length(mlogObs)
 o <- order(mlogObs, decreasing = TRUE)
 mlogExp <- -log10((1:m) / m)
 plot( y = mlogObs[o], x = mlogExp, col = 2, pch = "+",
       frame.plot = FALSE, main = main, xlab = "Expected -log10(p)", ylab = "Observed -log10(p)")
 abline(a = 0, b = 1)
}

#' Hypergeometric Test on Gene Sets with Customizable Output
#'
#' This function performs a hypergeometric test based on the provided parameters.
#' It retrieves feature sets from a database, constructs a design matrix, maps these sets,
#' and then applies a hypergeometric test. The output of the test can be customized
#' using the `output` argument.
#'
#' @param GAlist A list containing gene annotation data. This parameter cannot be NULL.
#' @param sets A list or other structure containing sets of data for analysis.
#'        This parameter cannot be NULL.
#' @param feature A character string specifying the type of feature to be used in the analysis.
#' @param featureIDs An optional vector of feature IDs to filter the analysis.
#' @param minsets A numeric value specifying the minimum number of sets to be considered.
#'        Must be a positive integer. Defaults to 5.
#' @param output A character string specifying the format of the output.
#'        Options include "p" for returning only p-values, "summary" for a summary
#'        including both p-values and effect sizes, or NULL for the default full output.
#'
#' @return Depending on the `output` parameter, returns either a matrix of hypergeometric
#'         test results, a matrix of p-values, or a list with p-values and effect sizes.
#'
#' @examples
#' \dontrun{
#'   # Example usage for default full output
#'   hgt_results <- hgtSets(GAlist=myGAlist, sets=mySets, feature="myFeature")
#'   # Example usage for p-values only
#'   p_values <- hgtSets(GAlist=myGAlist, sets=mySets, feature="myFeature", output="p")
#'   # Example usage for summary output
#'   summary_results <- hgtSets(GAlist=myGAlist, sets=mySets, feature="myFeature", output="summary")
#' }
#'
#' @export
hgtSets <- function(GAlist = NULL, sets = NULL, feature = NULL, featureIDs = NULL,
                  minsets = 5, output = NULL) {

 # Check for valid input
 if (!is.null(minsets) && (!is.numeric(minsets) || minsets < 1)) {
  stop("minsets must be a positive integer")
 }
 # Check for valid 'GAlist'
 if (is.null(GAlist)) stop("'GAlist' cannot be NULL")
 # Check for valid 'sets'
 if (is.null(sets)) stop("'sets' cannot be NULL")

 # Get feature sets from the database
 featureSets <- getFeatureSets(GAlist = GAlist, feature = feature, minsets = minsets)
 rowids <- unique(unlist(featureSets))

 # Generate the design matrix
 X <- designMatrix(sets = featureSets, rowids = rowids)

 # Map sets from the database
 sets <- mapSetsDB(sets = sets, featureID = rownames(X), index = TRUE)

 # Apply hypergeometric test
 hgtResults <- apply(X, 2, function(x) {
  # Using the exported function from the qgg package for hypergeometric testing
  hgtestDB(p = 1 - x, sets = sets, threshold = 0.5)
 })
 p <- sapply(hgtResults, function(x) {x[,"P-value"]})
 ef <- sapply(hgtResults, function(x) {x[,"Enrichment Factor"]})

 rownames(p) <- rownames(ef) <- names(sets)
 if(output=="p") hgtResults <- p
 if(output=="summary") hgtResults <- list(p=p,ef=ef)

 return(hgtResults)
}



#' # Add annotation to data frame based on Ensemlb IDs
#' addAnnotationDB <- function(df=NULL, ensg=NULL, ensr=NULL, hyperlinkEXCEL=FALSE) {
#'  annotation <- readRDS(file = file.path(GAlist$dirs["gsets"], "genesplus_annotation.rds"))
#'  if(!is.null(df)) {
#'   if(is.matrix(df)) df <- as.data.frame(df)
#'   ensg <- rownames(df)
#'  }
#'  if(!is.null(ensg)) {
#'   df <- as.data.frame(ensg)
#'   rownames(df) <- ensg
#'  }
#'
#'  if(is.null(ensg)) stop("Please provide rownames for df (should be EnSembl Ids")
#'  df <- cbind(ensg,GAlist$gsets[["ensg2sym"]][ensg], annotation[ensg,-1], df)
#'  colnames(df)[1:5] <- c("Ensembl Gene ID","Symbol", "Chr", "Start", "Stop")
#'  ensg2sym_list <- lapply(df[,"Symbol"], function(x){
#'   unlist(strsplit(x,split=" "))})
#'  gsym <- unlist(ensg2sym_list)
#'  rws <- rep(1:nrow(df),times=sapply(ensg2sym_list,length))
#'  df <- df[rws,]
#'  df[,"Symbol"] <- gsym
#'  if(hyperlinkEXCEL) {
#'   df <- cbind(df[,1],df[,1], df)
#'   res2hyperlink_ensembl <- paste0("http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",df[,1])
#'   res2hyperlink_opentarget <- paste0("https://platform.opentargets.org/target/",df[,1])
#'   res2hyperlink_ensembl <- paste0("=Hyperlink(",'"',res2hyperlink_ensembl,'"',";",'"',df[,1],'"',")")
#'   res2hyperlink_opentarget <- paste0("=Hyperlink(",'"',res2hyperlink_opentarget,'"',";",'"',df[,1],'"',")")
#'   df[,1] <- res2hyperlink_ensembl
#'   df[,2] <- res2hyperlink_opentarget
#'   colnames(df)[1:3] <- c("Ensembl","Open Target","Ensembl Gene ID")
#'  }
#'  return(df)
#' }
