#######################################################################################
# GACT database simulation functions
#######################################################################################
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
#' @param nreps an integer indicating the number of gene sets to generate.
#'   Default is 10.
#'
#' @return a list of gene sets, where each gene set is represented as a vector of
#'   SNP IDs.
#'
#'
#' @examples
#' \dontrun{
#' # Generate gene sets
#' sets <- gpath(GAlist = GAlist, ngenes = 10, ncgenes = 2, nreps = 5)
#' }
#'
#' @export
#'
gpath <- function(GAlist = NULL,
                  mcausal=500, causal=NULL, overlap=FALSE,
                  ngenes=10, ncgenes=2, nreps=10, format="markers"){
 # mcausal is the number of causal markers
 # causal is a vector of causal markers
 # ngenes is the number (or vector) of genes in pathway
 # ncgenes is the number (or vector) of causal genes in pathway
 # nreps is the number of pathways

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
 for (i in 1:nreps) {
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

