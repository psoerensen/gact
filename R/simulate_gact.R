#' Simulate gene sets for pathway analysis
#'
#' Generates gene sets for pathway analysis based on a given set of SNPs, using a
#' pathway database contained in a GAlist object. The function
#' samples causal SNPs and constructs gene sets of a given size and number of
#' causal genes.
#'
#' @param GAlist A list object containing SNP and pathway data.
#' @param mcausal Integer, number of causal SNPs to sample (default: 500).
#' @param causal Vector of SNP IDs for causal SNPs (default: NULL, SNPs are randomly sampled).
#' @param ngenes Integer or vector, number of genes per gene set (default: 10).
#' @param ncgenes Integer or vector, number of causal genes per gene set (default: 2).
#' @param nreps Integer, number of gene sets to generate (default: 10).
#' @param format String, format of the gene sets (default: "markers").
#' @return List of gene sets, each represented as a vector of SNP IDs.
#' @export
#' @examples
#' \dontrun{
#'   # Generate gene sets
#'   sets <- simPathwayDB(GAlist = GAlist, ngenes = 10, ncgenes = 2, nreps = 5)
#' }

simPathwayDB <- function(GAlist = NULL, mcausal = 500, causal = NULL, ngenes = 10,
                         ncgenes = 2, nreps = 10, format = "markers") {
 # Input validation
 if (is.null(GAlist)) stop("GAlist is required")
 if (min(ngenes) - max(ncgenes) < 0)
  stop("Number of causal genes (ncgenes) larger than number of genes in pathway (ngenes)")

 # Process marker sets
 rsids <- unlist(GAlist$rsids)
 causal <- if (is.null(causal)) sample(rsids, size = mcausal) else causal
 sets <- getMarkerSetsDB(GAlist, feature = "Genes")
 causal_genes <- names(Filter(function(x) any(x %in% causal), sets))
 non_causal_genes <- names(Filter(function(x) !any(x %in% causal), sets))

 # Generate gene sets
 return(generateGeneSets(causal_genes, non_causal_genes, ngenes, ncgenes, nreps, sets, format))
}

generateGeneSets <- function(causal_genes, non_causal_genes, ngenes, ncgenes, nreps, sets, format) {
 gsets <- list()
 setnames <- vector(mode = "character", length = nreps * length(ngenes) * length(ncgenes))
 ijk <- 0
 for (i in seq_len(nreps)) {
  for (j in seq_along(ngenes)) {
   for (k in seq_along(ncgenes)) {
    ijk <- ijk + 1
    csets <- sample(causal_genes, ncgenes[k])
    ncsets <- sample(non_causal_genes, ngenes[j] - ncgenes[k])
    gsets[[ijk]] <- if (format == "markers") unlist(sets[c(csets, ncsets)]) else c(csets, ncsets)
    setnames[ijk] <- paste(c(i, ngenes[j], ncgenes[k]), collapse = "_")
   }
  }
 }
 names(gsets) <- setnames
 gsets
}
