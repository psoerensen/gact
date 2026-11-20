\name{prepareLDSC}
\alias{prepareLDSC}
\title{Prepare a HapMap3 reference for ordinary LDSC}
\usage{
prepareLDSC(GAlist, ancestry = "EUR", nthreads = 1L, hm3_rsids = NULL,
            overwrite = FALSE)
}
\arguments{
\item{GAlist}{Database descriptor with aligned 1000G BED genotypes.}
\item{ancestry}{One of EUR, AFR, EAS, SAS or AMR.}
\item{nthreads}{Number of genotype computation threads.}
\item{hm3_rsids}{Optional HapMap3 marker IDs. The default uses rsIDs from the
checksum-pinned Pan-UKB EUR HapMap3 resource as a common selection list.}
\item{overwrite}{Recompute a previously prepared reference.}
}
\description{Computes two finite-sample corrected sums of squared correlations
within 1 cM on each ancestry panel. Regression scores include variable reference
SNPs with MAF above 1 percent. Weight scores include the HapMap3 regression
SNPs. Regression excludes the extended MHC (chr6, 25--34 Mb inclusive).
The heritability multiplier counts reference SNPs with MAF above 5 percent.}
\value{Updated GAlist with the reference path in ldscPanels. The file in marker/
stores both score columns, the common-SNP count, settings and input checksums.}
\seealso{\code{\link{prepareLD}}, \code{\link{getLDscores}}}
