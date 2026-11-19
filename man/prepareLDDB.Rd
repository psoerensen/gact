\name{prepareLDDB}
\alias{prepareLDDB}
\title{Prepare sparse LD and reference LD scores together}
\usage{
prepareLDDB(GAlist, ancestry = "EUR", nthreads = 1L,
  overwrite = FALSE, progress = TRUE)
}
\arguments{
\item{GAlist}{A database descriptor with a build, marker directory and rsIDs.}
\item{ancestry}{One or more of EUR, AFR, EAS, SAS or AMR.}
\item{nthreads}{Number of threads passed to gbase.}
\item{overwrite}{Replace existing outputs instead of checking reuse.}
\item{progress}{Report preparation, computation and reuse stages.}
}
\value{Updated GAlist with sparse resources in \code{ldPanels} and score tables
in \code{ldscorePanels}. GAlist is not saved automatically.}
\description{Uses the existing gsuite sparse-LD workflow on the database's
matching ancestry-specific BED/BIM/FAM panels. Scores are extracted from the
sparse resource produced by the same computation. Numerical computation
remains in gbase/gmat.}
\details{The settings match \code{\link{prepareLDscoresDB}}: a
1,000-selected-marker distance window within chromosomes, no physical window
or r-squared cutoff, no extra MAF filter, and uncorrected scores including the
unit diagonal. Signed correlations are stored in gbase's disk-backed CSR format.
Reference panels already aligned to the database backbone are used directly;
coordinates and BIM A1/A2 are checked against catalog EA/NEA. Input genotypes
are never modified.

Files use the existing marker layout:
\code{g1000_<ancestry>.subset.aligned.*},
\code{Glist_1000G_<ancestry>_filtered.rds}, and the score tables read by
\code{\link{getLDscoresDB}}. The checksum receipt covers the sparse resource,
its sidecars and descriptor. The returned GAlist registers sparse panels for
\code{\link{writeMarkerStat}}. Repeated calls verify and reuse unchanged
results. Verified scores-only tables can be extended with sparse LD without
\code{overwrite = TRUE}. Unknown or incomplete sparse outputs require explicit
overwrite. Sparse storage can be substantial with the zero r-squared cutoff.
Temporary alignment tables are released before genotype preparation;
the saved Glist is released before assembling the score output table.
The existing gsuite mapped reader reuses the immutable metadata index under
\code{marker/g1000_<ancestry>.metadata-<checksum>/}, avoiding the resident
Glist summary columns and native resident-panel reconstruction for sparse LD.}
\seealso{\code{\link{getLDscoresDB}}, \code{\link{writeMarkerStat}}}
