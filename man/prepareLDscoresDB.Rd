\name{prepareLDscoresDB}
\alias{prepareLDscoresDB}
\title{Prepare reference LD scores in the database marker directory}
\usage{
prepareLDscoresDB(GAlist, ancestry = "EUR", nthreads = 1L,
  overwrite = FALSE, progress = TRUE)
}
\arguments{
\item{GAlist}{A database descriptor with a build, marker directory and rsIDs.}
\item{ancestry}{One or more of EUR, AFR, EAS, SAS or AMR.}
\item{nthreads}{Number of threads passed to gbase.}
\item{overwrite}{Replace an existing score table instead of checking reuse.}
\item{progress}{Report preparation, computation and reuse stages.}
}
\value{Updated GAlist with paths in \code{ldscorePanels}. The descriptor is not
saved automatically. Tables and checksum/provenance JSON files are saved in
\code{marker/}; unchanged inputs reuse the prepared tables.}
\description{Uses the existing ancestry-specific
\code{g1000_<ancestry>.bed/.bim/.fam} files and gbase's direct BED scorer.
No pairwise LD resource is written.}
\details{The settings match the existing gact sparse-LD walkthrough: a distance
of 1,000 selected markers, no physical-distance bound, no r-squared cutoff,
and ordinary uncorrected scores including the unit diagonal. Chromosomes
are scored separately. No additional MAF filter is applied.
Coordinates and allele pairs must match the database catalog. Allele swaps
are allowed because squared correlations are unchanged. Source BED/BIM/FAM
files are never modified. The input checksums, marker selection, build,
ancestry, sample count, gbase version and settings are recorded.
This prepares \code{version = "1000G"} tables; HapMap3 tables remain separate
resources. An up-to-date gbase with \code{task = "ldscores"} is required.
Only required catalog columns are read. Alignment comparisons use bounded
chunks, and temporary catalog/BIM tables are released before genotype
preparation. Output metadata is read again after releasing Glist and unused
model arrays. Native gsuite allocations and numerical settings are unchanged.}
\seealso{\code{\link{getLDscoresDB}}}
