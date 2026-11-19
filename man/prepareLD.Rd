\name{prepareLD}
\alias{prepareLD}
\title{Prepare signed sparse LD and ordinary LD scores}
\usage{
prepareLD(GAlist, ancestry = "EUR", nthreads = 1L, window_cm = 1,
  r2 = 0, rsids = NULL, overwrite = FALSE, progress = TRUE, verify = FALSE)
}
\arguments{
\item{GAlist}{Database descriptor pointing to aligned BED panels. BIM A1 must be the database effect allele.}
\item{ancestry}{One or more of EUR, AFR, EAS, SAS or AMR.}
\item{nthreads}{Number of computation threads.}
\item{window_cm}{Inclusive genetic-distance radius within each chromosome, in cM. The BIM must provide an ordered genetic map.}
\item{r2}{Inclusive minimum squared correlation, between zero and one.}
\item{rsids}{Optional unique marker IDs. Selects both endpoints of each pair in BED order. Unknown IDs raise an error.}
\item{overwrite}{Recompute an existing result with these settings.}
\item{progress}{Report stages and chromosome progress.}
\item{verify}{Rehash all files on reuse. By default unchanged size and modification time reuse the recorded checksum.}
}
\description{Calls \code{gbase::gprep()} once per ancestry to prepare sparse signed
correlations and scores from the aligned panel in BED order. The default retains every
pair within 1 cM. Each score is one plus the sum of retained off-diagonal squared correlations.}
\value{Updated GAlist with the selected sparse panel in \code{ldPanels} and matching
score table in \code{ldscorePanels}. Save it to persist the selection.}
\details{The mapped BIM/FAM index is shared across selections. Outputs for distinct settings
and selections have separate paths in the marker directory. Requires gbase with sparse cM support.}
\seealso{\code{\link{getLDscores}}, \code{\link{writeMarkerStat}}}
