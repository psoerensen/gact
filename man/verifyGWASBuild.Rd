\name{verifyGWASBuild}
\alias{verifyGWASBuild}
\title{Verify the genome build of a GWAS preview}
\usage{
verifyGWASBuild(preview, map_dir, min_informative = 10L)
}
\arguments{
\item{preview}{Data frame with chr, pos, ref and alt columns.}

\item{map_dir}{Output directory from \code{\link{buildLiftoverMap}}.}

\item{min_informative}{Minimum number of distinct build-informative keys.}
}
\value{
List with the verified build and counts for GRCh37-only, GRCh38-only,
  both, and neither.
}
\description{
Compare distinct coordinate and allele-pair keys with both sides of the map.
A key supports a build when it matches that build only. Verification
requires at least \code{min_informative} such keys and 90 percent agreement.
An unresolved preview has an \code{NA} build and requires review.
}
