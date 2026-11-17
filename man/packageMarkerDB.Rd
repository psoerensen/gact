\name{packageMarkerDB}
\alias{packageMarkerDB}
\title{Package a GRCh38 marker database for release}
\usage{
packageMarkerDB(dbdir, out_dir = getwd())
}
\arguments{
\item{dbdir}{Existing GRCh38 database directory containing the marker directory.}
\item{out_dir}{Output directory outside the source database.}
}
\value{Invisibly, a list containing the archive path, checksum file path,
and archive SHA-256.}
\description{
Create a Zenodo-ready gzip tar archive containing the marker catalog,
coordinate-and-allele index, and provenance. Verify the catalog and index
against the build metadata, and include SHA-256 checksums for the archived
files and the archive itself.
}
\seealso{\code{\link{buildMarkerDB}}}
