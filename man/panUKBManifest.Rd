\name{panUKBManifest}
\alias{panUKBManifest}
\title{Prepare a Pan-UK Biobank manifest for harmonised Catalog files}
\usage{
panUKBManifest(
  selection,
  out_dir,
  phenotype_manifest = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz",
  metadata_dir = NULL
)
}
\arguments{
\item{selection}{Data frame with \code{filename}, \code{accession}, and
\code{analysis} (\code{EUR}, \code{hq}, or \code{all}). Each row identifies a
reviewed Catalog match.}
\item{out_dir}{New or empty directory for the manifest and format profiles.}
\item{phenotype_manifest}{Official compressed Pan-UKB manifest URL or local path.}
\item{metadata_dir}{Optional directory with Catalog per-file YAML and accession-named
\code{.md5} lists. Otherwise fetch them from the Catalog.}
}
\value{Invisibly, paths to the manifest and verified format profiles.}
\description{
Check reviewed Pan-UKB phenotype filename and Catalog accession pairs against
provider metadata. Keep cohort counts as study metadata.
}
