\name{scanGWAS}
\alias{scanGWAS}
\title{Scan Pan-UK Biobank and GWAS Catalog studies for ingestion}
\usage{
scanGWAS(
  selection = NULL,
  out_dir = NULL,
  phenotype_manifest = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz",
  preview_rows = 1000L,
  max_preview_bytes = 1024L * 1024L,
  recipe_dir = NULL,
  source = c("panukb", "catalog"),
  analysis = "per_ancestry",
  cohort_evidence = NULL,
  GAlist = NULL,
  repository = NULL,
  workers = 1L
)
}
\arguments{
\item{selection}{Pan-UKB data frame with filename, an ancestry code in the
analysis column, and optional accession, or a character vector of Catalog
accessions. Omit it to select all available phenotype and ancestry pairs.}
\item{out_dir}{New or empty scan report directory. Defaults to a new
directory under \code{recipe_dir/reviews}.}
\item{phenotype_manifest}{Official compressed Pan-UKB manifest URL or local path.}
\item{preview_rows}{Number of complete leading rows to inspect.}
\item{max_preview_bytes}{Maximum compressed bytes requested per file.}
\item{recipe_dir}{Persistent recipe directory. Defaults to the database's
recipes directory when \code{GAlist} is supplied, otherwise \code{out_dir/recipes}.}
\item{source}{\code{"panukb"} or \code{"catalog"}.}
\item{analysis}{\code{"per_ancestry"} (default) or one Pan-UKB ancestry code
for scans of the full phenotype manifest.}
\item{cohort_evidence}{Optional named character vector of reviewed
single-cohort evidence keyed by Catalog accession, when file metadata lacks
explicit evidence. Meta-analyses and mixed ancestry files are excluded.}
\item{GAlist}{Existing GRCh38 database descriptor.}
\item{repository}{\code{"pan-ukb"} or \code{"gwas-catalog"}. A Catalog scan
with no selection uses the published harmonised-file index.}
\item{workers}{Number of concurrent recipe scans. Defaults to 1.}
}
\value{
Invisibly, paths to the import manifest, format profiles and scan report, and a
per-study review table. The manifest includes only ready studies.
}
\description{
Match Pan-UK Biobank phenotypes or Catalog accessions to verified GRCh38
single-ancestry files using the C++ reader on a bounded gzip prefix. Verified recipes are
keyed by accession and input MD5. Use the same \code{recipe_dir} to reuse
them on later scans. Pan-UKB discovery uses phenotype codes and labels and
checks sex, ancestry and cohort counts against the Catalog metadata.
Ambiguous phenotype matches and shared checksums require review.
The scanner checks upstream metadata again and includes
unresolved studies in the review report. The importer verifies the full
file MD5.
}
