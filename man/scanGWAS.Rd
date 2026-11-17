\name{scanGWAS}
\alias{scanGWAS}
\title{Scan Pan-UK Biobank studies for GWAS ingestion}
\usage{
scanGWAS(
  selection,
  out_dir,
  phenotype_manifest = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz",
  preview_rows = 1000L,
  max_preview_bytes = 1024L * 1024L
)
}
\arguments{
\item{selection}{Data frame with Pan-UKB filename and analysis (\code{EUR},
\code{hq}, or \code{all}). An optional accession column supplies an explicit
Catalog match.}
\item{out_dir}{New or empty directory for recipes and scan reports.}
\item{phenotype_manifest}{Official compressed Pan-UKB manifest URL or local path.}
\item{preview_rows}{Number of complete leading rows to inspect.}
\item{max_preview_bytes}{Maximum compressed bytes requested per file.}
}
\value{
Invisibly, paths to the import manifest, format profiles and scan report, and a
per-study review table. The manifest includes only ready studies.
}
\description{
Match selected Pan-UK Biobank phenotypes to verified GRCh38 GWAS Catalog
files. Use the C++ reader on a bounded gzip prefix and prepare recipes and an
import manifest. The scanner does not open a database. The importer verifies
the full file MD5.
}
