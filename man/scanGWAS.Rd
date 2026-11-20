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
  workers = 1L,
  files = c("catalog", "original"),
  coverage_evidence = NULL,
  metadata_corrections = NULL,
  quality_filters = NULL
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
\item{cohort_evidence}{Optional named character vector of reviewed study-design
evidence keyed by Catalog accession, retained in recipes. Use
\code{metadata_corrections} to supply a reviewed classification.}
\item{GAlist}{Existing GRCh38 database descriptor.}
\item{repository}{\code{"pan-ukb"} or \code{"gwas-catalog"}. A Catalog scan
with no selection uses the published harmonised-file index.}
\item{workers}{Number of concurrent recipe scans. Defaults to 1.}
\item{files}{For Pan-UKB, use \code{"catalog"} (default) for harmonised Catalog
matches or \code{"original"} for the provider's GRCh37 files. Original files
require GAlist and its installed liftover map. Each recipe selects one
ancestry and records the preview build check and liftover index checksum.}
\item{coverage_evidence}{Optional named character vector of reviewed evidence
establishing genome-wide SNP coverage, keyed by Catalog accession.
Explicitly restricted metadata takes precedence. Unclear coverage needs review.
Proteins and metabolites are included when their files meet the same coverage rules.}
\item{metadata_corrections}{List keyed by Catalog accession. Each correction
includes \code{evidence} and any reviewed metadata, for example \code{type}.
For original files supply \code{mapping}, \code{build}, \code{strand = "forward"},
\code{coordinate_system = "1-based"}, \code{effect_scale = "beta"},
\code{se_scale = "beta"} and \code{p_scale}. Supply \code{md5} if the provider
does not publish it. The evidence should document the effect allele and scales.}
\item{quality_filters}{List keyed by accession. Each policy can define
\code{low_confidence = list(column, evidence)} and \code{info = list(column, min, evidence)}.
Pan-UKB uses the selected ancestry's flag when present. Missing values in a
selected filter fail it.}
}
\value{
Invisibly, paths to the import manifest, format profiles and scan report, and a
per-study review table. The manifest includes only ready studies.
\code{recipes_tsv} is the path to the refreshed inventory, and \code{recipe_dir}
can be passed to \code{ingestGWAS()}.
Coverage is \code{genome_wide}, \code{restricted}, or \code{unclear}. Restricted files
have status \code{excluded} and unclear files have status \code{needs_review}.
}
\description{
Match Pan-UK Biobank phenotypes or Catalog accessions to verifiable variant
files using the C++ reader on a bounded gzip prefix. Verified recipes are
keyed by accession and input MD5. Use the same \code{recipe_dir} to reuse
them on later scans. Pan-UKB discovery uses phenotype codes and labels and
checks sex, ancestry and cohort counts against the Catalog metadata.
Sample counts are compared numerically so JSON integer and double
representations of the same count reuse the same recipe.
Ambiguous phenotype matches and shared checksums require review.
The scanner checks upstream metadata again and includes
unresolved studies in the review report. The importer verifies the full
file MD5.
On completion, \code{recipes.tsv} in \code{recipe_dir} is refreshed from all
top-level verified JSON recipes, including those outside the current selection.
It contains one metadata row per recipe. JSONs are not modified. Unverified
drafts are omitted; unreadable or malformed JSONs are skipped with a warning.
The TSV is replaced only after the new file has been fully written. It is an
inventory of cached recipes; the scan manifest records current readiness.
Explicit Catalog case_count/control_count fields are retained as study metadata
and checked against sample_size. Missing breakdowns are not inferred.
Original Pan-UKB scans use the provider manifest's download paths, MD5s,
sex and ancestry-specific counts. Recipes select one ancestry's signed beta,
SE and negative log10 p-value columns. Cross-ancestry meta-analysis columns
remain unmapped. Binary case/control frequencies remain unmapped.
Prescription traits with no description use the provider phenocode as the
label and record that choice. Reserved filename characters are URL-encoded
and use a stable safe local cache filename.
Catalog scans accept standard and legacy harmonised columns after preview
validation. A GWAS-SSF label is not required. Legacy \code{hm_*} coordinates and
alleles take precedence over source columns. Conflicting beta/ratio effects
and missing or inferred mappings require review.
Coverage is screened before previews and recipe reuse. Explicit cis-only,
targeted, exome-only, gene-based and hits-only scope is excluded. When scope
wording is absent, a genome-wide assay plus at least one million successfully
harmonised variants in the provider's running log supplies screening evidence.
Neither an assay label nor a leading preview alone suffices. The report records
the count and log URL. The scan records \code{study_design}, its evidence,
and the ancestry code \code{acn}. Meta-analyses may have a single ancestry
or \code{MIXED}. Unspecified study designs are recorded as \code{unknown}.
Explicitly selected accessions discover available Catalog assets. Reviewed
original gzip TSVs require a successful preview build check.
A unique standard sample-size column (n, N, N_sample_size, sample_size,
n_total or total_sample_size) is mapped to per-variant n. Missing n stays
missing; observed values must be finite and positive. Multiple candidates
require review. Study-level sample counts are not assigned to variants.
}
