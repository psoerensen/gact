\name{ingestGWAS}
\alias{ingestGWAS}
\title{Ingest explicitly selected verified GWAS recipes}
\usage{
ingestGWAS(GAlist, selection, workers = 1L, import = TRUE,
  descriptor = NULL, data_dir = NULL, out_dir = NULL, scratch_dir = NULL)
}
\arguments{
\item{GAlist}{Database descriptor list, or path to its saved RDS.}
\item{selection}{Required accessions or recipe filenames, or selected rows from
recipes.tsv. Character selections must each match exactly one recipe. Use
recipe_file, ancestry and/or input_md5 to disambiguate data-frame rows.
TSV metadata edits do not override the JSON recipe.}
\item{workers}{Number of C++ staging workers. Windows uses one worker.}
\item{import}{Whether to import. FALSE checks selected recipes without
downloading full inputs or changing the database.}
\item{descriptor}{Saved GAlist RDS path. By default find the single GAlist*.rds
in the database root. The current saved descriptor is reloaded.}
\item{data_dir}{Optional input cache directory; defaults to download/.}
\item{out_dir}{Optional new report directory; defaults to recipes/reviews/.}
\item{scratch_dir}{Optional parent of temporary native working files.}
}
\value{
Invisibly, the per-study report from \code{\link{runGWASManifest}}, including
imported, already_imported, ready or needs_review status and local study IDs.
}
\description{
Reads recipes.tsv and uses selected verified JSON recipes as the source of study
metadata and mappings. Selection is mandatory. Existing imports are checked in
GWAS_information.csv by accession, ancestry and input checksum. Matching imports
are skipped; changed inputs or recipes require review. Local GWAS IDs are assigned
on publication. The CSV and saved descriptor are updated after successful imports.
JSON recipes and recipes.tsv remain unchanged. Study totals and case/control
counts populate study metadata; missing per-variant n remains missing.
}
