# AGENTS.md

## Ecosystem orientation

For cross-repository work, read the optional local
[agent overview](../gfactory/docs/agent-overview.md), then only the relevant
owner contract and current source. It gives ownership, reuse links and current
priorities; this repository retains its own interfaces and build instructions.
If private gfactory is absent, continue with local guidance; do not fetch it or
make builds, tests, runtime, public evidence or website rendering depend on it.
A cross-repository link does not authorize edits in another repository.

## Local responsibility

gact owns R workflows for acquiring, organizing and using genomic association
data and biological annotations. Keep shared numerical/statistical engines in
their owning native libraries. Existing gact analyses remain supported locally;
this ownership guidance does not authorize moving them.

## Local development

Read [README.Rmd](README.Rmd) and [DESCRIPTION](DESCRIPTION). README.md is
rendered from README.Rmd; change the source when editing that description.
For R code, parse changed files first. After an authorized package change,
install into an explicit isolated library with
`R CMD INSTALL --library="<isolated-library>" .`, then run the relevant test
in tests/testthat with that library first in .libPaths() and gact loaded.
Documentation-only edits do not require package compilation or data acquisition.

Check Git status before editing and preserve unrelated work. Do not download
GWAS/reference data, run ingestion campaigns, modify databases, render the
whole website or change dependencies merely to check documentation. Keep data
and credentials out of reports. Commit/push and publication require a request.
