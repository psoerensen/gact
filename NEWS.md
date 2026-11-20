# Development version

* `GAlist <- ingestGWAS(GAlist, selection)` returns the updated descriptor
  after import. Read the per-study report from `GAlist$ingest_report`.
  Dry runs with `import = FALSE` continue to return the readiness report.

# gact 0.0.0.4

* Complete study metadata, including checksums and optional recipe paths, before
  writing `GWAS_information.csv`. Successive studies and metadata-only entries
  retain equal column lengths; omitted sample sizes remain missing.
* Remove the unused inline Shiny application and its eager `readxl` attachment
  from package loading. The packaged ingestion application remains available.
* Require the repaired gbase source release. Installation still uses the public
  gsuite R repository followed by `remotes::install_github("psoerensen/gact")`
  and `gact::install_gsuite()`.
