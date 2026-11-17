test_that("gact builds both releases with their own backbone and coordinates", {
  skip_if(Sys.which("zip") == "", "zip is required to create the annotation fixture")
  root <- tempfile("gact-build-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  shared <- file.path(root, "shared")
  tables <- file.path(root, "reference", "grch38-prep.fixture", "tables")
  release_dir <- file.path(root, "release")
  dir.create(shared)
  dir.create(tables, recursive = TRUE)

  write_gz <- function(lines, file) {
    con <- gzfile(file, "wt")
    on.exit(close(con))
    writeLines(lines, con)
  }
  marker <- data.frame(rsids = c("rs37", "rs38"), chr = 1L, pos = c(100000L, 250000L),
                       ea = "C", nea = "A", eaf = 0.25)
  for (chr in 1:22) {
    path <- file.path(tables, paste0("chr", chr, ".tsv.gz"))
    data.table::fwrite(data.frame(rsids = if (chr == 1L) "rs38" else paste0("rs", chr),
      chr = chr, pos = 250000L, ref = "A", alt = "C", af = 0.25), path)
    jsonlite::write_json(list(build = "GRCh38", chromosome = chr,
      counts = list(retained = 1L), output = list(file = basename(path),
      md5 = unname(tools::md5sum(path)))), file.path(tables, paste0("chr", chr, ".json")),
      auto_unbox = TRUE)
  }

  saveRDS(list(GO = "ENSG1"), file.path(shared, "go.rds"))
  saveRDS(list(Drug = "ENSG1"), file.path(shared, "drug2ensg.rds"))
  saveRDS(list(stale = "rs37"), file.path(shared, "ensg2rsids.rds"))
  archive <- file.path(root, "hsa.0.0.1.zip")
  local({
    previous <- setwd(shared)
    on.exit(setwd(previous))
    expect_equal(utils::zip(archive, c("go.rds", "drug2ensg.rds", "ensg2rsids.rds"),
                            flags = "-q"), 0)
  })
  write_gz(c("gene_stable_id\ttranscript_stable_id\tprotein_stable_id\txref",
             paste("ENSG1", "ENST1", paste0("ENSP", 1:5), "1", sep = "\t")),
           file.path(shared, "Homo_sapiens.GRCh38.109.entrez.tsv.gz"))
  write_gz(c("protein1 protein2 combined_score", paste("9606.ENSP1", paste0("9606.ENSP", 1:5), 999L)),
           file.path(shared, "9606.protein.links.v12.0.txt.gz"))
  write_gz(c("chemical\tprotein\tcombined_score",
             paste("chemical1", paste0("9606.ENSP", 1:5), 999L, sep = "\t")),
           file.path(shared, "9606.protein_chemical.links.v5.0.tsv.gz"))
  writeLines("ENSG1\tR-HSA-1", file.path(shared, "Ensembl2Reactome.txt"))
  writeLines("R-HSA-1\tPathway", file.path(shared, "ReactomePathways.txt"))

  versions <- gact:::gact_versions
  for (version in names(versions)) {
    release <- versions[[version]]
    build <- release$build
    pos <- if (build == "GRCh37") 100000L else 250000L
    sources <- release$sources
    for (resource in c("ensembl", "reactome", "string", "stitch")) {
      paths <- file.path(shared, basename(strsplit(sources[[resource]], ";", fixed = TRUE)[[1L]]))
      if (resource == "ensembl") {
        matching <- paths[grepl(paste0(".", build, "."), basename(paths), fixed = TRUE)]
        gtf <- matching[grepl("gtf.gz$", matching)]
        gff <- matching[grepl("gff.gz$", matching)]
        expect_length(gtf, 1L)
        expect_length(gff, 1L)
        expect_match(basename(gtf), build)
        expect_match(basename(gff), build)
        write_gz(c(paste("#!genome-build", build),
          paste(1, "ensembl_havana", "gene", pos, pos + 1L, ".", "+", ".",
                'gene_id "ENSG1";', sep = "\t")), gtf)
        write_gz(paste(1, "reg", "enhancer", pos, pos, ".", ".", ".", "ID=enhancer:ENSR1",
                       sep = "\t"), gff)
        paths <- paths[file.exists(paths)]
      }
      paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
      file_prefix <- if (.Platform$OS.type == "windows") "file:///" else "file://"
      sources[[resource]] <- paste0(file_prefix, paths, collapse = ";")
    }
    versions[[version]]$sources <- sources
  }

  testthat::local_mocked_bindings(
    .gact_version_metadata = function(version) versions[[version]],
    download_zenodo = function(doi, path, parallel) {
      if (identical(doi, versions[["hsa.GRCh37.0.0.1"]]$sources$marker)) {
        data.table::fwrite(marker, file.path(path, "markers.txt.gz"))
      } else {
        expect_identical(doi, versions[["hsa.GRCh37.0.0.1"]]$sources$annotation)
        expect_true(file.copy(archive, path))
      }
    }, .package = "gact"
  )
  for (version in names(versions)) {
    messages <- capture.output(ga <- gact(version = version, task = "build",
                                          dbdir = release_dir), type = "message")
    expect_true(any(grepl(paste0("Database ready at ", ga$dbdir, " (",
                                 length(ga$rsids), " markers)."), messages, fixed = TRUE)))
    expected <- if (ga$assembly == "GRCh37") "rs37" else "rs38"
    expect_identical(ga$assembly, versions[[version]]$build)
    expect_identical(ga, readRDS(file.path(release_dir, version, paste0("GAlist_", version, ".rds"))))
    expect_true(all(dir.exists(ga$dirs)))
    expect_true(all(startsWith(ga$downloads$dest_file, paste0(ga$dbdir, "/"))))
    expect_true(file.exists(file.path(ga$dirs[["gsets"]], "map.rds")))
    expect_length(list.files(ga$dirs[["ensembl"]], pattern = "\\.(gtf|gff)\\.gz$"), 2L)
    expect_identical(getMarkerSets(ga, "Genesplus")$ENSG1, expected)
    expect_identical(getMarkerSets(ga, "Regulatory Regions")$ENSR1, expected)
    expect_identical(getMarkerSets(ga, "GO")$GO, expected)
    expect_identical(getMarkerSets(ga, "Pathways")$"R-HSA-1", expected)
    expect_identical(getMarkerSets(ga, "DrugGenes")$Drug, expected)
    expect_error(gact(version = version, task = "build", dbdir = release_dir), "already exists")
    if (ga$assembly == "GRCh37") expect_true(file.exists(ga$ldRegions$GRCh37$EUR$file))
    if (ga$assembly == "GRCh38") {
      expect_true(any(grepl("Chromosome 1/22: checking and reading", messages, fixed = TRUE)))
      expect_true(any(grepl("Chromosome 22/22: checking and reading", messages, fixed = TRUE)))
      manifest <- jsonlite::read_json(file.path(ga$dirs[["marker"]], "marker_manifest.json"))
      expect_equal(manifest$n_markers, 22)
      expect_length(manifest$inputs, 22L)
    }
  }
  retry <- file.path(root, "retry")
  unlink(file.path(shared, "ReactomePathways.txt"))
  expect_warning(
    expect_error(gact(version = "hsa.GRCh38.0.0.1", task = "build", dbdir = retry),
                 "Could not download reactome"),
    "ReactomePathways.txt"
  )
  expect_length(list.files(retry, all.files = TRUE, no.. = TRUE), 0L)
})
