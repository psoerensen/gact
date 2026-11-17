test_that("the same set builders and readers follow each release's coordinates", {
  scratch <- tempfile("gact-sets-")
  dir.create(scratch)
  write_gz <- function(path, lines) {
    con <- gzfile(path, "wt")
    on.exit(close(con))
    writeLines(lines, con)
  }
  for (build in c("GRCh37", "GRCh38")) {
    ga <- createDB(paste0("hsa.", build, ".0.0.1"), dbdir = scratch)
    ensembl <- ga$dirs[["ensembl"]]
    sets <- ga$dirs[["gsets"]]
    for (assembly in c("GRCh37", "GRCh38")) {
      pos <- if (assembly == "GRCh37") 100L else 200L
      gtf <- if (assembly == "GRCh37") "GRCh37.87.gtf.gz" else "Homo_sapiens.GRCh38.109.gtf.gz"
      gff <- if (assembly == "GRCh37") "GRCh37.Regulatory_Build.regulatory_features.gff.gz" else
        "homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz"
      write_gz(file.path(ensembl, gtf), c(paste("#!genome-version", assembly), "# comment",
        paste(1, "ensembl_havana", "gene", pos, pos + 10L, ".", "+", ".",
              'gene_id "ENSG1"; gene_name "GENE";', sep = "\t")))
      write_gz(file.path(ensembl, gff), c("##gff-version 3",
        paste(1, "reg", "enhancer", pos, pos, ".", ".", ".",
              "ID=enhancer:ENSR1;description=test", sep = "\t")))
    }
    markers <- data.frame(chr = 1L, pos = c(100L, 110L, 200L, 210L),
                          rsids = paste0("rs", 1:4))
    data.table::fwrite(markers, file.path(ga$dirs[["marker"]], "markers.txt.gz"))
    for (input in c("go", "reactome2ensg", "string2ensg", "stitch2ensg", "drug2ensg")) {
      saveRDS(list(set1 = "ENSG1"), file.path(sets, paste0(input, ".rds")))
    }
    ga$map <- list(ensg2ensp = list(ENSG1 = "ENSP1"))
    expected <- if (build == "GRCh37") c("rs1", "rs2") else c("rs3", "rs4")
    ga <- createSetsDB(ga, what = c("gene", "regulatory", "marker"),
                       upstream = 0, downstream = 0, overwrite = TRUE)
    expect_identical(ga$sets_metadata$build, build)
    expect_identical(getMarkerSets(ga, "Genesplus")$ENSG1, expected)
    for (feature in c("Pathways", "GO", "ProteinComplexes", "ChemicalComplexes", "DrugGenes",
                      "DrugComplexes")) {
      expect_identical(getMarkerSets(ga, feature)$set1, expected)
    }
    expect_identical(getMarkerSets(ga, "Regulatory Regions")$ENSR1, expected[1L])
    expect_identical(getMarkerSets(ga, "Regulatory Categories")$enhancer, expected[1L])
    expect_identical(readRDS(file.path(sets, "ensr2ensg.rds"))$ENSR1, "ENSG1")
    expect_identical(readRDS(file.path(sets, "ensg2ensr.rds"))$ENSG1, "ENSR1")
    other <- if (build == "GRCh37") "GRCh38" else "GRCh37"
    expect_error(createRegulatorySets(ga, build = other), "consistent")
    ga$assembly <- other
    expect_error(createSetsDB(ga, what = "gene"), "consistent")
  }
})
