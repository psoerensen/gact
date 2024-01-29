
<!-- README.md is generated from README.Rmd. Please edit that file -->

### An R Package for Creating a Database of Genomic Association of Complex Traits

The ***gact*** package is designed for establishing and populating a
comprehensive database focused on genomic associations with complex
traits. The package serves two primary functions: infrastructure
creation and data acquisition. It facilitates the assembly of a
structured repository that includes single marker associations, all
rigorously curated to ensure the high quality of data. Beyond individual
genetic markers, the package integrates a broad spectrum of genomic
entities, encompassing genes, proteins, and an array of biological
complexes (chemical and protein), as well as various biological
pathways.

**gact** provides an infrastructure for efficient processing of
large-scale genomic association data, including core functions for:

- Establishing and populating a database for genomic association.
- Downloading and processing a range of biological databases.
- Downloading and processing summary statistics from genome-wide
  association studies (GWAS).
- Conducting bioinformatic procedures to link genetic markers with
  genes, proteins, metabolites, and biological pathways.
- Finemapping of genomic regions using Bayesian Linear Regression
  models.
- Performing advanced gene set enrichment analysis utilizing a variety
  of tools and methodologies.

**gact** constructs gene and genetic marker sets from a range of
biological databases including:

- `"GO"`: Gene Ontology sets from the [GO](https://geneontology.org)
  database.
- `"Pathways"`: Pathway sets from the [Reactome](https://reactome.org)
  and [KEGG](https://www.genome.jp/kegg/pathway.html) databases.
- `"ProteinComplexes"`: Protein complex sets from the
  [STRING](https://string-db.org) database.
- `"ChemicalComplexes"`: Chemical complex sets from the
  [STITCH](https://stitch.embl.de) database.
- `"DrugGenes"`: Drug-gene interaction sets the
  [DrugBank](https://go.drugbank.com) database.
- `"DrugATCGenes"`: Drug ATC gene sets based on the
  [ATC](https://www.whocc.no/atc_ddd_index/) and
  [DrugBank](https://go.drugbank.com) databases.
- `"DrugComplexes"`: Drug complex sets combining information from
  [STRING](https://string-db.org) and
  [DrugBank](https://go.drugbank.com).
- `"DiseaseGenes"`: Disease-gene sets based on experiments, textmining
  and knowledge base from the
  [DISEASE](https://diseases.jensenlab.org/Search) database.
- `"GTEx"`: GTEx project eQTL sets from the
  [GTEx](https://www.gtexportal.org/home/downloads/adult-gtex/overview)
  database.
- `"GWAScatalog"`: GWAS catalog sets from the
  [GWAScatalog](https://www.ebi.ac.uk/gwas/) database.
- `"Ensembl Regulation"`: Regulatory genomic feature sets from the
  Ensembl Regulation
  [Regulation](https://www.ensembl.org/info/genome/funcgen/index.html)
  database.

**gact** is designed to aid in the biological interpretation of genomic
associations, shedding light on their complex relationships in the
context of genomic associations of complex traits.

![](Images/Concept.png)

### Installation of the ***gact*** package

To install the most recent version of the gact package from GitHub, use
the following commands in R:

``` r
library(devtools)
devtools::install_github("psoerensen/gact")
```

### Tutorials for downloading and installing the gact database

Below is a set of tutorials used for the gact package:

Download and set up the gact database, which is focused on genomic
associations for complex traits:  
[Download and install gact
database](Document/Download_and_install_gact_database.html)

Downloading and processing genome-wide association summary statistic and
ingest into database:  
[Download and process new gwas summary
statistics](Document/Download_and_process_gwas.html)

Download and process genotype data from the 1000 Genomes Project (1000G)
for different ancestries (European, East Asian, South Asian) used in
different genomic analysis:  
[Download and process of 1000G data](Document/Process_1000G.html)

Computing sparse Linkage Disequilibrium (LD) matrices for 1000 Genomes
Project (1000G) data across different ancestries and exploring the LD
data which is used in a number of genomic analysis (LD score regression,
Vegas gene analysis, Bayesian Linear Regression models):  
[Compute sparse LD matrices for 1000G
data](Document/Compute_sparseLD_1000G.html)

### Tutorials for various types of genomic analysis using the gact database

Gene analysis using the VEGAS (Versatile Gene-based Association Study)
approach using the 1000G LD reference data processed above:  
[Gene analysis using VEGAS](Document/Gene_analysis_vegas.html)

Gene set enrichment analysis (GSEA) based on BLR (Bayesian Linear
Regression) model derived gene-level statistics and MAGMA (Multi-marker
Analysis of GenoMic Annotation).  
[Gene set analysis using
BLR-MAGMA](Document/Gene_set_analysis_blr_magma.html)

Pathway prioritization using a BLR-MAGMA model and gene-level statistics
derived from VEGAS.  
[Pathway prioritization using
BLR-MAGMA](Document/Pathway_prioritization_blr_magma.html)

Finemapping of gene regions using single trait Bayesian Linear
Regression models.  
[Finemapping of gene regions using BLR
models](Document/Finemapping_gene_regions_blr.html)

Finemapping of LD regions using single trait Bayesian Linear Regression
models.  
[Finemapping of LD regions using BLR
models](Document/Finemapping_ld_regions_blr.html)

LD score regression for estimating genomic heritability and
correlations.  
[LD score regression](Document/LD_score_regression.html)

#### References

1.  Rohde PD, Sørensen IF, Sørensen P. 2020. qgg: an R package for
    large-scale quantitative genetic analyses. *Bioinformatics* 36:8.
    doi.org/10.1093/bioinformatics/btz955

2.  Rohde PD, Sørensen IF, Sørensen P. 2023. Expanded utility of the R
    package, qgg, with applications within genomic medicine.
    *Bioinformatics* 39:11. doi.org/10.1093/bioinformatics/btad656

3.  Shrestha et al. 2023. Evaluation of Bayesian Linear Regression
    Models as a Fine Mapping Tool. *Submitted*
    doi.org/10.1101/2023.09.01.555889

4.  Bai et al. 2024. Evaluation of multiple marker mapping methods using
    single trait Bayesian Linear Regression models. *In preparation*

5.  Gholipourshahraki et al. 2024. Evaluation of Bayesian Linear
    Regression Models for Pathway Prioritization. *In preparation*
