
<!-- README.md is generated from README.Rmd. Please edit that file -->

### An R Package for Creating a Database of Genomic Association of Complex Traits

The ***gact*** package is designed for establishing and populating a
comprehensive database focused on genomic associations with complex
traits. This function serves two primary roles: infrastructure creation
and data acquisition. It facilitates the assembly of a structured
repository encompassing single marker associations, rigorously curated
to ensure high-quality data. Beyond individual genetic markers, the
function integrates a broad spectrum of genomic entities, encompassing
genes, proteins, and an array of biological complexes (chemical,
protein) as well as various biological pathways. This integration aims
to provide a holistic view of genomic associations and their
multifaceted relationships in the context of complex traits.

**gact** provides an infrastructure for efficient processing of
large-scale genomic association data, including core functions for:

- Establishing and populating a database for genomic association.
- Downloading and processing a range of biological databases.
- Downloading and processing summary statistics from genome-wide
  association studies (GWAS).
- Conducting bioinformatic procedures to link genetic markers with
  genes, proteins, metabolites, and biological pathways.
- Performing gene set enrichment analysis.

### Installation

To install the most recent version of the gact package from GitHub, use
the following commands in R:

``` r
library(devtools)
devtools::install_github("psoerensen/gact")
```

### Tutorials

Below is a set of tutorials used for the gact package:

Download and set up the gact database, which is focused on genomic
associations for complex traits:  
[Download and install gact
database](Document/Download_and_install_gact_database.html)

Processing and preparing genotype data from the 1000 Genomes Project
(1000G) for different ancestries (European, East Asian, South Asian)
used for different genomic analysis:  
[Processing of 1000G data](Document/Process_1000G.html)

Computing sparse Linkage Disequilibrium (LD) matrices for 1000 Genomes
Project (1000G) data across different ancestries and exploring the LD
data which is used in a number of genomic analysis (LD score regression,
Vegas gene analysis, Bayesian Linear Regression models):  
[Compute sparse LD matrices for 1000G
data](Document/Compute_sparseLD_1000G.html)

Downloading and processing summary statistics from genome-wide
association studies (GWAS) and ingest into database:  
[Download and process new gwas summary
statistics](Document/Download_and_process_gwas.html)

Gene analysis using the VEGAS (Versatile Gene-based Association Study)
approach using the 1000G LD reference data processed above:  
[Gene analysis using VEGAS](Document/Gene_analysis_vegas.html)

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
