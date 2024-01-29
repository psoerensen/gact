
<!-- README.md is generated from README.Rmd. Please edit that file -->

### An R Package for Creating a Database of Genomic Association of Complex Traits

The ***gact*** package is designed for establishing and populating a
comprehensive database focused on genomic associations with complex
traits. This function serves two primary roles: infrastructure creation
and data acquisition. It facilitates the assembly of a structured
repository encompassing single marker associations, rigorously curated
to ensure high-quality data. Beyond individual genetic markers, the
function integrates a broad spectrum of genomic entities, encompassing
genes, proteins, and an array of complexes (chemical, protein) as well
as various biological pathways. This integration aims to provide a
holistic view of genomic associations and their multifaceted
relationships in the context of complex traits.

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
Project (1000G) data across different ancestries and exploring this
data:  
[Compute sparse LD matrices for 1000G
data](Document/Compute_sparseLD_1000G.html)

Downloading and processing summary statistics from genome-wide
association studies (GWAS):  
[Incorporate genetic association data into
database](Document/Download_and_process_gwas.html)

Gene analysis using the VEGAS (Versatile Gene-based Association Study)
approach using the 1000G LD reference data processed above:  
[Gene analysis using VEGAS](Document/Gene_analysis_vegas.html)

#### References

1.  Edwards SM, Thomsen B, Madsen P, Sørensen P. 2015. Partitioning of
    genomic variance reveals biological pathways associated with udder
    health and milk production traits in dairy cattle. *Genet Sel Evol*
    47:60. <doi:10.1186/s12711-015-0132-6>  
2.  Edwards SM, Sørensen IF, Sarup P, Mackay TFC, Sørensen P. 2016.
    Genomic prediction for quantitative traits is improved by mapping
    variants to gene ontology categories in *Drosophila melanogaster*.
    *Genetics* 203:1871–1883. <doi:10.1534/genetics.116.187161>
3.  Ehsani A, Janss L, Pomp D, Sørensen P. 2015. Decomposing genomic
    variance using information from GWA, GWE and eQTL analysis. *Anim
    Genet* 47:165–173. <doi:10.1111/age.12396>
4.  Fang L, Sahana G, Ma P, Su G, Yu Y, Zhang S, Lund MS,
    Sørensen P. 2017. Exploring the genetic architecture and improving
    genomic prediction accuracy for mastitis and milk production traits
    in dairy cattle by mapping variants to hepatic transcriptomic
    regions responsive to intra-mammary infection. *Genet Sel Evol*
    49:1–18. <doi:10.1186/s12711-017-0319-0>
5.  Fang L, Sahana G, Su G, Yu Y, Zhang S, Lund MS, Sørensen P. 2017.
    Integrating sequence-based GWAS and RNA-seq provides novel insights
    into the genetic basis of mastitis and milk production in dairy
    cattle. *Sci Rep* 7:45560. <doi:10.1038/srep45560>
6.  Fang L, Sørensen P, Sahana G, Panitz F, Su G, Zhang S, Yu Y, Li B,
    Ma L, Liu G, Lund MS, Thomsen B. 2018. MicroRNA-guided
    prioritization of genome-wide association signals reveals the
    importance of microRNA-target gene networks for complex traits in
    cattle. *Sci Rep* 8:1–14. <doi:10.1038/s41598-018-27729-y>
7.  Ørsted M, Rohde PD, Hoffmann AA, Sørensen P, Kristensen TN. 2017.
    Environmental variation partitioned into separate heritable
    components. *Evolution* (N Y) 72:136–152. <doi:10.1111/evo.13391>
8.  Ørsted M, Hoffmann AA, Rohde PD, Sørensen P, Kristensen TN. 2018.
    Strong impact of thermal environment on the quantitative genetic
    basis of a key stress tolerance trait. *Heredity* (Edinb).
    <doi:10.1038/s41437-018-0117-7>
9.  Rohde PD, Krag K, Loeschcke V, Overgaard J, Sørensen P, Kristensen
    TN. 2016. A quantitative genomic approach for analysis of fitness
    and stress related traits in a *Drosophila melanogaster model
    population*. *Int J Genomics* 2016:1–11.
10. Rohde PD, Demontis D, Cuyabano BCD, The GEMS Group, Børglum AD,
    Sørensen P. 2016. Covariance Association Test (CVAT) identify
    genetic markers associated with schizophrenia in functionally
    associated biological processes. *Genetics* 203:1901–1913.
    <doi:10.1534/genetics.116.189498>
11. Rohde PD, Gaertner B, Ward K, Sørensen P, Mackay TFC. 2017. Genomic
    analysis of genotype-by-social environment interaction for
    *Drosophila melanogaster*. *Genetics* 206:1969–1984.
    <doi:10.1534/genetics.117.200642/-/DC1.1>
12. Rohde PD, Østergaard S, Kristensen TN, Sørensen P, Loeschcke V,
    Mackay TFC, Sarup P. 2018. Functional validation of candidate genes
    detected by genomic feature models. *G3 Genes, Genomes, Genet*
    8:1659–1668. <doi:10.1534/g3.118.200082>
13. Sarup P, Jensen J, Ostersen T, Henryon M, Sørensen P. 2016.
    Increased prediction accuracy using a genomic feature model
    including prior information on quantitative trait locus regions in
    purebred Danish Duroc pigs. *BMC Genet* 17:11.
    <doi:10.1186/s12863-015-0322-9>
14. Sørensen P, de los Campos G, Morgante F, Mackay TFC,
    Sorensen D. 2015. Genetic control of environmental variation of two
    quantitative traits of *Drosophila melanogaster* revealed by
    whole-genome sequencing. *Genetics* 201:487–497.
    <doi:10.1534/genetics.115.180273>
15. Sørensen IF, Edwards SM, Rohde PD, Sørensen P. 2017. Multiple trait
    covariance association test identifies gene ontology categories
    associated with chill coma recovery time in *Drosophila
    melanogaster*. *Sci Rep* 7:2413. <doi:10.1038/s41598-017-02281-3>
