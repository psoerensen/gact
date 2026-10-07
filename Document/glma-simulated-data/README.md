# Simulated-human BIM metadata

`human.bim` accompanies the original simulated-human BED and FAM files from
[qgdata](https://github.com/psoerensen/qgdata/tree/main/simulated_human_data).
It retains all 50,000 rows in their original order, with unchanged chromosome,
marker ID, physical position and alleles. Thirteen negative genetic-map values
are replaced with zero for the gbase metadata contract. These association
examples use physical positions; the genetic-map column is not used by glma.

The phenotype and causal-marker files come from the recorded gact example.
Direct gbase filtering at MAF 0.05, missingness 0.05 and HWE 1e-12, excluding
ambiguous SNPs and indels, reproduces its 37,991 selected markers exactly.
