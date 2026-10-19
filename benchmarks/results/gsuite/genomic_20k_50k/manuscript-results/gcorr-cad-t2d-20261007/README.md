# CAD/T2D sample size and website filtering, 7 October 2026

The original adapter replaced stat$n with total participants and applied
an extra QC mask. This correction restores extracted n and follows the website
reference filter: BED selected to GAlist$rsids, then every marker in that Glist.
It uses 7,221,923 common reference/summary/LD-score markers, exactly matching
the cached website input universe. No additional MAF, palindromic-allele,
indel, duplicate-position, HWE or missingness exclusion is applied. Full BED
provenance remains 9,985,562 markers. Scores are reused without recalculation.
No new reference/LD computation or qgg timing campaign was run.

| check | count |
| --- | --- |
| reference_markers | 9,985,562 |
| common_database_score_rows | 7,221,923 |
| matched_reference_summary_score_rows | 7,221,923 |
| retained_maf_le_0.05 | 2,358,470 |
| retained_palindromic | 1,083,607 |
| retained_duplicated_positions | 0 |

## Verified source semantics

The raw GWAS1/GWAS2 database files have no n column. getMarkerStat supplies
study metadata neff, constant at 49071.273294961 and 40743.152404981.
Binary-study metadata is Ncase*Ncontrol/(Ncase+Ncontrol). Ingestion can also
infer a missing summary n from round(median(1/(2*af*(1-af)*seb^2))).
The adapter now retains extracted marker-specific n without averaging or
substitution. A focused fixture verified unequal n across markers, allele
swapping, the common marker mask and stored-score reuse.

## Corrected ordinary LDSC

| quantity | trait1 | trait2 | estimate | standard_error |
| --- | --- | --- | --- | --- |
| h2 | GWAS1 |  | 0.469755 | 0.0246419 |
| h2 | GWAS2 |  | 0.310603 | 0.0192085 |
| cov | GWAS1 | GWAS2 | 0.146811 | 0.0119046 |
| rg | GWAS1 | GWAS2 | 0.384344 | 0.028498 |

## Saved qgg context

The original qgg results were T2D h2 0.317088743 and CAD h2 0.271063762.
Their chromosome sums were 0.328760832 and 0.275251463. These are existing
5 October outputs, not a new qgg timing campaign. The website-filtered gcorr
input now uses the same 7,221,923-row base marker universe and stat$n. Ordinary
LDSC applies qgg's default squared-Z threshold (80 here), and set LDSC applies
no Z cutoff, as in qgg's set branch. Set normalization uses matched set lengths.
Remaining estimator differences: gcorr is weighted and uses a coordinated
block jackknife; qgg is unweighted. The gcorr ordinary joint fit has one common
pre-filter marker-count parameter; qgg uses separate post-filter marginal/pair
counts. Native Z filtering is inclusive; qgg is strict. Neither agreement of
inputs nor these h2 values by themselves establish a biological trait scale.

## Remaining limits

Old GNOVA, SumHer and HESS numerical results are superseded. The sample-size-only
7 October diagnostic retaining the extra QC/full-reference counts is also
superseded; its successful runs remain preserved locally. The remaining methods
have now been refitted under stat$n and the website marker selection; see the
[completed remaining analyses](remaining-analyses/README.md).
The historical full HESS projection retains 83728 components, exceeding both
information sizes. Its correction requires 1-total_rank/n > 0, so the historical
50-component cap is inadmissible. The corrected genome illustration explicitly
uses a 20-component cap, retains 33971 components, and passes rank admission.
Its CAD total is nevertheless negative, genome rg is undefined, and conditional
uncertainty is unavailable. These HESS results do not support interpretable
CAD/T2D genetic-correlation inference under that specification.
Do not bypass guards or rescale old HESS points as corrected results. Previous
source and timing archives remain intact; their resource records describe those runs.
The native/R source baseline is preserved in the neighboring 6 October
source-and-evidence.zip and chromosome-sets-source-and-evidence.zip archives.

## Completed corrected timings

| analysis | fit_seconds | analysis_seconds | process_wall_seconds | peak_working_set_GiB |
| --- | --- | --- | --- | --- |
| gcorr_inputs | unavailable |  377.66 | 434.658 | 5.09188 |
| gcorr_ldsc |  353.26 |  370.55 |  478.33 |  6.4394 |
| gcorr_sets |     383 |  412.19 | 510.646 | 7.56313 |

Completed run: build/benchmarks/manuscript-20261007-003431
The frozen successful workers and compact measurements accompany this record.
No raw marker-level summaries, genotypes, n matrix or installed binaries are exported.


## Corrected chromosome-set LDSC

| Set | T2D h2 | T2D SE | CAD h2 | CAD SE | rg | rg SE |
| --- | --- | --- | --- | --- | --- | --- |
| chr1 | 0.03251 | 0.005642 | 0.02872 | 0.004129 | 0.6071 | 0.08725 |
| chr2 | 0.04199 | 0.008042 | 0.03304 | 0.004827 | 0.4674 | 0.08525 |
| chr3 | 0.04189 | 0.006333 | 0.02131 | 0.005017 | 0.6509 | 0.06792 |
| chr4 | 0.02636 | 0.005503 | 0.0133 | 0.003774 | 0.6419 | 0.1587 |
| chr5 | 0.02524 | 0.006452 | 0.01129 | 0.002269 | 0.5833 | 0.1693 |
| chr6 | 0.04417 | 0.01143 | 0.02908 | 0.005188 | 0.4172 | 0.06901 |
| chr7 | 0.02083 | 0.003106 | 0.01397 | 0.0036 | 0.7241 | 0.1476 |
| chr8 | 0.02903 | 0.007399 | 0.007274 | 0.003444 | 1.208 | 0.3434 |
| chr9 | 0.02569 | 0.005309 | 0.02091 | 0.007846 | 0.5497 | 0.113 |
| chr10 | 0.04108 | 0.02039 | 0.02071 | 0.005112 | 0.3198 | 0.06986 |
| chr11 | 0.03066 | 0.003407 | 0.01957 | 0.002429 | 0.5008 | 0.07992 |
| chr12 | 0.02821 | 0.004672 | 0.01828 | 0.002508 | 0.6103 | 0.0787 |
| chr13 | 0.01127 | 0.002215 | 0.006964 | 0.003109 | 0.4508 | 0.1693 |
| chr14 | 0.00948 | 0.003218 | 0.007679 | 0.002181 | 1.022 | 0.1813 |
| chr15 | 0.02171 | 0.005343 | 0.01147 | 0.00472 | 0.5049 | 0.1386 |
| chr16 | 0.02172 | 0.006993 | 0.008354 | 0.003284 | 0.6818 | 0.1863 |
| chr17 | 0.01744 | 0.004029 | 0.01462 | 0.0039 | 0.7755 | 0.103 |
| chr18 | 0.01238 | 0.00353 | 0.004471 | 0.001211 | 0.7948 | 0.4394 |
| chr19 | 0.008318 | 0.001134 | 0.01363 | 0.006147 | 0.2608 | 0.1058 |
| chr20 | 0.007753 | 0.003317 | 0.009448 | 0.001506 | 0.5481 | 0.1867 |
| chr21 | 0.0008346 | 0.0007212 | 0.004773 | 0.002876 |  1.09 | 0.734 |
| chr22 | 0.007865 | 0.001273 | 0.005257 | 0.002314 | 0.3415 | 0.1294 |


Joint model totals:

| quantity | trait1 | trait2 | estimate | standard_error |
| --- | --- | --- | --- | --- |
| h2 | GWAS1 |  | 0.5064 | 0.02925 |
| h2 | GWAS2 |  | 0.3241 | 0.02089 |
| cov | GWAS1 | GWAS2 | 0.2303 | 0.009884 |
| rg | GWAS1 | GWAS2 | 0.5685 | 0.02602 |


Set-score preparation took 9.54 seconds; fitting took 383.00 seconds. The analysis clock was 412.19 seconds, process wall 510.65 seconds, and peak R working set 7.56 GiB. These are single observed runs; the configured budget is four threads, but this fitter exposes no thread argument.


0 of 22 set correlations have unavailable point estimates, 0 have unavailable standard errors, and 3 finite raw ratios lie outside [-1,1]. Raw negative h2 and invalid ratios remain visible; intervals in the plot are estimate +/- 1.96 SE, without clipping or a calibrated local-significance claim.


![Chromosome-set heritability and correlation](gcorr-chromosome-set-results.png)


The 22 h2/covariance contributions sum to their joint totals; rg is unchanged by constant trait-specific n scaling.
