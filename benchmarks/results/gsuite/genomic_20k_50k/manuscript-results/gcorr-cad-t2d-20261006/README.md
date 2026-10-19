# CAD/T2D gcorr workflow, 6 October 2026

**Sample-size correction, 7 October 2026:** numerical results below used total
participants in place of extracted `stat$n`. They are superseded and preserved
as historical evidence. Audited information sizes are 49071.273294961 (GWAS1)
and 40743.152404981 (GWAS2). See `../gcorr-cad-t2d-20261007/README.md` for corrected
LDSC results and qgg convention differences. The 83728-component full-genome
HESS projection exceeds both corrected information sizes, making its joint
correction inadmissible at the historical spectral settings.

This evidence records the requested gact tutorial using the standalone gcorr
R package. It reuses the existing CAD/T2D summaries, ordinary EUR LD scores
and 503-person BED reference. No GWAS/reference download, database ingestion,
genome-wide LD matrix or persistent eigenvector cache was created.

The runnable page is maintained in the sibling gact repository at
`Document/Genetic_correlation_CAD_T2D_gcorr.Rmd`, with compact downloads in
`Document/data/`. The protocol and benchmark worker are maintained two levels
above this directory. `source-and-evidence.zip` captures the relevant public
sources and these compact measurements; the accompanying SHA-256 inventory
records their exact contents. Raw inputs, installed libraries and logs are not
included in that archive.

## Inputs and model

- GWAS1: T2D without UK Biobank, PMID30297969, recorded total N = 456,236.
- GWAS2: CAD, PMID26343387, recorded total N = 184,305.
- Reference: 9,985,562 markers; 503 EUR individuals; GRCh37 coordinates.
- Shared QC: 4,129,274 markers, MAF > 0.05, nonpalindromic nucleotide alleles,
  unique positions, database/BIM coordinates checked and Z oriented to BIM
  allele 1. The existing three-study GWAS1/GWAS2/GWAS6 availability mask is
  retained; only GWAS1 and GWAS2 are fitted.
- Regional partition: pinned Berisa/Pickrell EUR boundaries, 1,702 nonempty
  disjoint regions, 24 on chromosome 22, no uncovered retained markers;
  largest region 7,801 markers.
- HESS spectral controls: at most 50 components, eigenvalues > 1, relative
  tolerance 1e-10. The smaller Gram is at most 503 by 503 here. Both traits
  share each eigendecomposition; only identities, spectra and moments survive.

The working Z/sqrt(N) scale uses recorded total participants rather than the
database's n_case*n_control/N field. It is not liability-scale heritability or
a matched numerical/performance comparison with the earlier qgg timings.
Per-variant sample sizes were not independently verified. CAD's multiancestry
meta-analysis is imperfectly matched by the EUR reference. Legacy LD-score
window/correction provenance remains unresolved.

LDSC estimates marginal and cross intercepts. GNOVA conditions on the fitted
LDSC cross intercept and uses its supplied regression-marker universe and a
unit marginal noise baseline. SumHer is an illustrative equal-marker baseline,
not a LDAK architecture/software comparison. These are conditional workflow
illustrations, not newly qualified manuscript estimates.

## Completed timings

| Analysis | Preparation (s) | Fits (s) | Analysis total (s) | R process wall (s) | Peak working set (GiB) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Chromosome 22 HESS/rho-HESS | 193.29 | 0.17 | 195.47 | 234.16 | 3.86 |
| LDSC | included in fit | 234.92 | 248.89 | 315.48 | 5.19 |
| GNOVA | included in fit | 122.47 | 129.00 | 194.45 | 6.09 |
| SumHer baseline | included in fit | 220.70 | 223.41 | 289.33 | 5.83 |
| Full-genome HESS/rho-HESS | 649.43 | 49.81 | 702.19 | 742.34 | 5.86 |

The full regional fit column includes two marginal fits (11.46 s), the
zero-overlap pair fit (8.07 s), and two additional overlap sensitivity fits
(30.28 s). R adapters and result formatting are included. Regional preparation
consists of R validation/matching 67.42 s, BED reader setup 108.09 s, genotype
decoding 75.58 s, Gram products 97.74 s, eigendecomposition 290.38 s,
both-trait projection 2.74 s, and other checks/adapters 7.45 s.

Hardware: Windows 11 Enterprise, Intel Core i7-1365U, 10 physical/12 logical
cores, 16 GB RAM, R 4.4.1, compatible portable development installation with
native/R gcorr 0.1.1. The budget is four threads; regional decoding/products
receive nthreads=4, while the score fitters expose no thread argument.
Measurements concern the actual R process. Analysis clocks exclude saved-input
loading and output serialization; process wall clocks include them. Working
set and sampled private memory are separate metrics, both retained in CSV.
These are single observed runs, not a broad benchmark campaign.

## Interpretation limits

Genome-wide correlation estimates are LDSC 0.3844 (SE 0.0345), GNOVA 0.4696
(SE 0.0393), and the illustrative SumHer baseline 0.3852 (SE 0.0343),
conditional on the distinct model/intercept conventions above.

Shared-sample count and phenotypic correlation are not recorded in the
database. The regional zero-overlap scenario is an explicit assumption; an
LDSC cross intercept is not a phenotypic correlation. The sensitivity scenarios
(25,000 shared, rP = 0.1) and (50,000 shared, rP = 0.3) are also assumed values.

**The zero-overlap full-genome HESS correlation ratio is 1.6578, outside
[-1,1], with an indefinite estimated total genetic covariance matrix. It
cannot be interpreted as a biological genetic correlation.** There are 1,397
negative local heritability estimates; all 6,812 regional/total standard errors
in this pair fit are unavailable under the conditional plug-in checks. The
50,000/rP=0.3 scenario gives a raw total ratio of 0.8718 but still has no valid
conditional uncertainty. Raw points and explicit reasons are retained; no PSD
repair or clipping was applied. Actual overlap, reference compatibility and
spectral-model assumptions require review before biological interpretation.
The model also assumes negligible inter-region LD and treats supplied regions
as the complete genetic universe. The chr22 total is a chromosome-only pilot.

Focused native dense-projection checks passed for marker/sample spaces and
1/4 threads; a small R BED-derived dense oracle passed, including allele swaps,
interleaved region mapping and RDS save/reload. Standalone load/DLL checks
passed. This is numerical equivalence evidence, not statistical calibration or
complete source-package/ecosystem qualification.

## Local provenance and reproduction

Paths below are relative to the gsuite checkout:

- Completed LDSC/GNOVA/SumHer stages:
  `build/benchmarks/manuscript-20261006-213249`.
- Completed chr22 and full-genome regional stages:
  `build/benchmarks/manuscript-20261006-220304`.
- Independently validated input cache:
  `build/benchmarks/manuscript-20261006-211254/gcorr_inputs-4/output.rds`.

Input extraction saved its cache before a later parser failure caused by
editing a running worker. Its markers, region membership, coordinates,
orientation, sample sizes and stored scores were independently validated before
reuse; that failed input process is excluded from completed timings. The
launcher now freezes each worker before starting R. An earlier stopped regional
run repeatedly hashed the full marker universe in an R per-region lookup;
the mapping now runs once. That incomplete run is also excluded. Existing
accepted input caches remain in place. Failed-run logs and provenance are
retained in `build/benchmarks/final-summaries/`.

To regenerate compact tables and plots without repeating any fits:

```sh
Rscript --vanilla benchmarks/genomic_20k_50k/summarize_gcorr_cad_t2d.R \
  build/benchmarks/manuscript-20261006-220304 \
  build/benchmarks/manuscript-20261006-213249 \
  build/benchmarks/manuscript-20261006-220304 \
  build/benchmarks/manuscript-20261006-211254/gcorr_inputs-4 \
  build/gact-revision-cleanup/tutorial-output
```

The protocol documents bounded reruns through the existing process-measured
launcher. Neither the compact archive nor the public tutorial contains the
original GWAS/genotype inputs required for a full rerun.

## Chromosome marker-set analysis

The user stated their understanding that CAD and T2D participants do not
overlap. The working model therefore uses no shared samples and fixes the
cross-trait LDSC intercept to zero; cohort membership has not been independently
audited. A free cross intercept can also absorb shared confounding, so this
constraint is a model choice beyond the sample-overlap assertion. Marginal
intercepts are estimated jointly across all 22 chromosome marker sets. No qgg
comparison was requested for this extension.

Chromosomes enter the fit as an ordinary named `sets` list. The generic helper
`gcorr_set_scores.R` masks the existing ordinary scores into disjoint set
columns, using the within-chromosome LD approximation. Arbitrary feature sets
with LD across their boundaries require true annotation LD scores instead.
Original reference counts were read from the BIM file (one file contains all
22 chromosomes), preserving the earlier reference normalization. The 22
heritability/covariance contributions sum to the joint model totals within
1e-10; correlation ratios are not additive. No new eigendecomposition or LD
files were created, and the set-score matrix was not saved to disk.

| Quantity | Trait 1 | Trait 2 | Estimate | SE |
| --- | --- | --- | ---: | ---: |
| h2 | GWAS1 |  | 0.084820 | 0.006064 |
| h2 | GWAS2 |  | 0.106882 | 0.009444 |
| cov | GWAS1 | GWAS2 | 0.051056 | 0.002700 |
| rg | GWAS1 | GWAS2 | 0.536221 | 0.034867 |

Set-score preparation: 10.01 s; fit: 207.11 s; analysis: 236.82 s; R process wall: 299.20 s; peak working set: 5.64 GiB.

Undefined chromosome correlation points: none. Unavailable correlation SEs: none. Raw chromosome ratios outside [-1,1]: chr8, chr21.

These set-model results differ from the earlier freely estimated cross
intercept fits. They remain conditional workflow illustrations, with the same
binary-trait scaling, legacy-score and reference-ancestry limitations. The
earlier zero-overlap HESS ratio of 1.6578 remains unsuitable for biological
interpretation; missing phenotypic correlation is not a blocker when Ns=0.

The completed set run is `build/benchmarks/manuscript-20261006-224711`.
An earlier metadata-preparation attempt (`manuscript-20261006-224455`)
stopped before fitting because the initial count adapter assumed one chromosome
per reference file. It is excluded from measurements; its logs and provenance
remain in `build/benchmarks/final-summaries/`.
The numerical engines and installed native/R gcorr 0.1.1 were reused unchanged.
The focused set-helper oracle checked reordered disjoint sets, masked-score
identity, reference counts and known dense-fit heritability contributions.

The chromosome tables/plot are the `gcorr-chromosome-set-*` files alongside
this README. `chromosome-sets-source-and-evidence.zip` and its SHA-256 inventory
capture the extended sources, frozen successful worker and compact evidence.
The original `source-and-evidence.zip` and checksum remain preserved.

To regenerate the added results without fitting again:

```sh
Rscript --vanilla benchmarks/genomic_20k_50k/summarize_gcorr_sets.R \
  build/benchmarks/manuscript-20261006-224711 \
  build/gact-revision-cleanup/tutorial-output
```
