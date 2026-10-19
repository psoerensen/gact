# Consolidated GACT timing evidence

Saved tutorial evidence and authorized gsuite ecosystem measurements from 19 October 2026 are identified separately. The tables describe representative workloads, including a bounded simulated scoring grid, rather than timings of the original UK Biobank manuscript analyses.

## Available measurements

Hardware context: Intel Core i7-1365U, approximately 16 GiB RAM, Windows and R 4.4.1. Exact package versions and hardware records belong to each evidence collection; these runs were not all performed on one software version.

| Operation | Workload | Threads | Wall time | Resident memory | Interpretation | Source |
| --- | --- | --- | --- | --- | --- | --- |
| Database download and setup | 2.57 GiB transfers; 3.04 GiB basic database | 1 requested | 544.20 s | 2.68 GiB | Published earlier-version measurement; network dependent | W; H |
| VEGAS complete initial analysis | CAD/T2D; 17,862 genes; 503 EUR reference people | 4 decoding; serial spectra | 899.94 s | 5.88 GiB | Complete analysis; file input loading separate | V |
| LDSC | CAD/T2D; 7,221,923 common markers | No fitter thread argument | 9.28 s | 3.02 GiB | Sampled analysis peak; inputs loaded before clock | C |
| GNOVA | CAD/T2D; 7,221,923 common markers | No fitter thread argument | 6.39 s | 3.09 GiB | Sampled analysis peak; inputs loaded before clock | C |
| SumHer equal-marker baseline | CAD/T2D; 7,221,923 common markers | No fitter thread argument | 12.87 s | 3.53 GiB | Sampled analysis peak; inputs loaded before clock | C |
| Chromosome-set LDSC preparation and fit | CAD/T2D; 7,221,923 common markers | No fitter thread argument | 29.91 s | 3.27 GiB | Sampled analysis peak; inputs loaded before clock | C |
| HESS regional preparation | 1,702 regions; 18,696 components; moments reused for refits | 4 | 1128.67 s | 3.18 GiB | Sampled analysis peak; inputs loaded before clock | C |
| HESS raw pair refit | 1,702 regions; 18,696 components; moments reused for refits | No fitter thread argument | 3.61 s | 0.41 GiB | Sampled analysis peak; inputs loaded before clock | C |
| HESS partial pair refit | 1,702 regions; 18,696 components; moments reused for refits | No fitter thread argument | 3.31 s | 0.43 GiB | Sampled analysis peak; inputs loaded before clock | C |
| Linear regression | Simulated panel; 5,000 people; 37,991 filtered markers | 4 | 1.64 s | 153 MiB | Single recorded tutorial run | W |
| Gaussian LOCO | Simulated panel; 5,000 people; 37,991 filtered markers | 4 | 57.00 s | 386 MiB | Single recorded tutorial run | W |
| Prior-selected LOCO | Simulated panel; 5,000 people; 37,991 filtered markers | 4 | 124.12 s | 392 MiB | Single recorded tutorial run | W |
| Four regional BayesC fits and posterior sets | Simulated panel; 5,000 people; 37,991 filtered markers | 1 | 5.48 s (5.22–6.28) | 0.56 GiB | Median of 3 serial repeats; timed calls differ from two-chain tutorial calls | F |
| Genome-wide BayesC fit from prepared LD | Simulated panel; 5,000 people; 37,991 filtered markers | 1 | 2.36 s (2.03–2.38) | 0.81 GiB | Median of 3 serial repeats; timed calls differ from two-chain tutorial calls | F |
| Genome-wide credible-set grouping | Simulated panel; 5,000 people; 37,991 filtered markers | 1 | 3.95 s (3.83–4.16) | 0.52 GiB | Median of 3 serial repeats; timed calls differ from two-chain tutorial calls | F |
| Total plus 20 marker-set score application | 2,000 held-out people; 50,000 markers; gsim-based 20K panel | 1 | 14.88 s | Not recorded in these tables | Earlier package snapshot; combined scoring cost; single run | S |
| Total plus 20 marker-set score application | 2,000 held-out people; 50,000 markers; gsim-based 20K panel | 4 | 5.31 s | Not recorded in these tables | Earlier package snapshot; combined scoring cost; single run | S |

Memory definitions differ: C is a sampled analysis-only peak with required inputs loaded; F is median whole-process peak working set; V is complete-worker peak including input loading. W retains the website definition. The scoring tables supply no measured peak. Do not treat these as incremental solver memory.

Fine-mapping controls: regional fits have one chain, 1,000 post-burn sweeps and 100 burn-in sweeps; genome fit has one chain and 100 sweeps without burn-in. Regional marker counts are 850, 785, 1,091 and 1,399. Timings exclude simulation, QC and whole-genome sparse-LD preparation. Completion and short runtimes do not establish convergence.

The scoring snapshot records gscore 0.0.1 and gsim 0.16.0. It uses total plus 20 artificial gene-marker sets, not Reactome pathways. Score application is separate from model fitting. A four-thread limit is not proof that every suboperation used four threads.

## VEGAS preparation and reuse

These are component clocks and alternative reuse measurements, not extra rows to add to the complete 899.94-second workflow. Per-component memory peaks were not recorded.

| Operation | Seconds | Scope |
| --- | ---: | --- |
| marker_set_spectrum_preparation | 643.36 | Once for all 17862 sets; reused across two traits; includes validation, setup, decoding and LD products |
| vegas_prepared_object_conversion | 79.52 | Uses stored spectra; no LD read or eigendecomposition |
| vegas_fit_initial | 97.69 | Both traits; includes statistic validation, gene statistics, tails and result assembly |
| vegas_fit_reused | 71.44 | Separate repeat for both traits; prepared-object loading excluded |
| input_loading | 110.45 | Shown separately from complete analysis and fit |
| marker_and_gene_alignment | 8.59 | Shown separately from spectrum preparation and fit |
| summary_statistic_conversion | 12.71 | Shown separately from spectrum preparation and fit |
| prepared_object_loading | 17.26 | Separate reuse measurement; excluded from repeated fit |

The historical qgg combined workflow took 29,876.32 seconds. The newer workflow shares preparation across traits and uses a different tail evaluator. These are dated measurements, not a contemporaneous identical-method speedup comparison.

## Recovered separate scoring measurements

These saved 4 October records supply separate total and pathway score clocks and whole-process peak R working set. They fill the general scoring-memory gap. They are different workloads from the combined total-plus-20-set run above; their memory must not be attached to that run.

| Operation | Workload | Threads requested | Operation wall time | Process peak memory | Source |
| --- | --- | --- | --- | --- | --- |
| Total score | Artificial-LD panel; 2,000 recipients; 50,000 SNPs | 1 | 8.03 s | 192.3 MiB | A |
| Three synthetic pathway scores | Artificial-LD panel; 2,000 recipients; 50,000 SNPs | 1 | 6.30 s | 192.8 MiB | A |
| Total score | Artificial-LD panel; 2,000 recipients; 50,000 SNPs | 4 | 3.23 s | 192.9 MiB | A |
| Three synthetic pathway scores | Artificial-LD panel; 2,000 recipients; 50,000 SNPs | 4 | 2.66 s | 192.7 MiB | A |
| Total score | 1000G-derived gsim panel; 2,000 recipients; 5,360 SNPs | 1 | 1.03 s | 157.2 MiB | R |
| Three synthetic pathway scores | 1000G-derived gsim panel; 2,000 recipients; 5,360 SNPs | 1 | 2.41 s | 157.8 MiB | R |

A uses binomial/block-copy genotypes and gsim phenotypes. R uses gsim reference-conditioned founders derived from 503 EUR 1000G donors over chr22; 5,360 markers remain after filtering. Neither uses UK Biobank recipients. R uses deterministic sparse simulated weights, not fitted posterior effects. The three pathway memberships are synthetic. Preparation, fitting and serialization are excluded from operation clocks; process memory includes loading and saving. Each setting is one successful execution.

The R collection superseded A for the real-LD benchmark objective. A remains explicitly historical synthetic evidence, including the one/four-thread scoring comparison. Neither establishes current-version scaling or biological prediction performance.

R also supplies a joint Bayesian gene-set fit: **0.10 s and 102.9 MiB**, two synthetic traits and three pathways, one thread. Four chains use 200 burn-in and 800 sampling sweeps with fixed covariance and pattern settings. It is a bounded historical standalone model, not the exact current website call or a manuscript Reactome reproduction.

## Current Bayesian MAGMA and single full-data PoPS models

These reruns use **17,955 saved gene scores and 1,324 Reactome sets** from the already installed hsa.0.0.1 database, gsea 0.1.2 and gbayes 0.1.6. **Gene-association statistics were not recomputed**: their upstream cost is excluded, not measured as zero. Cached fixture loading is separate from the calls below; no new database acquisition, genotype processing or simulation was performed.

| Full-data model | Preparation (s) | BLR fit (s) | Whole call (s) | Process peak (GiB) |
| --- | ---: | ---: | ---: | ---: |
| Bayesian MAGMA: single-trait BayesC | 1.33 | 4.75 | 6.07 | 0.16 |
| Bayesian MAGMA: joint GWAS1/GWAS2 BayesC | 37.81 | 30.99 | 68.86 | 0.52 |
| PoPS: one full-data BayesC model | 1.02 | 6.77 | 7.82 | 0.15 |
| PoPS: one full-data BayesR model | 1.81 | 84.88 | 86.72 | 0.15 |

Preparation includes input checks, set mapping and model-specific summaries. Scalar independent-gene set models compute cross-products directly without forming X. The joint-trait model retains its general projected preparation. Native BLR time includes checks, sampling, posterior accumulation and original-unit coefficient back-transformation. Whole-call time also includes posterior summaries, prediction where relevant and remaining interface overhead; serialization and process startup are excluded. Memory is the complete R worker's lifetime peak, including loading and saving, and cannot be assigned to an individual phase.

All rows use two chains, seeds 10/11 and one requested thread. Single-trait MAGMA uses BayesC with pi=.01 and 5,000 sampling plus 1,000 burn-in sweeps; joint MAGMA uses 1,000 plus 200. MAGMA retains its default fixed hyperparameters; PoPS learns slab variance, residual variance and mixture weights. The two single full-data PoPS models use their current PoPS priors and learning flags, pi=.001, 5,000 sampling plus 1,000 burn-in sweeps. For BayesR that scalar pi is divided among three active components; it is not the .95/.02/.02/.01 mixture from the earlier qgg-oriented probe. Pathway residual normalization remains n-1 for an intercept-only fit. These are representative current workflows, not original manuscript model reruns or convergence evidence.

## Separate 22-fold chromosome-held-out validation

| Validation workflow | Preparation (s) | Sum of BLR fits (s) | Prediction (s) | Whole call (s) | Process peak (GiB) |
| --- | ---: | ---: | ---: | ---: | ---: |
| PoPS BayesC: 22 held-out chromosome folds | 33.64 | 156.90 | 1.27 | 192.97 | 0.24 |
| PoPS BayesR: 22 held-out chromosome folds | 35.82 | 1555.23 | 1.30 | 1593.50 | 0.26 |

These are **22 training models**, each with two chains and 5,000 sampling plus 1,000 burn-in sweeps. They are validation costs and are excluded from the single-fit table. Saved gene scores are reused in every fold; preparation uses only training genes. The native fitting column is the sum over folds, with preparation and held-out prediction recorded separately. Each setting has one execution; there is no repeated-run timing uncertainty estimate.

The updated independent set path needs no dense gene-by-feature matrix. The 1,536 MiB preparation limit is retained for consistency; it is a matrix-payload guard, not measured process memory. The joint model still uses projected preparation. Earlier preparation-included runs and the separate one-chain 14.21-second native probe remain historical evidence, not current rows or matched speedup baselines. A direct sufficient-statistic R interface is now available as `gsea::blr()`, backed by native gbayes. The qgg R/native vb/vg positional mismatch is not reproduced; these timings make no external qgg or MAGMA/FUMA comparison.

See the [rerun report and phase records](../results/20261019-direct-cross-products/README.md) for settings, per-fold timing, binary/input hashes, memory definitions and the excluded duplicate-membership harness diagnostic.

## Historical supplementary evidence

These completed results support the evidence archive but require version and workload checks before inclusion as current website timings. They are not a matched external MAGMA executable comparison.

| Operation | Implementation | Threads requested | Operation seconds | Whole-process peak GiB | Source |
| --- | --- | ---: | ---: | ---: | --- |
| qgg_joint | original gact/qgg | 4 | 0.80 | 4.19 | H |
| qgg_mapping | original gact/qgg | 4 | 95.71 | 4.16 | H |
| suite_joint | standalone counterpart | 4 | 3.11 | 2.93 | H |
| suite_mapping | standalone counterpart | 4 | 40.48 | 3.78 | H |
| website_pops_rr | original gact/qgg | 1 | 2.04 | 0.40 | H |
| website_pops_c | original gact/qgg | 1 | 3.21 | 0.49 | H |
| website_pops_r | original gact/qgg | 1 | 7.71 | 0.49 | H |
| website_magma_b_joint | original gact/qgg | 1 | 0.05 | 0.26 | H |
| website_magma_d_joint | original gact/qgg | 1 | 0.04 | 0.26 | H |
| website_magma_v_joint | original gact/qgg | 1 | 0.05 | 0.26 | H |
| website_magma_b_bayesr | original gact/qgg | 1 | 0.50 | 0.27 | H |
| website_magma_d_bayesr | original gact/qgg | 1 | 0.44 | 0.27 | H |
| website_magma_v_bayesr | original gact/qgg | 1 | 0.48 | 0.27 | H |

Historical joint Bayesian gene-set fits used KEGG and short tutorial controls; they do not establish the manuscript Reactome 10,000-iteration workload. Standalone mapping and joint models differ from their qgg counterparts. PoPS produces gene rankings, not individual polygenic scores.

## Remaining gaps

- Separate per-phase VEGAS memory remains unavailable; complete-worker memory is already measured.
- Scoring has a controlled 40-cell grid, but each cell has one execution; repeated-run uncertainty, larger cohorts and broader workload extrapolation remain unmeasured.

Current Reactome gene-set fits, chromosome-held-out PoPS, separate scoring memory and a bounded scoring grid are now measured. External-method comparisons remain separate from this ecosystem benchmark.

## Sources

All paths are relative to `benchmarks/`; imported-file provenance is in [provenance.csv](../provenance.csv).
- W: [Tutorial_timings_and_memory.Rmd](../scripts/website-snapshot/Tutorial_timings_and_memory.Rmd)
- V: [README.md](../results/gsuite/genomic_20k_50k/manuscript-results/manuscript-vegas-complete-20261005/README.md)
- C: [gcorr-current-timings-memory.csv](../results/website-data/gcorr-current-timings-memory.csv)
- F: [README.md](../results/gsuite/genomic_20k_50k/manuscript-results/manuscript-simulated-finemapping-dense-20261006/README.md)
- S: [summary.csv](../results/gsuite/genomic_20k_50k/results/summary.csv)
- S four threads: [summary-4cores.csv](../results/gsuite/genomic_20k_50k/results/summary-4cores.csv)
- S one thread versions: [provenance.csv](../results/gsuite/genomic_20k_50k/results/provenance.csv)
- S four thread versions: [provenance-4cores.csv](../results/gsuite/genomic_20k_50k/results/provenance-4cores.csv)
- H: [completed-stage-measurements-current.csv](../results/gsuite/genomic_20k_50k/manuscript-results/20261005-structured/completed-stage-measurements-current.csv)

- A: [Artificial-LD measurements](../results/gsuite/genomic_20k_50k/manuscript-results/20261004/measurements.csv) and [scope](../results/gsuite/genomic_20k_50k/manuscript-results/20261004/README.md).
- R: [Reference-derived simulation measurements](../results/gsuite/genomic_20k_50k/manuscript-results/20261004-real-1000g/measurements.csv) and [scope](../results/gsuite/genomic_20k_50k/manuscript-results/20261004-real-1000g/README.md).
