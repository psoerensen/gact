# Bayesian MAGMA and PoPS: direct-cross-product reruns

Measured on 19 October 2026 from the same saved hsa.0.0.1 gene statistics and Reactome annotation used by the earlier benchmark collection. No gene statistics, database contents or genotypes were recomputed or downloaded. The cached fixture is hashed in provenance.json.

## Single full-data models

| Workflow | Preparation s | Native BLR s | Posterior summaries s | Whole call s | R process peak GiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| Bayesian MAGMA: single-trait BayesC | 1.33 | 4.75 | 0.01 | 6.07 | 0.16 |
| Bayesian MAGMA: joint GWAS1/GWAS2 BayesC | 37.81 | 30.99 | 0.01 | 68.86 | 0.52 |
| PoPS: one full-data BayesC model | 1.02 | 6.77 | 0.02 | 7.82 | 0.15 |
| PoPS: one full-data BayesR model | 1.81 | 84.88 | 0.02 | 86.72 | 0.15 |

## Separate chromosome-held-out validation

| Workflow | Preparation s | Sum BLR s | Posterior summaries s | Prediction s | Whole call s | R process peak GiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| PoPS BayesC: 22 held-out chromosome folds | 33.64 | 156.90 | 0.34 | 1.27 | 192.97 | 0.24 |
| PoPS BayesR: 22 held-out chromosome folds | 35.82 | 1555.23 | 0.31 | 1.30 | 1593.50 | 0.26 |

## Inputs, settings and boundaries

All fits reuse 17,955 saved genes and 1,324 Reactome sets. Per-fold training universes and varying features are recorded separately; 22 folds do not imply 22 recomputations of gene evidence. No phenotype or gene-score simulation is performed. The separately timed 17,862-gene VEGAS workflow has different inputs and must not be added to these clocks.

Every model uses two chains, seeds 10/11 and one requested thread. Single-trait Bayesian MAGMA is BayesC with pi=.01, 5,000 sampling sweeps and 1,000 burn-in. Joint GWAS1/GWAS2 BayesC uses 1,000 sampling sweeps and 200 burn-in. Both retain the interface's fixed hyperparameters. Joint traits retain general projected preparation and independent default trait errors.

The full-data PoPS models and each validation fold use pi=.001, 5,000 sampling sweeps, 1,000 burn-in, thinning one and learned slab/residual variances and mixture weights. BayesC uses active probability .001; BayesR starts with weights (.999,.001/3,.001/3,.001/3) and multipliers (0,.01,.1,1). Slab variance starts at .5/(feature_count*weighted_active_mass), residual variance .5, variance priors df=4 with scale equal to the corresponding starting variance, and unit Dirichlet concentrations. Variable feature counts can recalibrate a fold's starting slab variance. Effective fold settings are preserved with raw fit objects and summarized in validation-checks.json.

Scalar independent-gene binary-set models use native direct cross-products; joint traits retain projected preparation. Model-specific n-1 intercept-only residual normalization, original-unit coefficients and existing priors are preserved. This is not the standalone qgg-compatible BLR convention n. The earlier one-chain 14.21-second native probe and preparation-included general runs remain earlier evidence with different controls/binaries. qgg's vb/vg R/native positional mismatch is not reproduced, and no matched qgg speedup or external MAGMA/FUMA comparison is claimed.

## Timing and memory definitions

Gene-association statistics are caller supplied: their upstream computation time is unmeasured (NA/null), not zero. Cached fixture loading is measured separately in timing.json and excluded from operation time. Preparation includes set mapping/deduplication, R checks and native summaries or general projection. Native BLR time includes checks, sampling, posterior accumulation and original-unit coefficient back-transformation. Subsequent posterior conversion/diagnostics are recorded separately. PoPS also records held-out prediction. Total call time includes remaining interface/assembly overhead; components need not sum exactly to the total. Output serialization, process startup and input loading are separate.

PoPS validation records every fold's preparation and fitting and explicitly labels the aggregate as a sum of 22 models. A process-local exit trace records completed-fold progress; it does not change inputs or controls, and its small overhead is included in workflow wall time. The same tracing is not used in full-data fits.

Memory is the measured R worker PID lifetime high-water working set observed at 100 ms polling, including loading and serialization. Sampled private memory is a separate process.json field. Peaks cannot be assigned to preparation or fitting independently and are not incremental sampler allocation. The 1,536 MiB preparation guard is not total process memory. The workstation is the same Intel Core i7-1365U/approximately 16 GiB Windows host; runs are serial, but the workstation is not an isolated performance laboratory. Each row is one execution, with no timing uncertainty estimate. Completion is not convergence evidence.

## Provenance and diagnostics

Small stage records, phase tables, session information and effective settings are preserved here. Original workers and orchestration scripts are under scripts/. Input, source and loaded native binary hashes are in provenance.json. Raw fits and gene-level predictions remain in ignored local work storage; retained-artifacts.json records their locations, sizes and hashes. The earlier evidence collection is preserved unchanged.

The first full-data PoPS harness attempt passed raw saved memberships to strict gsea and stopped before sampling because some memberships were duplicated. The successful harness applies the same unique/intersection preparation already performed by the magma/pops interfaces. Its 0.01-second failed attempt is retained under diagnostics and excluded from all successful clocks. The diagnostic caused no model or dataset change.

These are gsuite ecosystem benchmarks, not external MAGMA, FUMA, UK Biobank scoring, biological concordance or general scaling measurements. No commit, push or publication is part of the measurement harness.
