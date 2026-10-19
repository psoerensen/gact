# Current gsuite revision benchmarks

Measured on 19 October 2026 using existing installed gsuite packages and frozen inputs. Raw inputs, fit objects and score matrices remain in ignored local storage or the original working directories; no individual data are part of this evidence collection.

## Gene set and single-model BLR fits

The existing hsa.0.0.1 database supplies 17,955 saved gene scores and 1,324 Reactome sets with at least 50 annotated genes. Inputs are reused without recomputing VEGAS or rebuilding the database. The scalar/joint gene-set calls use the current gsea interface with Reactome instead of the KEGG annotation in the gene-set tutorial. The cached gene universe is distinct from the 17,862-gene timed VEGAS workflow.

Single-trait BayesC uses pi=.01, 5,000 sampling sweeps, 1,000 burn-in and two chains. Joint BayesC uses pi=.01, 1,000 sampling sweeps, 200 burn-in and two chains. Bayesian PoPS uses pi=.001, 5,000 sampling sweeps, 1,000 burn-in, two chains and 22 chromosome-held-out folds. Ridge PoPS applies its p<.05 feature screen. All these runs request one thread. Current defaults otherwise apply.

The full Bayesian Reactome workload requires an explicit 1,536 MiB preparation budget; the 512 MiB default guard attempt is preserved as a failed diagnostic. A rejected 2 GiB integer-control attempt and a rejected integer chromosome-label PoPS attempt are also excluded from the successful timings. The successful harness supplies character chromosome labels. A preliminary scoring attempt stopped before score application because its assertion checked the full-file sample count rather than selected study IDs; the corrected harness checks the selected ID vector and every returned score row. This is a resource/type adaptation, not a reduction in genes, sets, sweeps or folds.

| Operation | Operation seconds | Process peak GiB | Output RDS bytes |
| --- | --- | --- | --- |
| Saved gene evidence and annotation loading | 17.95 | 0.99 | 1382110 |
| Single-trait Reactome BayesC | 50.37 | 0.51 | 413867 |
| Joint CAD/T2D Reactome BayesC | 91.36 | 0.52 | 481890 |
| Native BayesR BLR from precomputed summaries; one chain | 14.21 | 0.02 (native only) | Not RDS |

## Single-model BLR from precomputed summaries

Sparse feature construction took 0.46 s; centered/scaled X'X, X'y and y'y preparation took 0.20 s. The existing native gsuite gbayes engine then fitted one BayesR model in **14.21 s**, from these prepared summaries. This uses 17,955 cached gene scores, 1,324 Reactome features, one chain, 1,000 burn-in and 5,000 sampling sweeps, with qgg's actual BayesR starting priors and learning settings. No gene-association statistics were recomputed. Native summary loading (0.0153 s), binding (0.0037 s) and output writing (0.0424 s) are separate. The 0.02 GiB peak covers only the native worker, not R preparation.

The [full report and R-interface implementation plan](summary-blr-qgg-settings/README.md) document qgg matching, output/trajectory differences and the native benchmark driver. The sufficient-statistic sampler exists, but a direct XX/Xy/yy/n R wrapper is not yet exposed. The earlier 162.87-second general gsea call is retained as preparation-included evidence with different controls; it is excluded from the main BLR benchmark. These are not timings of a single unchanged binary or a matched qgg speedup comparison. The separately timed 17,862-gene VEGAS workload has different inputs and must not be added to the BLR timing.

## Separate chromosome-held-out validation

| Validation workflow | Wall seconds | Process peak GiB |
| --- | --- | --- |
| Ridge, 22 chromosome-held-out folds | 3.72 | 0.70 |
| BayesC BLR, 22 chromosome-held-out folds | 977.20 | 1.61 |
| BayesR BLR, 22 chromosome-held-out folds | 1908.66 | 1.63 |

These optional validation timings are not timings of one full-data BLR fit. Bayesian validation uses two chains per fold with the same sweep and burn-in counts as above.

## Scoring grid

The frozen 20K/50K artificial-LD fixture was generated earlier with binomial/block-copy genotypes and gsim effects. The current gbase/gscore workers select 500/2,000 held-out people and 10K/50K available SNPs, using fixed effects. They separately time total scores or 3/20 synthetic sets with 500 SNPs each and zero/50% adjacent overlap. The set-score union, not the available panel size, determines the number of contributing SNPs. There is no model fitting or new simulation in the scoring clock.

Eight total-score and 32 set-score cells are single executions, with one or four requested threads. Inputs are loaded and recipient descriptors are prepared before the operation clock. Whole-process peaks include these steps and final score serialization. Standardized genotype coding and unscaled completed scores are used. These are computational probes, not tests of prediction quality or uncertainty estimates. Do not infer general asymptotic scaling from this bounded grid.

| People | Available SNPs | Sets | Adjacent overlap | Union SNPs | Seconds 1 / 4 threads | Process MiB 1 / 4 threads |
| --- | --- | --- | --- | --- | --- | --- |
| 500 | 10000 | 0 | 0 | 10000 | 1.77 / 0.72 | 129.3 / 129.6 |
| 500 | 10000 | 3 | 0 | 1500 | 0.47 / 0.28 | 128.3 / 129.4 |
| 500 | 10000 | 3 | 0.5 | 1000 | 0.36 / 0.21 | 129.7 / 129.4 |
| 500 | 10000 | 20 | 0 | 10000 | 3.57 / 1.11 | 134.3 / 135.4 |
| 500 | 10000 | 20 | 0.5 | 5250 | 1.66 / 0.64 | 132.2 / 136.1 |
| 500 | 50000 | 0 | 0 | 50000 | 13.62 / 3.77 | 155.7 / 155.8 |
| 500 | 50000 | 3 | 0 | 1500 | 0.53 / 0.28 | 155.7 / 156.2 |
| 500 | 50000 | 3 | 0.5 | 1000 | 0.42 / 0.27 | 155.5 / 155.7 |
| 500 | 50000 | 20 | 0 | 10000 | 4.17 / 1.19 | 158.8 / 159.0 |
| 500 | 50000 | 20 | 0.5 | 5250 | 2.10 / 0.75 | 167.9 / 159.2 |
| 2000 | 10000 | 0 | 0 | 10000 | 1.75 / 0.68 | 129.6 / 125.4 |
| 2000 | 10000 | 3 | 0 | 1500 | 0.48 / 0.25 | 129.5 / 129.7 |
| 2000 | 10000 | 3 | 0.5 | 1000 | 0.37 / 0.19 | 129.9 / 130.1 |
| 2000 | 10000 | 20 | 0 | 10000 | 3.83 / 1.04 | 134.8 / 135.4 |
| 2000 | 10000 | 20 | 0.5 | 5250 | 1.67 / 0.66 | 134.0 / 132.2 |
| 2000 | 50000 | 0 | 0 | 50000 | 13.67 / 3.78 | 156.2 / 155.7 |
| 2000 | 50000 | 3 | 0 | 1500 | 0.59 / 0.31 | 156.3 / 156.3 |
| 2000 | 50000 | 3 | 0.5 | 1000 | 0.45 / 0.27 | 156.3 / 155.8 |
| 2000 | 50000 | 20 | 0 | 10000 | 4.75 / 1.35 | 159.0 / 166.4 |
| 2000 | 50000 | 20 | 0.5 | 5250 | 2.39 / 0.73 | 159.2 / 159.7 |

## Remaining tutorial modules

| Operation | Status | Seconds | Process peak GiB |
| --- | --- | --- | --- |
| Chromosome-component Gaussian LOCO | completed | 90.26 | 0.35 |
| Chromosome-component prior-selected LOCO | completed | 67.07 | 0.37 |
| Two-chain genome-wide BayesC from cached LD | completed | 2.91 | 0.67 |
| Two-chain regional BayesC | completed | 6.13 | 0.34 |

The two-chain regional run supplies the existing synthetic assembly label. The genome-wide run reuses validated cached sparse LD with an explicit max_ld_entries=1e8 budget, matching the earlier preparation protocol. Initial assembly/entry-guard failures are preserved as diagnostics and excluded from timings.

## Interpretation and provenance

Hardware, package/session records, settings, process counters, input hashes and output sizes are preserved. The preparation budget is not a total process-memory limit. Peaks are Windows R PID lifetime high-water working sets observed at 100 ms polling; sampled private memory is a different field. Model and scoring operation clocks exclude startup, input loading and final serialization. The separately labeled saved-input preparation clock includes loading the installed descriptor, saved scores and annotations; model-specific design preparation remains within each model clock. The workers run serially but the workstation is not an isolated performance laboratory.

Returned Bayesian fits do not establish convergence; retain their diagnostics before biological interpretation. Raw model objects remain at recorded local paths. The benchmark is restricted to the gsuite ecosystem; it does not run the external MAGMA executable or establish FUMA/MAGMA concordance. Scripts retain measured-path assumptions and must be configured deliberately before future use.
