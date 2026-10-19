# Simulated-human fine-mapping: qgg and standalone R packages

This initial report measured older installed gmap/gbayes 0.0.1 wrappers. The
[current-package comparison](../manuscript-simulated-finemapping-current-20261005/README.md)
uses refreshed 0.1.0 packages and supersedes its primary model timing table.
The original source-and-evidence archive remains unchanged. Repository sources
inside that archive provide context; installed binaries and session logs identify
the runtime actually measured.

Client date: 2026-10-05. The [published gact example](https://psoerensen.github.io/gact/Document/Finemapping_bayesian_linear_regression_simulated_data.html)
has regional BayesC mapping and genome-wide summary BayesC followed by LD-based
posterior sets. We reused its frozen genotype files, QC and simulated statistics:
5,000 genotype samples, 50,000 physical markers, 37,991 retained markers and four
causal-centered regions. After intersecting with the retained statistics, regional
sizes are 850, 785, 1,091 and 1,399 markers. qgg reports 4,998 in its stat$n column
(regression residual degrees of freedom); this is preserved, not relabelled as
the actual number of people. The frozen simulation seed is 20261005.

## Timings

Fresh one-thread descriptive executions; input loading and RDS serialization are
outside operation clocks. Regional standalone analysis includes full regional LD
construction/persistence, eigendecomposition, fitting and native posterior sets.
Both genome fits reuse already constructed genome-wide LD. Download, QC, phenotype
simulation and marginal GWAS were reused identically, not rerun for either engine.

| Operation | qgg | Standalone packages |
| --- | ---: | ---: |
| Four regional BayesC fits and posterior sets | 15.87 s | 73.77 s |
| Genome-wide BayesC fit | 13.51 s | 7.47 s |
| Genome posterior-set post-processing | 4.71 s | 12.37 s using **qgg** |
| Standalone selected-reference preparation | — | 7.50 s |
| 1,000-marker-window LD construction (standalone includes sampled audit) | 190.53 s | 198.31 s |

The earlier qgg LD operation was 95.44 s and is preserved as historical evidence;
the table uses a fresh qgg repeat after the standalone runs, without concurrent
timing workers. Standalone LD retains r2=0, a 1,000-marker
distance window and no base-pair window. The saved float32 band and reopened CSR
were compared at 32 sampled markers (the exact audit is saved in the Glist).
This is a sampled identity check, not an exhaustive LD qualification.
The maximum sampled correlation difference was 1.86e-09.
Maximum standalone process working set across these stages: 0.81 GiB.
Process clocks/private allocation and operation clocks are separately recorded
in [stage-measurements.csv](stage-measurements.csv).

## Model matching and limits

Regional: both use BayesC, 1,000 post-burn sweeps plus 100 burn-in sweeps, one
chain, inclusion prior .01, h2-based starting residual/slab values, four-degree
scaled inverse-chi-square variance priors and learning of slab variance,
residual variance and inclusion probability. Standalone gmap uses qgg's supplied
working wy/ww and response normalization, full regional signed LD, and retained
positive eigenvalue mass 1 with fallback disabled. It does not apply the default
99.5% truncation. It still performs a native eigendecomposition; qgg's mcmc route
uses its dense correlation matrix directly. qgg samples/returns genetic variance
separately; gmap has no equivalent independent genetic-variance draw.
Posterior-set algorithms are different marginal-PIP grouping procedures; they
are not component credible sets or a coverage calibration experiment.

Genome: both use BayesC, 100 sweeps, no burn-in, one chain, inclusion prior .001,
matching normalized initial variance/prior scales and learning flags. Standalone
gbayes uses its standardized summary likelihood and fixed reference-LD adjustment
coefficient .9. qgg uses the supplied variable ww/wy and its marker-specific
LD-neighborhood adjustment involving genetic variance. These are corresponding
workflows, not the same likelihood or transition kernel. Effect scales and RNG
trajectories differ. No identical-model speedup or inferential equivalence is
claimed. These short single-chain examples do not establish convergence;
standalone regional fits explicitly remain provisional because chain agreement
cannot be assessed from one chain.

The standalone public R interface currently lacks LD grouping of a genome-wide
**summary** gbayes fit. gmap_bed() supports saved individual-BED fits only.
The genome posterior-set comparison therefore uses the same qgg getCredibleSets()
helper on each implementation's PIPs. Its standalone-side clock is a shared qgg
postprocessor, not a gsuite implementation; the result records this explicitly.

## Numerical outputs

| Region | Markers | qgg causal PIP | gmap causal PIP | All-marker PIP correlation |
| --- | ---: | ---: | ---: | ---: |
| 19:55551662 | 850 | 1.0000 | 1.0000 | 0.99899 |
| 11:935343 | 785 | 1.0000 | 1.0000 | 0.99808 |
| 13:92838400 | 1091 | 1.0000 | 1.0000 | 0.99779 |
| 11:92296987 | 1399 | 1.0000 | 1.0000 | 0.99635 |

| Causal marker, genome fit | qgg PIP | gbayes PIP |
| --- | ---: | ---: |
| 19:55551662 | 1.0000 | 1.0000 |
| 11:935343 | 1.0000 | 1.0000 |
| 13:92838400 | 1.0000 | 1.0000 |
| 11:92296987 | 1.0000 | 1.0000 |

Across all genome markers, PIP correlation was 0.97751;
mean absolute PIP difference was 0.00038088.
The common qgg genome postprocessor selected the same four singleton sets.

These output comparisons characterize this example; they are not tests of equal
seeded draws. See [regional-comparison.csv](regional-comparison.csv),
[genome-causal-comparison.csv](genome-causal-comparison.csv) and the saved fits for
posterior estimates, diagnostics, hyperparameters and posterior-set memberships.

## Input correction and evidence

The original tutorial BIM contains 13 invalid negative/nonfinite
genetic-map distances. gmat rejected it before computation. An owned BIM copy
replaces only those genetic-map entries with zero. BED bytes, FAM, physical
coordinates, markers and alleles are unchanged. These model/LD calls do not use
genetic-map distance. The original BIM and failed preparation remain intact.
The detailed changes are archived as genetic-map-adjustments.csv.

This run uses installed development packages. It is not a fresh R CMD check,
calibration study, release qualification or gact workflow modification. The
benchmark scripts are local uncommitted work. Existing manuscript evidence is
preserved. [Source and evidence](source-and-evidence.zip) contains the exact
driver, worker, installed qgg function captures, fitted objects, timing/session
logs, metadata adjustment and installed DLLs. SHA256: `8d4e70cc8d6d4c4c2c6022ab717d78053f5d5dbdd2e009587021102178fba00c`.
Large immutable genotype/LD inputs remain in their original cache; no original
data or earlier timing report was overwritten.
