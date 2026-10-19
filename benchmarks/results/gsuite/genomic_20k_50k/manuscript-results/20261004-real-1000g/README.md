# Real-1000G-derived standalone workflow benchmark

Measured on 4 October 2026 at one thread. These replace the artificial-LD
analysis timings as the current benchmark. Earlier results, input caches and
failed diagnostic runs are preserved. Database installation evidence is reused
from the [completed database audit](../20261004/README.md).

## Measurements

| Workload | Analysis wall time (s) | Peak R working set (MiB) | Output RDS (bytes) |
| --- | ---: | ---: | ---: |
| VEGAS: two traits, 1,000 sets of 274 markers | 44.78 | 274.7 | 290,906 |
| Joint Bayesian pathways: two traits, three pathways | 0.10 | 102.9 | 12,589 |
| BayesR mcmc-eigen: 20 regions per trait, 100 markers each | 195.65 | 157.3 | 1,195,338 |
| Total scoring: 2,000 recipients, 5,360 markers | 1.03 | 157.2 | 121,159 |
| Three pathway scores: same recipients | 2.41 | 157.8 | 180,722 |

Database download/install took 544.20 s, with peak R working set 2.677 GiB.
Installed files total 3,261,689,414 bytes; retained download archives total
440,170,912 bytes. Saving the database descriptor added 12.48 s and 52,318,792
bytes. All seven component categories passed the recorded file/checksum audit.

Auxiliary real-reference setup took 113.86 s for import, 112.32 s for founder
generation, 404.91 s for pedigree/BED conversion and 75.51 s for phenotype,
GWAS and LD preparation. The BED is 240,440,003 bytes; BIM/FAM add 2,875,943
bytes. The sparse-LD resource is 11,695,176 bytes across its five files.
Fixed scoring-weight construction took 0.01 s and produced a 72,116-byte RDS.
These setup stages are excluded from the five analysis times.

[Measurements](measurements.csv) also record whole-process elapsed time,
serialization time, sampled private memory, stage-directory size and warnings.
Analysis wall time surrounds the public operation; dependency loading and
input loading are excluded. Mapping statistic preparation and recipient Glist
preparation occur before the analysis clock. Peak memory covers the complete R
process lifetime, including loading and serialization, using Windows
`PeakWorkingSet64` for the actual R PID. Private memory is sampled every 100 ms.
MiB/GiB use powers of 1024. Each workload was executed once successfully; these
are illustrative measurements without repeated-run uncertainty estimates.

## Inputs and models

The cached official 1000G Phase 3 chr22 VCF, GRCh37 genetic map and population
panel were read without modification. All 503 EUR donors were imported over
22:20,000,000–22,200,000. Public gsim reference-conditioned simulation produced
19,998 founders and two offspring, preserving real physical coordinates and
map distances. They are simulated individuals derived from real reference LD.
The original 48,088 sites became 5,360 markers after cohort MAF >= .01 filtering.
This is a 20K-person workload, **not a post-QC 50K-marker workload**.

Two continuous A/B traits have sparse fixed effects (causal probability .02,
h2 .3/.4), with separate 6,000-person GWAS cohorts, an 18,000-person LD reference
and 2,000 held-out score recipients. These are performance proxies, not measured
CAD/T2D analyses. Input hashes and the exact window membership are recorded in
[the input manifest](input-manifest.csv) and [gene-set manifest](gene-set-manifest.csv).

The 1,000 overlapping contiguous windows each contain 274 markers, matching
the installed gact Genesplus **raw** median before reference/GWAS filtering.
They are size-controlled timing sets, not biological gene annotations or an
estimate of the post-QC gene-size distribution. The three pathway memberships
are synthetic. LD stores signed float32 correlations with no r2 pruning and a
273-marker distance window, giving complete pairwise coverage within each
274-marker timing set and each 100-marker mapping region.

VEGAS includes eigenspectrum preparation and both-trait saddlepoint evaluation.
Float32 LD required `relative_psd_tolerance=1e-7`; the default 1e-10 rejected
rounding-level negative eigenvalues. Mapping uses `verification_tolerance=1e-7`.
These numerical controls are recorded; no broad PSD projection was performed.
The global window-truncated LD is not guaranteed positive semidefinite.

The joint gsea BayesC model uses independent gene covariance, identity trait
and effect covariance, null/A/B/both pattern probabilities .985/.005/.005/.005,
four seeds 11/29/47/71, 200 burn-in and 800 sampling sweeps. This is the current
standalone joint Bayesian pathway workflow corresponding to the intended
BLR-MAGMA workload; no equivalence to the original qgg model is asserted.

gmap BayesR uses residual variance 1, effect variance .04, mixture probabilities
.7/.1/.1/.1 and variance multipliers 0/.01/.1/1. Four chains use 500 burn-in and
700 sampling sweeps. Retained eigen masses .995/.99/.95 are attempted under
the default convergence-selection rule. The measured time includes retries.
Only **3 of 40 region–trait fits were accepted; 37 are provisional**.
[Diagnostics](mapping-diagnostics.csv) and [all eigen attempts](mapping-eigen-attempts.csv)
are preserved. This measures the bounded computation, not a converged
inferential fine-mapping analysis.

Scoring uses deterministic sparse simulation-effect weights, normalized to unit
L2 norm, standardized genotype coding and unscaled output. Three pathway scores
use deduplicated unions of the timing sets. These are fixed weights, not fitted
posterior weights. An auxiliary genome-wide BayesC fit on window-truncated LD
failed with finite-range arithmetic overflow; its output was excluded. No
prediction-quality claim is made from the scoring measurements.

## Provenance, validation and interpretation

Hardware: Intel Core i7-1365U, 10 physical/12 logical cores, 16,794,288,128 bytes
RAM, Windows 11 Enterprise 26200. Software: R 4.4.1, GCC 13.2.0, existing Release
native builds, OpenMP configured to one thread and BLAS to one thread.
[Hardware](hardware.json), package versions/library paths, native ABI/version
reports and installed DLL hashes are stored alongside this report. gsim is
0.16.0. The standalone packages are the point of departure for future qgg/gact
integration; this is not a timing of the original gact/qgg implementations.
The reused database was installed through gact's published download interface.

Dimensions and finite VEGAS, joint pathway, mapping and score outputs passed
the recorded lightweight checks. This does not certify model calibration or
MCMC convergence. Four-thread comparisons were performed only for the earlier
artificial fixture; their speedups must not be assigned to this real-LD run.
The older gbase build manifest differs from four installed DLL hashes; exact
measured hashes are recorded, with this source-attribution limitation retained.

A linear gene-count extrapolation at fixed 274-marker set size gives about
818.5 s (13.6 min) for 18,279 sets from the measured 44.78 s/1,000 sets. This is
an estimate, not a genome-wide measurement: spectrum cost depends strongly on
set size, the large-gene tail is unmeasured, and these windows overlap.
Scoring here uses 5,360 markers and cannot be extrapolated by gene count alone.

For Reviewer 2, a useful bounded comparison is one official standalone MAGMA
gene-analysis run using identical A-trait summaries, markers, 18K reference
samples and 274-marker sets, followed by a matched VEGAS run. Fix executable
version, gene model, allele/QC alignment and computation boundaries first.
Such a comparison measures workload performance, not equivalence to the joint
Bayesian pathway model. No standalone MAGMA comparison has been executed.

## Reproduction and preserved runs

Drivers: [measurement launcher](../../run_manuscript.ps1),
[real-reference worker](../../real_ld_stage.R), [summary](../../summarize_real_ld.R).
The [protocol](../../real-1000g-protocol.md) records the final controls.
Run the launcher with `-RealLD` for a fresh workflow, or use `-ReferenceRds`,
`-CohortRds`, `-FixtureRds` and `-StartStage` to reuse completed stages.
Every launch creates a new directory and refuses overwrites.

Successful stages were reused across these preserved directories under
`build/benchmarks/`: import `manuscript-20261004-160735`, founders/BED
`manuscript-20261004-161106`, preparation `manuscript-20261004-162110`, VEGAS/joint
`manuscript-20261004-162349`, mapping `manuscript-20261004-162613`, and
weights/scoring `manuscript-20261004-163047`. Failed stages remain diagnostic
evidence and are excluded from the table. The file inventory includes all
successful stage resources; source VCF/map/panel and the database are separate.
