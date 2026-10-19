# LD-swap timings on the frozen simulated-human example

3 serial paired repeats with alternating off/on order. One thread for OpenMP and BLAS,
the same inputs, seeds, priors and iteration counts in each pair. Swaps off uses probability
zero; swaps on uses probability 0.05, minimum r2 0.8, maximum 50 friends and up to one proposal
per selected sweep. Regional fits have 100 burn-in and 1000 retained sweeps (seed 10);
the genome fit has 100 sweeps and no burn-in (seed 20261005). Each fit has one chain.
The frozen inputs contain 5,000 people and 50,000 physical markers, of which
37,991 pass the existing QC; summary-statistics n=4998 represents regression
residual degrees of freedom. Download, simulation and QC are reused unchanged.

| Workflow | Swaps | Median operation | Range | Median peak RSS | Attempts / accepted |
| --- | --- | ---: | ---: | ---: | ---: |
| gmap, four dense regions | off | 15.890 s | 14.920–22.950 s | 0.493 GiB | 0 / 0 |
| gmap, four dense regions | 5% | 20.060 s | 15.770–20.730 s | 0.491 GiB | 0 / 0 |
| gbayes, sparse genome | off | 6.610 s | 6.200–9.740 s | 0.814 GiB | 0 / 0 |
| gbayes, sparse genome | 5% | 8.170 s | 6.080–8.200 s | 0.814 GiB | 0 / 0 |

Runtime varied across repeats, including the LD phases that precede swapping.
Whole-operation median differences therefore do not isolate the swap cost.
paired.csv preserves each matched comparison; region-summary.csv includes native-fit phases.
Do not infer a speed improvement or a precise percentage overhead from these short pairs.

No swap proposals were attempted. The untimed LD audit found no off-diagonal pair
at r2 >= 0.8 in any of the four full regional matrices or the stored genome LD.
Maximum regional r2 values are 0.00444–0.00602; the genome maximum is 0.00631.
These weak-LD inputs measure enabling/preparation overhead, rather than the cost
or mixing benefit of accepted swaps. A mixing assessment needs a workload with
strongly correlated alternative markers. The audit CSV files preserve the counts.

| Region | Markers | Off median | On median | On attempts / accepted |
| --- | ---: | ---: | ---: | ---: |
| 19:55551662 | 850 | 3.136 s | 3.945 s | 0 / 0 |
| 11:935343 | 785 | 2.629 s | 3.427 s | 0 / 0 |
| 13:92838400 | 1091 | 3.605 s | 4.590 s | 0 / 0 |
| 11:92296987 | 1399 | 4.880 s | 5.796 s | 0 / 0 |

Operation time excludes process startup, package/input loading and RDS serialization.
Regional operation time includes reading the existing BED, computing full regional LD,
preparing its dense likelihood, native fitting, posterior tables and credible groups.
Genome sparse LD preparation is reused and excluded. Neither workflow recomputes genotype LD
for swap proposals. The regional native-fit phase includes one-time LD-friend preparation
when enabled; it is not a separate pure Gibbs-loop measurement. Phase totals omit some
R wrapper work, so use whole-operation timings for workflow comparisons.
Process peak RSS includes package loading and is measured from the Windows R process.
Full per-repeat timings, memory, phases and swap counts are in the adjacent CSV files.

The four known causal-marker PIPs range from 1.000 to 1.000 across both workflows and settings.
Fixed seeds reproduce estimates and swap counts exactly within each setting across repeats.
Enabling swaps changes RNG consumption and can change other PIPs. Full PIP comparisons are
recorded separately. Native single-chain diagnostic statuses and ESS/MCSE values are retained
in diagnostics.csv and the fit RDS files. These short examples do not establish a mixing
improvement or convergence, especially the 100-sweep genome fit. Swaps remain opt-in.

These are fresh controls against the current installed packages and the corrected native
reverse-edge check. Earlier qgg comparisons and source/evidence archives are unchanged.
This report does not replace historical timings or assert numerical equivalence to qgg.

Reproduce from the gsuite checkout with:
`powershell.exe -NoProfile -ExecutionPolicy Bypass -File benchmarks/genomic_20k_50k/time_ld_swap.ps1`.
Then run `Rscript --vanilla benchmarks/genomic_20k_50k/summarize_ld_swap.R <runs.csv> <new-report-directory>`.
The scripts reuse the preserved frozen example and prepared sparse-LD cache under build/benchmarks.
Run `Rscript --vanilla benchmarks/genomic_20k_50k/audit_ld_swap_example.R <report-directory>`
after the timing processes finish to reproduce the untimed LD audit.
The build/evidence directory records source commits/diffs, installed-package SHA256 hashes,
benchmark scripts and the run inventory; each run includes hardware, input hashes, session
information, effective controls and full fitted results. No download, QC or simulation is repeated.
`source-and-evidence.zip` preserves the six runs, this report and its CSV files,
the exact installed gbase/gmap/gbayes packages and DLLs, source provenance and
benchmark scripts, frozen statistics and corrected BIM copy. Original genotype
files and the large sparse-LD payload remain in place; their input references
and genotype checksums are retained. The ZIP integrity check and adjacent SHA256
file support later recovery without relying only on the build directory. The
archive is retained locally and is excluded from the source commit because it
contains installed packages and DLLs. Its SHA256 is
`e7044ff1788eb7d8f1a7e174186cd66bf32d1a53304b52c7dd6361e8e4780920`.
Retained local measurements and logs are in `build/benchmarks/final-summaries/`,
with the study prefix `ld-swap-20261006-110710--`.

Measurements were collected before committing the LD-swap update. Source
commits and diffs in the preserved archive identify the measured development
state; the source commit includes the scripts, report and compact timing/LD
evidence, rather than generated binaries, logs or full fitted-result RDS files.
