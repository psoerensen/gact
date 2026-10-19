# Complete regional VEGAS workflow

Client date: 2026-10-05. Complete cached European CAD/T2D website workload:
**17,862 retained Genesplus genes**, **503 reference
samples**, and **9,985,562 physical reference markers**.
Both trait inputs are the original qgg-checked marker-level GWAS data. No gene
X2 sums from the previous benchmark are used as input to the new analysis.
Gene names, marker counts and reconstructed statistic sums are checked against
the previous complete qgg output. Download, installation, QC and reference-Glist
creation remain separate, just as in the original qgg analysis timing.

## Measured workflow

Report marker-set eigenvalue preparation separately from the VEGAS fit. For
the same 17,862 marker sets and two traits, the existing clocks give:

| Operation | R call | Seconds | Minutes |
| --- | --- | ---: | ---: |
| Compute marker-set LD spectra | `gbase::gprep(G, task="eigen", sets=sets)` | **643.36** | **10.72** |
| Convert spectra to the reusable VEGAS object | `glma::prepare_vegas(G)` | **79.52** | **1.33** |
| Initial VEGAS fit for both traits using the prepared object | `glma::vegas(stat, prepared=prepared)` | **97.69** | **1.63** |
| Repeat VEGAS fit using the restored prepared object | `glma::vegas(stat, prepared=restored)` | **71.44** | **1.19** |

The marker-set preparation clock includes reference validation, reader setup,
genotype decoding, LD products and eigenvalue computation. Summing the native
per-gene clocks gives **102.10 s** decoding, **160.49 s** matrix products and
**108.40 s** in the eigensolver; these are subphases of the 643.36 s preparation
clock, not additional workflow time. The spectra are computed once and shared
by both traits. Converting them with `prepare_vegas()` performs no new LD
computation or eigendecomposition.

The fit clock measures the complete `glma::vegas()` call with supplied spectra,
including statistic validation, gene-statistic construction, tail evaluation
and result assembly. The repeat fit is a separate reuse measurement; its
17.26 s prepared-object loading time is excluded from the 71.44 s fit clock.
Input loading, marker alignment and summary-statistic conversion remain
separate in the detailed table below.

The **15.00-minute** figure is the complete analysis, including preparation,
adapter work, alignment and checkpoint overhead. It must not be labelled as
the VEGAS fit time. The historical 8.30-hour qgg result covers its combined
workflow; separate qgg eigenvalue and fit clocks are not available from that
measurement. [Phase summary](phase-summary.csv) is derived from the preserved
stage clocks; no new benchmark run or change to the original evidence archive
was required.

| Stage | Seconds |
| --- | ---: |
| input_loading | 110.45 |
| common_marker_and_gene_alignment | 8.59 |
| regional_eigen_preparation | 643.36 |
| vegas_spectrum_adapter | 79.52 |
| summary_alignment_and_statistic_conversion | 12.71 |
| vegas_evaluation | 97.69 |
| prepared_RDS_loading | 17.26 |
| repeated_vegas_evaluation | 71.44 |

Complete analysis wall clock, from common-marker mapping through both traits'
VEGAS outputs: **899.94 seconds
(15.00 minutes)**. This includes reference validation,
full regional BED eigenvalue preparation, conversion to glma's serializable
prepared object, alignment of marker statistics and reconstruction of gene
quadratic sums, and native tail evaluation. It also includes intermediate
checkpoint/cleanup overhead outside the individually timed stages. Input-file
loading and the additional persistence/reuse verification are shown separately.
The earlier qgg analysis took 29,876.32 seconds (8.30 hours); it is dated evidence,
not a contemporaneous repeated comparison, and its tail algorithm differs.
The original website used a multi-trait P-value interface. This benchmark
converts those marker P-values to sqrt(qchisq(p, df=1, lower.tail=FALSE)), applying
the same 1e-10 P-value floor. Their signs do not affect the diagonal quadratic
statistic; all marker identities and alleles are still checked explicitly.

The restored prepared object returned identical per-trait results without any
LD read or eigendecomposition. It contains **3,192,403 union
markers**, in the exact marker/allele order required for each trait. Each gene
retains its complete original membership. The smaller sample-space Gram matrix
preserves the full nonzero spectrum for genes larger than the sample count.
Local BED decoding uses four threads; portable products and Eigen spectra remain
serial. No OpenBLAS change or global library/thread preference was introduced.

## Remaining numerical differences

| Trait | Available | Unavailable | Maximum X2 difference | Median P difference | Maximum P difference | Changed Bonferroni calls |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| GWAS1 | 17,862 | 0 | 9.095e-13 | 2.417e-08 | 0.01838 | 0 |
| GWAS2 | 17,862 | 0 | 4.547e-13 | 3.098e-08 | 0.01962 | 0 |

glma retains its existing default saddlepoint evaluator, with explicit reporting
of its moment fallback and unavailable cases. These are approximate tails.
qgg uses a different evaluator and floors eigenvalues below 1e-7. The adapter
does not apply that floor: it rejects material negative eigenvalues, discards
rounding-level negatives under the existing PSD control, omits exact zeros and
retains all positive eigenvalues. Consequently the timing gain does not prove
identical P-values, calibrated inference, or a drop-in replacement for qgg in the
paper. The complete fitted outputs and algorithm labels are archived.

As a separate benchmark oracle, qgg's installed original tail evaluator was
applied to the **newly reconstructed** statistic sums and new BED spectra with
qgg's same 1e-7 eigenvalue floor. This took **50.63 seconds** for
both traits. Median absolute P difference from the historical qgg result was
**1.97065e-15**; maximum difference was
**1.71336e-05**. Changed calls: **0**
at 0.05 and **0** at 0.05/17,862. This comparison
isolates LD/statistic reconstruction from the different glma tail approximation.
It uses qgg's internal function only as a benchmark oracle; production glma
does not depend on qgg. Its P-values are archived separately.

## Validation and preserved evidence

Canonical isolated gbase/glma source builds, Windows binary dependency audits,
clean loads and focused test suites passed. The gsea VEGAS re-export/handoff
check passed. New checks compare independent-marker analytic chi-square tails,
constant-weight rescaling, overlapping and reordered sets, exact set membership,
RDS reuse after removing temporary reference files, and cross-file family/marker
identity checks including unused components. This is partial validation, not a
fresh complete R CMD check or an ecosystem qualification campaign.

The redundant BedPanel construction in Glist validation was removed. The same
ordered FAM identity, unique individual-ID and globally unique physical-marker
checks now use metadata already parsed by the checked BedMatrix reader. The
subsequent native regional preparation still builds one checked panel; no
validation rule was waived. Current shared validation/preparation times are in
[shared-eigen-times.csv](shared-eigen-times.csv).

[All per-gene clocks](per-gene-times.csv), [50 slowest genes](slowest-genes.csv),
[stage times](stage-times.csv), and [qgg numerical comparison](qgg-comparison.csv).
Peak working set for the full worker, including input loading:
**5.88 GiB**; sampled peak private allocation:
**12.43 GiB**. These are process metrics,
not the size of each gene matrix.

The [source and evidence archive](source-and-evidence.zip) contains current
gbase/glma R/native integration source, the approved gmat regional helper,
scripts, spectra, serializable preparation, fitted outputs, timing and session
metadata. SHA256: `e0acae22de201d89cbe749cf3b5847744d0522ff9da027c188b3e0f3788ccfe3`. Large input caches and earlier evidence remain in
their original locations; no previous timing report was replaced.
Publication status updated 2026-10-05: implementation committed and pushed in
gsuite (`98681ba`) and gmat (`36dc14b`). The validated public documentation
snapshot was synchronized through the guarded updater and committed and pushed
in gtools (`4428541`). The evidence ZIP remains the original measured-source
snapshot; its contents and checksum are unchanged.
