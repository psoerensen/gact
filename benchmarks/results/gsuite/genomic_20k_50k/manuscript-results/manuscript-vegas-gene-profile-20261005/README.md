# Website VEGAS gene profile

Client date: 2026-10-05. Complete cached European CAD/T2D workload from the
[gact website example](https://psoerensen.github.io/gact/Document/Gene_analysis_vegas.html):
18,279 raw Genesplus sets; **17,862 retained genes**, identical memberships and
marker counts to the saved qgg output after its common-marker QC and minimum size.
Reference: **503 samples; 9,985,562 physical markers**. No new download, change to
gene boundaries, marker trimming, or sparse-LD construction.

## Observed bottlenecks

The initial gbase helper rebuilt a 10-million-marker lookup for every gene.
Its 20-gene pilot took **14.32 s**.
After the single-pass fix, all 17,862 genes matched in **5.51 s**.
For one real gene, qgg's ID-based genotype call took
**0.75 s**; its offset-based call returned identical
genotypes in less than the 0.01 s timer resolution. This proves an avoidable
identifier lookup cost; it does not assign an exact fraction of the historical
whole qgg run to that cost.

| Stage | Seconds |
| --- | ---: |
| validation | 232.74 |
| matching | 6.90 |
| panel_setup | 209.86 |
| decode | 163.88 |
| matrix_product | 231.94 |
| eigen | 175.49 |
| other / cleanup | 4.97 |
| two-trait tails | 81.22 |

Whole `gprep(task="eigen")`: **1025.78 s**.
Two-trait qgg tail calculation using the new spectra and cached exact X2 sums:
**81.22 s**. Eigen preparation plus tails:
**1107.00 s**.
Historical qgg VEGAS total: **29,876.32 s**, from the prior completed run,
not rerun concurrently. These scopes differ: the new numerical benchmark reuses
the prior X2 sums; the former qgg total included its marker/statistic alignment.

The full reference audit and native panel setup are measured separately.
The validator constructs a full BedPanel, followed by another full BedPanel
for the eigenspectra. This duplicate setup is a concrete next optimization
target. No validation rule was removed during this run.

## Per-gene evidence

[All per-gene clocks](per-gene-times.csv), [50 slowest genes](slowest-genes.csv),
[size groups](gene-size-summary.csv), and [stage totals](stage-times.csv).
Median native time per gene: **12.92 ms**;
95th percentile: **141.25 ms**;
99th percentile: **246.72 ms**.
Native total excludes shared reference validation/panel setup and tail evaluation.
Product time includes matrix allocation, transpose where needed, the product,
Gram symmetry construction and genotype buffer release. Decomposition time
includes eigenvalues-only solving, zero padding and sorting. Other/cleanup is
the remainder of the complete gprep wall clock.

| Gene | Markers | Decode s | Product s | Eigen s | Native total s |
| --- | ---: | ---: | ---: | ---: | ---: |
| CSMD1 (ENSG00000183117) | 10,385 | 0.4033 | 0.9856 | 0.0593 | 1.4492 |
| RBFOX1 (ENSG00000078328) | 9,246 | 0.3453 | 0.8716 | 0.0701 | 1.2885 |
| PTPRD (ENSG00000153707) | 7,241 | 0.2364 | 0.6445 | 0.0810 | 0.9622 |
| CNTNAP2 (ENSG00000174469) | 6,435 | 0.2886 | 0.5928 | 0.0641 | 0.9460 |
| WWOX (ENSG00000186153) | 4,669 | 0.1839 | 0.4176 | 0.0736 | 0.6777 |
| CTNNA2 (ENSG00000066032) | 3,644 | 0.1795 | 0.4148 | 0.0732 | 0.6679 |
| MACROD2 (ENSG00000172264) | 4,808 | 0.1482 | 0.4347 | 0.0727 | 0.6578 |
| DLG2 (ENSG00000150672) | 4,483 | 0.1664 | 0.4222 | 0.0657 | 0.6546 |

![Gene timing profile](gene-profile.png)

## Numerical checks and limits

All gene identities/counts agree with the former qgg output. Maximum absolute
P-value difference: **1.27423e-05**; maximum log10
difference: **1.4226e-05**.
Median absolute difference is 1.97e-15; the 95th percentile is 1.65e-13.
No classifications changed at P < 0.05 or P < 0.05/17,862.
All returned spectra and P-values were finite; P-values lay within [0,1].
Both calculations use qgg's original 1e-7 eigenvalue floor and its installed
tail implementation; the latter is accessed only as a benchmark oracle.
This does not add a consumer to glma VEGAS.

A focused follow-up recomputed qgg's original marker-space eigendecomposition
for five genes among the largest discrepancies, comparing their CAD tails.
Maximum eigenvalue difference was **1.42109e-12**.
The same qgg tail routine turns these roundoff-level spectrum differences into
P-value differences of approximately 1e-5 in a few cases. Even repeating qgg's
original spectrum gave a maximum difference of
**1.22611e-05** from its saved CAD result.
This supports numerical sensitivity in the tail evaluation; it does not establish
an exact tail-error bound or bitwise equivalence. The numerical audit is included
in the ZIP, and the maximum discrepancy is retained above rather than hidden.

The first attempt stopped on a strict symmetry check: general matrix products
accumulate the two triangles independently. The corrected Gram construction
retains one triangle and reflects it, checking the raw discrepancy against a
dot-product rounding-error bound before solving. The solver tolerance is unchanged. Maximum measured raw relative
triangle discrepancy: **2.27596e-14**.
The deterministic regression checks agreement with full R `cor()` eigenvalues
under a tighter 1e-18 symmetry control. All **125 focused gbase checks passed**.

Four local BED decoder threads; portable product and Eigen spectrum solver
remain serial. No OpenBLAS installation/change. Per-region steady-clock timing
and a 20 ms R profiler were enabled. One development run, no confidence intervals;
Windows process peak working set: **6.36 GiB**;
sampled peak private allocation: **12.19 GiB**,
including input loading and setup. Those metrics are not regional matrix sizes.

The reusable raw cache is `build/gbase-eigen/website-genes`; the 787 MB lean input
RDS is intentionally excluded from the evidence ZIP. The ZIP contains exact
changed source, scripts, spectra, timings, session information, and failed-run
logs. Its SHA256 is `30dad04a9460548782104f741f607785f806802755ffad5ce467447098b3aa12`.
Publication status updated 2026-10-05: the completed implementation is committed
and pushed in gsuite (`98681ba`) and gmat (`36dc14b`), with public documentation
in gtools (`4428541`). This earlier profiling evidence and its ZIP remain
unchanged; they describe the source snapshot measured during that run.
