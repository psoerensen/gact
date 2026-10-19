# Signed moment uncertainty and automatic-rank CAD/T2D confirmation

This subsequent analysis uses gcorr R 0.1.3 and native gcorr 0.1.2 (ABI 0).
It preserves the original cap-20 estimates/timings and follows the earlier
[input/conditioning update](../hess-robustness/README.md).

## Statistical contract

Compact fits default to `uncertainty="moment"`. The Gaussian fourth-moment
formulas are evaluated as a plug-in covariance estimator at signed local/total
moment estimates, then every signed block is propagated through the joint
correction. A negative h2 alone does not remove SEs elsewhere. Only finite
nonnegative resulting scalar variances receive SEs. Negative variances stay
flagged; no points or covariance matrices are clipped/repaired. Correlation
inference also requires positive h2, an interior correlation and a PSD estimated
local/total sampling covariance. Heritability P values use the upper normal tail.

This is approximate uncertainty conditional on the supplied reference LD,
sample sizes, overlap and retained subspaces. It has no broad calibration,
unbiasedness or exact upstream variance-parity claim. It does not assert a
Gaussian model with negative genetic variance. `uncertainty="strict"` preserves
the previous conditional admissibility checks; legacy prepared-LD R/native
interfaces retain that policy. Existing native symbols and struct layouts remain.

## CAD statistic provenance

Database metadata identifies GWAS2 as `CARDIoGRAMplusC4D.txt.gz`, Nikpay 2015,
PMID 26343387. The audited cached Z agrees with stored b/SE. The current
`getMarkerStat()` extraction uses b/SE and introduces no genomic-control
adjustment. Ingestion/extraction source inspection found no GC correction in
that path; it does not establish the original release file's complete history.

The [primary paper](https://discovery.dundee.ac.uk/ws/files/7268983/CAD1000G2015.pdf)
documents study-level GC before meta-analysis and a second correction afterward
(Online Methods). Figure 4 gives the final additive meta-analysis factor 1.18.
Multiplying CAD Z-squared by 1.18 is therefore a **partial final-layer undo /
sensitivity**, conditional on the stored release using that layer. It cannot
reconstruct cohort-specific factors and weights from the aggregate statistics.
No uncorrected source or full cohort correction history has been verified.
The default scale remains unchanged; no factor is inferred from LDSC.

## Saved cap-20 check

Strict re-fitting reproduces the saved cap-20 point estimates to 1e-12 comparison
tolerance. Changing only the uncertainty policy preserves every point exactly.
New SE availability:

| component | quantity | trait1 | trait2 | rows | available_points | available_se | negative_variance |
| --- | --- | --- | --- | --- | --- | --- | --- |
| region | cov | GWAS1 | GWAS2 | 1702 | 1702 | 1701 | 1 |
| region | h2 | GWAS1 |  | 1702 | 1702 | 1701 | 1 |
| region | h2 | GWAS2 |  | 1702 | 1702 | 1476 | 226 |
| region | rg | GWAS1 | GWAS2 | 1702 | 350 | 259 | 0 |
| total | cov | GWAS1 | GWAS2 | 1 | 1 | 1 | 0 |
| total | h2 | GWAS1 |  | 1 | 1 | 1 | 0 |
| total | h2 | GWAS2 |  | 1 | 1 | 1 | 0 |
| total | rg | GWAS1 | GWAS2 | 1 | 0 | 0 | 0 |

Saved totals under the moment policy:

| quantity | trait1 | trait2 | estimate | standard_error |
| --- | --- | --- | --- | --- |
| h2 | GWAS1 |  | 0.680078 | 0.0147577 |
| h2 | GWAS2 |  | -0.187175 | 0.0360318 |
| cov | GWAS1 | GWAS2 | 0.106081 | 0.00500492 |
| rg | GWAS1 | GWAS2 | unavailable | unavailable |

The negative CAD total remains visible and genome rg remains undefined.
Obtaining marginal SEs does not validate the statistic/reference scale.

## Fresh automatic-rank confirmation

The same 7,221,923 website-common GWAS1/GWAS2/GWAS6-available markers and 1,702
disjoint regions are retained. GWAS6 supplies availability only. There are no
new MAF, palindrome, indel, missingness, HWE or duplicate-position exclusions.
Extracted `stat$n` remains 49,071.2732949614 (T2D) and 40,743.1524049809 (CAD).
The database convention is Ncase*Ncontrol/(Ncase+Ncontrol), without a factor-four
or total-participant substitution. User-supplied sample overlap remains zero,
without an independent overlap audit. Four preparation threads are requested.

Preparation streams the existing 503-person EUR BED reference and discards
LD/eigenvectors after each region. Only compact projections are saved privately.
The automatic cap is 11, selected from N and region count before trait projection.
Actual conditioning:

| trait | sample_size | retained_rank | rank_fraction | correction_denominator | moment_amplification |
| --- | --- | --- | --- | --- | --- |
| GWAS1 | 49071.2732949614 | 18696 | 0.380996838774096 | 0.619003161225904 | 1.61550063495565 |
| GWAS2 | 40743.1524049809 | 18696 | 0.458874654915371 | 0.541125345084629 | 1.84800066950035 |

Unadjusted automatic-rank totals:

| quantity | trait1 | trait2 | estimate | standard_error |
| --- | --- | --- | --- | --- |
| h2 | GWAS1 |  | 0.273841 | 0.00797956 |
| h2 | GWAS2 |  | 0.0131915 | 0.00890362 |
| cov | GWAS1 | GWAS2 | 0.0733679 | 0.00367812 |
| rg | GWAS1 | GWAS2 | 1.2207 | unavailable |

**Unadjusted genome rg is 1.2207, outside [-1,1]. Its delta SE is unavailable;
this raw ratio does not provide interpretable biological correlation inference.**

Partial CAD final-meta GC adjustment (T2D factor 1, CAD factor 1.18):

| quantity | trait1 | trait2 | estimate | standard_error |
| --- | --- | --- | --- | --- |
| h2 | GWAS1 |  | 0.273841 | 0.00797956 |
| h2 | GWAS2 |  | 0.168206 | 0.0100067 |
| cov | GWAS1 | GWAS2 | 0.079698 | 0.00370513 |
| rg | GWAS1 | GWAS2 | 0.371345 | 0.0190625 |

Unadjusted SE availability:

| component | quantity | trait1 | trait2 | rows | available_points | available_se | negative_variance |
| --- | --- | --- | --- | --- | --- | --- | --- |
| region | cov | GWAS1 | GWAS2 | 1702 | 1702 | 1685 | 17 |
| region | h2 | GWAS1 |  | 1702 | 1702 | 1670 | 32 |
| region | h2 | GWAS2 |  | 1702 | 1702 | 1468 | 234 |
| region | rg | GWAS1 | GWAS2 | 1702 | 585 | 327 | 0 |
| total | cov | GWAS1 | GWAS2 | 1 | 1 | 1 | 0 |
| total | h2 | GWAS1 |  | 1 | 1 | 1 | 0 |
| total | h2 | GWAS2 |  | 1 | 1 | 1 | 0 |
| total | rg | GWAS1 | GWAS2 | 1 | 1 | 0 | 0 |

A smaller cap improves conditioning but discards more signal. These are
conditional estimates under the database scale and finite EUR reference, not
validated liability-scale heritabilities or calibrated biological rg inference.
The partial GC scenario does not resolve the complete source-scale limitation.
The final CAD h2 adjustment matches an independent aggregate-Q formula;
T2D points are unchanged and covariance scales by sqrt(1.18). Regional h2/cov
contributions sum to their corresponding totals (all checks within 1e-11).

## Timing, memory and validation

Automatic preparation: 1473.86 s. Unadjusted pair fit: 24.03 s. Partial-GC pair fit: 24.95 s.

Peak R working set: 4.946 GiB; sampled peak private memory: 14.334 GiB (100 ms sampling).
Memory covers the entire R process, including saved-moment verification and
input loading. Monitoring began during the first cache read. Input reads and
compact serialization are separate phase clocks. This is one confirmation run;
documentation rendering and focused native compilation ran concurrently, so it is not a controlled speed
comparison. Earlier measured timing evidence is unchanged.

Three focused native tests passed, including an independently assembled dense
Jacobian/signed fourth-moment oracle. The installed R wrapper passed 43 focused
assertions across five small cases. They check point preservation, valid SEs at
negative local estimates, directional P values and confined negative-variance
failures. This is bounded numerical validation, not full package qualification.

## Reproduction and files

[Verification driver](../../../verify_hess_moment_uncertainty.R), run from the
gsuite checkout after installing gcorr in the existing isolated library:

```text
Rscript --vanilla benchmarks/genomic_20k_50k/verify_hess_moment_uncertainty.R build/gact-revision-cleanup/hess-moment
```

The driver reuses the existing corrected input and cap-20 caches. It acquires
no data and creates no genome LD/eigenvector cache. Memory monitoring is an
external process measurement; the driver records preparation/fitting clocks.

- [analysis-input-audit.csv](analysis-input-audit.csv)
- [automatic-moment-conditioning.csv](automatic-moment-conditioning.csv)
- [automatic-moment-estimates.csv](automatic-moment-estimates.csv)
- [automatic-moment-traits.csv](automatic-moment-traits.csv)
- [automatic-momentinference-counts.csv](automatic-momentinference-counts.csv)
- [automatic-partial-meta-gc-conditioning.csv](automatic-partial-meta-gc-conditioning.csv)
- [automatic-partial-meta-gc-estimates.csv](automatic-partial-meta-gc-estimates.csv)
- [automatic-partial-meta-gc-traits.csv](automatic-partial-meta-gc-traits.csv)
- [automatic-partial-meta-gcinference-counts.csv](automatic-partial-meta-gcinference-counts.csv)
- [partial-scaling-checks.csv](partial-scaling-checks.csv)
- [phases.csv](phases.csv)
- [process-memory.csv](process-memory.csv)
- [regional-phases.csv](regional-phases.csv)
- [saved-cap20-moment-conditioning.csv](saved-cap20-moment-conditioning.csv)
- [saved-cap20-moment-estimates.csv](saved-cap20-moment-estimates.csv)
- [saved-cap20-moment-traits.csv](saved-cap20-moment-traits.csv)
- [saved-cap20-momentinference-counts.csv](saved-cap20-momentinference-counts.csv)
- [scaling-provenance.csv](scaling-provenance.csv)
