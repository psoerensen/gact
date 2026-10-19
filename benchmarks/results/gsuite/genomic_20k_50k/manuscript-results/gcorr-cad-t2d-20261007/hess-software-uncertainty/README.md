# HESS software uncertainty and saved CAD/T2D comparison

gcorr R 0.1.4 / native 0.1.3 (ABI 0) changes compact fits to
`uncertainty="hess"`. The original software's full-regional-LD-rank h2 and
covariance variance systems replace the previous conditional SE default.
Genome rg uses the original fixed-contribution delete-one-region jackknife,
centered on the full-sample ratio. Diagonal-plus-rank-one variance solves use
linear storage and arithmetic. No LD/eigenvectors were prepared again.

Point estimates are unchanged from the previous
[automatic-rank confirmation](../hess-moment-uncertainty/README.md), within
1e-12 of its decimal exports. The same 7,221,923 website-common markers,
1,702 regions, 18,696 retained components (cap 11), 503-person EUR reference,
zero-overlap assumption and extracted `stat$n` remain. N is
49,071.2732949614 for T2D and 40,743.1524049809 for CAD, using the database's
Ncase*Ncontrol/(Ncase+Ncontrol) convention. No additional QC mask or sample-size
transformation was applied. Binary/liability-scale validity is not established.

## Current totals

Unchanged statistic scale:

| quantity | trait1 | trait2 | estimate | standard_error |
| --- | --- | --- | --- | --- |
| h2 | GWAS1 |  | 0.273841 | 0.018285 |
| h2 | GWAS2 |  | 0.0131915 | 0.0283976 |
| cov | GWAS1 | GWAS2 | 0.0733679 | 0.0149692 |
| rg | GWAS1 | GWAS2 | 1.2207 | 0.332661 |

Raw CAD h2 is close to zero and poorly determined. The raw rg above one has
only a descriptive jackknife SE; Wald Z/P are withheld. It is not an
interpretable biological correlation. Larger SEs correct the uncertainty
convention; they do not recover the source statistic scale or discarded signal.

Partial final-meta CAD factor 1.18, explicitly a sensitivity:

| quantity | trait1 | trait2 | estimate | standard_error |
| --- | --- | --- | --- | --- |
| h2 | GWAS1 |  | 0.273841 | 0.018285 |
| h2 | GWAS2 |  | 0.168206 | 0.0241534 |
| cov | GWAS1 | GWAS2 | 0.079698 | 0.0127451 |
| rg | GWAS1 | GWAS2 | 0.371345 | 0.0209133 |

The [CAD paper](https://discovery.dundee.ac.uk/ws/files/7268983/CAD1000G2015.pdf)
documents cohort-level and final-meta genomic control. The 1.18 factor can
only address the final layer, conditional on this database release using it.
No uncorrected source or complete earlier correction history was verified.
LD-reference and ancestry mismatch and eigenvalue truncation remain possible
contributors. A finite SE does not validate this sensitivity as a full correction.

## Independent original-code comparison

The pinned [HESS estimation source](https://github.com/huwenboshi/hess/blob/c7a23308096059ffa4df1071eab8a3b77c3c1ad6/src/estimation.py)
has LF SHA256 d53aa10941746d7af8c164404aabf4bd3ea8e9a614078fb1e5d6bafd6b9cb5ee.
Its unmodified second-stage function definitions were executed on the saved
moments, with loaders/projection sums adapted to those moments and Python-2
`xrange`/removed NumPy aliases supplied. The same cap 11 was held fixed in
both traits. This is a second-stage comparison, not an upstream CLI run with
its original MAF filter or separate automatic rank selections.

The original column-rank point convention also gives raw CAD h2 0.0131514,
compared with native 0.0131915. Its partially adjusted h2 is 0.168218,
compared with native 0.168206. That discrepancy does not explain the low
heritability. Native retains its row-rank moment equations and their
independent mean-correction checks; the point convention is explicit.

Original variance functions were then executed at the **same native points**:

| scenario | quantity | trait | max_variance_difference | negative_variances |
| --- | --- | --- | --- | --- |
| raw | h2 | GWAS1 | 2.6456544711896037e-19 | 4 |
| raw | h2 | GWAS2 | 8.16063305120315e-19 | 4 |
| raw | cov | GWAS1/GWAS2 | 9.595295105615121e-20 | 0 |
| partial_gc | h2 | GWAS1 | 2.6456544711896037e-19 | 4 |
| partial_gc | h2 | GWAS2 | 6.347453148486913e-19 | 4 |
| partial_gc | cov | GWAS1/GWAS2 | 8.605325348511658e-20 | 0 |

All finite nonnegative regional variance estimates agree to 1e-12 comparison
tolerance; negative variances stay unavailable. Total variance is the signed
sum, including negative regional entries. The independently resummed R
delete-region ratios and jackknife SE agree within 1e-11.

This explicitly named software convention reproduces its variance equations,
including recorded full LD rank and overlap/residual polynomial terms. It is
distinct from the conditional Gaussian fourth moments. No joint sampling
covariance or regional rg SE is supplied. `uncertainty="moment"` retains the
previous signed conditional covariance/delta route; `"strict"` retains its
admissibility checks. Native legacy/prepared-LD defaults remain strict.

Two differences from original inference are deliberate: negative h2 uses
the directional upper tail instead of absolute-Z one-sided P, and invalid
jackknife ratios remain missing and invalidate the SE instead of being
replaced with zero. Out-of-range full rg has no Wald inference. None of these
changes clips point estimates or repairs a covariance matrix. Statistical
calibration under finite-reference, binary-scale or GC uncertainty remains
unestablished.

## Checks and timing

Three focused native checks passed (7.92 s), including independently assembled
dense full-rank variance systems with overlap and explicit jackknife deletion.
The installed R wrapper passed 50 focused assertions across six cases.
The current refits took 45.42 s (raw) and 37.23 s (partial GC). Website rendering
and installation verification overlapped part of this run; these clocks are
not a controlled comparison with earlier timings. Numerical variance/jackknife
computation remains linear in region count; validation of marker identities
and R/native conversion remains part of fit time. Existing preparation/timing
evidence is preserved. No new peak-memory measurement was taken.

## Reproduction

From gsuite with the current isolated installation, run
[the saved-moment driver](../../../verify_hess_software_uncertainty.R):

```text
Rscript --vanilla benchmarks/genomic_20k_50k/verify_hess_software_uncertainty.R build/gact-revision-cleanup/hess-software
```

[The independent Python comparison](../../../compare_hess_software_reference.py)
requires NumPy, pandas, SciPy, a caller-supplied copy of the pinned original
source, the exported saved moment CSV and the native fit outputs:

```text
python benchmarks/genomic_20k_50k/compare_hess_software_reference.py PINNED_ESTIMATION_PY MOMENTS_CSV OUTPUT_DIR
```

Original source and dependencies are validation inputs, not package/runtime or
public-render dependencies. No upstream source is vendored in the package.
This is focused numerical comparison, not full-package or broad statistical
qualification. Changes remain local; no commit, push or deployment is included.

- [fit-timings.csv](fit-timings.csv)
- [native-totals.csv](native-totals.csv)
- [original-point-convention-totals.csv](original-point-convention-totals.csv)
- [partial_gc-estimates.csv](partial_gc-estimates.csv)
- [partial_gc-jackknife.csv](partial_gc-jackknife.csv)
- [point-convention-comparison.csv](point-convention-comparison.csv)
- [raw-estimates.csv](raw-estimates.csv)
- [raw-jackknife.csv](raw-jackknife.csv)
- [upstream-variance-parity.csv](upstream-variance-parity.csv)
