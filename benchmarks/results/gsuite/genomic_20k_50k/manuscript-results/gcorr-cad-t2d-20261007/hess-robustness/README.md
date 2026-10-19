# Robust compact HESS input and conditioning workflow

This records the earlier R 0.1.2/native 0.1.1 validation step. The subsequent
[signed-moment and automatic-rank confirmation](../hess-moment-uncertainty/README.md) describes the current
R 0.1.3/native 0.1.2 workflow. Reproduction below targets the earlier installation;
use the new verification driver for the current version.

The standalone gcorr R adapter is version 0.1.2. The native gcorr 0.1.1
estimator and strict conditional Gaussian uncertainty contract are unchanged.
This update follows the [focused investigation](../hess-investigation/README.md).

## Changes

- Compact preparation defaults to a common rank cap of
  `min(50, floor(min(N)/(2*number_of_regions)))`. The conservative total-rank
  bound is at most half the smallest supplied N. This is checked before
  reference validation/reading and uses no trait projections. A budget below
  one component per region stops with an actionable error.
- Explicit numeric caps remain reproducible. A high conservative rank bound
  warns before preparation, and actual rank/N above 0.5 warns at fitting.
  Actual ranks, correction denominators and moment amplification are reported.
- Explicit trait-named genomic-control factors scale unchanged compact moments
  at fit time. Each factor multiplies Z-squared and cross moments receive the
  geometric mean. Nonunit factors require source or sensitivity provenance.
  Every fit starts from the original moments; repeated fits cannot accumulate
  scaling. Sample sizes are unchanged and no factor is inferred from LDSC.
- Inference diagnostics identify negative regions by trait, total h2 outside
  [0,1], indefinite regional genetic matrices, invalid residual matrices and
  missing-SE counts. Raw estimates and native uncertainty guards are preserved.

Stronger truncation can discard genetic signal. Better conditioning does not
establish unbiased estimates or calibration. A hypothetical genomic-control
factor is a sensitivity assumption, not a reconstructed correction.

## Focused validation

The installed isolated package passed the focused `test-hess-projection.R`
file, including four small deterministic cases. Checks cover the original
BED/dense-LD oracle and RDS/allele behavior, automatic lower-rank selection
against a dense joint solve, named factor scaling against directly rescaled
Z scores, unchanged raw projections/repeated fits, rejected invalid factors,
and negative-region diagnostics while retaining valid marginal uncertainty for
the other trait. This is focused validation, not full package qualification
or statistical calibration.

The [saved-moment verification driver](../../../verify_hess_robust.R) reuses the
original 1,702-region compact genome object. It fits no new regional
eigensystems and acquires no data. It confirms:

- All saved cap-20 pair estimates, missingness, uncertainty and reasons are
  unchanged to the stated numerical comparison tolerance.
- The conservative automatic cap is 11, with total-rank bound 18,722 and
  smallest supplied N 40,743.1524049809. Its worst-case denominator is about
  0.5405. This is a preflight bound, not a measured rank after decomposition.
- The diagnostics recover one negative T2D region and 1,352 negative CAD
  regions, preserving their raw estimates and the missing-SE qualification.
- The explicit CAD Z-squared factor 1.18 reproduces the independent aggregate
  moment calculation. Its provenance states that this is a sensitivity and
  the complete cohort correction history remains unverified.

The original cap-20 estimates and timings remain the measured reference. No
automatic-cap genome projection or new timing campaign has been run.

Saved checks:

- [Automatic rank preflight](automatic-genome-rank-preflight.csv)
- [Saved cap-20 conditioning](saved-cap20-conditioning.csv)
- [Saved cap-20 inference summary](saved-cap20-inference.csv)
- [Negative regions by trait](saved-cap20-negative-regions.csv)
- [Explicit CAD 1.18 sensitivity](saved-cap20-cad-1.18-sensitivity.csv)
- [Recorded scaling provenance](scaling-provenance.csv)

## Reproduction

From the gsuite checkout with the updated gcorr package installed first in
`build/task-packages/library`:

```text
Rscript --vanilla benchmarks/genomic_20k_50k/verify_hess_robust.R
  build/benchmarks/manuscript-20261007-070903/gcorr_genome-4/output.rds
  build/gact-revision-cleanup/hess-robust
```

The command is shown on separate lines for readability; pass all three
arguments in one shell invocation. Its private cache is the existing corrected
workflow output and contains compact moments, identities and spectra, without
persistent LD/eigenvectors. Public outputs contain only aggregate diagnostics
and regional blocker identities.

The gcorr guide and existing CAD/T2D tutorial describe the new interface and
retain the original result/timing qualification. Their rendering does not run
analyses. The tutorial rendered successfully and its 44 local references
passed inspection. The gsuite allowlisted publisher checked 7,471 internal
references across 92 HTML pages with zero failures and synced the local gtools
snapshot, preserving a backup. The final isolated install and installed help
match the updated API. Commits/pushes remain pending.
