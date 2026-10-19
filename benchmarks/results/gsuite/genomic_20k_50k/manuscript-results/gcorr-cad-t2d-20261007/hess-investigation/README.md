# Focused CAD/T2D HESS investigation

The compact projection is numerically correct on the three checked real regions.
The negative CAD genome estimate is a shortfall of projected energy relative to
the unit-noise correction, amplified by the high retained rank. This investigation
does not replace the reported results with a chosen positive answer. Production
gcorr 0.1.1, `stat$n`, the website marker selection and the published fits are
unchanged.

## Independent regional calculation

Regions were chosen by marker count, without using either trait's effects:
the smallest region, the closest to 300 markers, and the closest to 3,000 markers.
An independent SNP-major PLINK decoder reads the existing BED offsets, centers
observed allele-1 dosages, imputes missing centered values as zero and normalizes
each genotype column to unit norm. Base R SVD computes marker-space projections;
the two small regions also use a direct dense-LD eigendecomposition. These are
bounded numerical oracles, not another qgg comparison or calibration campaign.

| Region | Markers | Native space | Maximum absolute moment difference |
| --- | ---: | --- | ---: |
| chr1:178944309-178954470 | 25 | markers | 3.93e-19 |
| chr10:134336329-134494673 | 451 | markers | 2.17e-18 |
| chr6:47311898-48391125 | 3,000 | samples | 5.53e-18 |

The dense marker-space checks differ by at most 1.52e-18. Both marginal
quadratic forms and the cross-trait moment agree. See
[regional-projection-oracle.csv](regional-projection-oracle.csv).

The cached CAD Z scores equal the database effect/SE values exactly for all
3,476 selected markers. For the 3,182 informative finite P-value comparisons,
the median ratio `(b/SE)^2 / Z(P)^2` is 1.00000022; its range is
0.99998655-1.00001123, consistent with rounding. There is no evidence here
that the cache substituted a different P-based or effect/SE statistic. See
[cad-z-and-p-provenance.csv](cad-z-and-p-provenance.csv).

## Why CAD is negative

Write `Q = N * sum(q_region)` and let `K` be total retained rank. The current
native joint total is `(Q-K)/(N-K)`. The independently reconstructed totals
agree with the saved native fit within 1e-10.

| Trait | Extracted N | K | Projected energy Q | Q-K | h2 |
| --- | ---: | ---: | ---: | ---: | ---: |
| T2D | 49,071.2733 | 33,971 | 44,240.3604 | 10,269.3604 | 0.680078 |
| CAD | 40,743.1524 | 33,971 | 32,703.4203 | -1,267.5797 | -0.187175 |

CAD energy is about 3.7% below the subtracted baseline. The correction
denominator is only 0.166216, so that shortfall is amplified. Rescaling N alone
cannot change the sign while N remains above K: it rescales q and K/N together.
The diagnostic factor-four N scenario therefore stays negative; it is not used
as a replacement information-size convention. See
[genome-moment-decomposition.csv](genome-moment-decomposition.csv) and
[sample-size-sensitivity.csv](sample-size-sensitivity.csv).

The pinned upstream code's different row/column rank convention is too small
to explain the result. At the same measured moments and ranks its CAD total
would be -0.188073 instead of -0.187175; T2D would be 0.681014 instead of
0.680078. This is an algebraic matrix-convention comparison, not execution
of the complete upstream pipeline. See
[joint-rank-convention-sensitivity.csv](joint-rank-convention-sensitivity.csv).

## Genomic control and reference uncertainty

The [CAD study](https://discovery.dundee.ac.uk/ws/files/7268983/CAD1000G2015.pdf)
documents genomic control within cohorts and again after meta-analysis;
the additive meta-analysis factor was 1.18. Its participants were 77% European,
whereas this analysis uses an EUR-only reference. The downloaded database Z/P
agreement is consistent with a common corrected scale, but does not reconstruct
the complete cohort-specific correction history or the GWAS LD covariance.

[HESS guidance](https://huwenboshi.github.io/hess/local_hsqg/) recommends undoing
a known applied genomic-control factor explicitly. At these fixed moments,
multiplying CAD Z-squared by 1.05 gives a total of 0.05428; a factor of 1.10
gives 0.29573; the paper's 1.18 factor gives 0.68206. The break-even factor is
1.03876. These are sensitivity calculations, not established corrected
heritabilities. They show why a small variance-scale change matters here.
An inverse LDSC intercept is not automatically the study's genomic-control
factor. No reinflation is silently applied to either trait. See
[reinflation-sensitivity.csv](reinflation-sensitivity.csv).

One deterministic disjoint 251/252-person reference split also shows the
limitation of assuming the fitted LD is the GWAS noise covariance. In the
three regions, held-out inverse-LD noise energy per retained component is
0.980, 0.935 and 0.942, rather than exactly one. This single reference-only
experiment is illustrative; it is neither a GWAS correction nor a sampling
distribution. Correlation is undefined for columns monomorphic in either half:
only this split experiment omits one of 451 and 19 of 3,000 markers. All original
markers remain in the native/direct projection checks and published fits. See
[reference-split-noise-diagnostic.csv](reference-split-noise-diagnostic.csv).

Genomic control, finite-reference error, ancestry mismatch and truncated signal
can contribute. This bounded investigation does not isolate their respective
genome-wide contributions or establish a binary-trait liability scale.

## Rank policy and uncertainty

The 20-component cap was chosen to pass `K < N`, which establishes numerical
admission but not good statistical conditioning. The inspected upstream HESS
code defaults to `gwse_thres=2`, giving a cap
`min(50, ceiling(N/(2 * number_of_regions)))`: 15 for T2D and 12 for CAD across
1,702 regions. A common two-trait policy would use 12 here. This would reduce
the worst-case rank ratio to about one-half, but changes the retained signal
and does not guarantee a positive or unbiased total. No genome-wide cap-12
run was made. The direct three-region cap comparisons show that estimates
can change sign and are not monotonic as components are added. See
[upstream-rank-policy-comparison.csv](upstream-rank-policy-comparison.csv) and
[regional-truncation-sensitivity.csv](regional-truncation-sensitivity.csv).

The native conditional uncertainty rule requires every regional heritability
to be nonnegative and the total to be within [0,1]. T2D has one negative region
(chr15:20001200-21131604, h2 = -1.667e-5);
this blocks the joint conditional standard errors for all 1,702 regions,
even though 1,701 raw marginal plug-in moment variances are nonnegative.
CAD has 1,352 negative regions and a negative total. The bivariate gate also
checks regional genetic and total residual covariance validity. Missing
standard errors therefore reflect explicit model-admissibility checks; they
are not missing projected information or a numerical decomposition failure.
No guard was relaxed or point estimate clipped. See
[uncertainty-gates.csv](uncertainty-gates.csv).

## Practical conclusion

Keep LDSC and chromosome-set LDSC as the current CAD/T2D correlation examples.
For HESS, the next analysis needs an explicit genomic-control/input-scale
contract and a preselected, better-conditioned rank policy, while preserving
the extracted information sizes and website marker selection. Any conditional
uncertainty policy allowing raw negative local estimates needs its own declared
statistical contract. The present raw HESS results remain diagnostic and do not
support paper-level regional correlation inference.

## Reproduction and scope

The driver is [investigate_gcorr_hess.R](../../../investigate_gcorr_hess.R).
From the gsuite root, using the existing corrected input and saved genome run:

```powershell
& 'C:/Program Files/R/R-4.4.1/bin/x64/Rscript.exe' --vanilla `
  benchmarks/genomic_20k_50k/investigate_gcorr_hess.R `
  build/benchmarks/manuscript-20261007-003431/gcorr_inputs-4/output.rds `
  build/benchmarks/manuscript-20261007-070903/gcorr_genome-4/output.rds `
  build/gact-revision-cleanup/hess-investigation `
  build/benchmarks/manuscript-20261004-143135/database-1/database/hsa.0.0.1/gstat/GWAS2.txt.gz
```

The successful focused run reused a private reduced descriptor of the same
3,476 markers to avoid repeatedly deserializing the full input. The driver
can create that descriptor from the original input. Original caches and prior
evidence remain intact. No data were downloaded or ingested, no whole-genome
analysis was rerun, and no LD/eigenvector files were stored. Only aggregate
CSV diagnostics are published; the private descriptor, raw summary statistics,
genotypes and logs stay under ignored `build/`.

Sources inspected: native gcorr e807f74, gsuite f7fae9f, gact 16aa3a5, and
upstream HESS commit c7a23308096059ffa4df1071eab8a3b77c3c1ad6 (`src/estimation.py`).
This is focused numerical agreement and input/model diagnosis, not broad
statistical calibration. The native and R gcorr interfaces are unchanged.
