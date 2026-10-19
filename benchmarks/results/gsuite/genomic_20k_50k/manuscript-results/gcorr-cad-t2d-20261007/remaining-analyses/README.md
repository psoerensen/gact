# Corrected remaining CAD/T2D gcorr analyses, 7 October 2026

All methods use the same 7,221,923-marker website availability mask and extracted stat$n.
Stored LD scores are reused by GNOVA and SumHer. HESS streams BED regions and discards LD/eigenvectors.
These recorded HESS runs use K <= 50 for chr22 and an explicit K <= 20 for the genome. They predate the automatic rank policy in the compact R adapter.
The lower genome cap is set from input rank admission, not fitted trait signal; it changes the spectral model.
Results remain conditional on the database information-size convention and statistical restrictions described in the tutorial.
The subsequent [original HESS software uncertainty comparison](../hess-software-uncertainty/README.md)
refits the saved automatic-cap projections, corrects h2/covariance SEs and adds
the genome rg region jackknife. It preserves the point estimates and keeps the
partial CAD genomic-control adjustment explicitly labelled as a sensitivity.
The unchanged native/R gcorr 0.1.1 baseline is preserved in the neighboring 6 October source-and-evidence.zip; this record supplies corrected workers, inputs audits and results.

| method | quantity | trait1 | trait2 | estimate | standard_error |
| --- | --- | --- | --- | --- | --- |
| ldsc | h2 | GWAS1 |  | 0.469755 | 0.0246419 |
| ldsc | h2 | GWAS2 |  | 0.310603 | 0.0192085 |
| ldsc | cov | GWAS1 | GWAS2 | 0.146811 | 0.0119046 |
| ldsc | rg | GWAS1 | GWAS2 | 0.384344 | 0.028498 |
| gnova | h2 | GWAS1 |  | 0.498234 | 0.0276627 |
| gnova | h2 | GWAS2 |  | 0.110616 | 0.0191808 |
| gnova | cov | GWAS1 | GWAS2 | 0.141822 | 0.0105536 |
| gnova | rg | GWAS1 | GWAS2 | 0.604114 | 0.0586559 |
| sumher | h2 | GWAS1 |  | 0.487726 | 0.0289254 |
| sumher | h2 | GWAS2 |  | 0.361629 | 0.0229225 |
| sumher | cov | GWAS1 | GWAS2 | 0.16261 | 0.0129741 |
| sumher | rg | GWAS1 | GWAS2 | 0.387194 | 0.0272514 |

Completed corrected HESS/rho-HESS totals:

| analysis | quantity | trait1 | trait2 | estimate | standard_error |
| --- | --- | --- | --- | --- | --- |
| gcorr_chr22 | h2 | GWAS1 |  | 0.00387591 | unavailable |
| gcorr_chr22 | h2 | GWAS2 |  | -0.00351761 | unavailable |
| gcorr_chr22 | cov | GWAS1 | GWAS2 | 0.00177085 | unavailable |
| gcorr_chr22 | rg | GWAS1 | GWAS2 | unavailable | unavailable |
| gcorr_genome | h2 | GWAS1 |  | 0.680078 | unavailable |
| gcorr_genome | h2 | GWAS2 |  | -0.187175 | unavailable |
| gcorr_genome | cov | GWAS1 | GWAS2 | 0.106081 | unavailable |
| gcorr_genome | rg | GWAS1 | GWAS2 | unavailable | unavailable |

Chr22 uses K <= 50 and defines a chromosome-only modeled universe. Genome uses K <= 20 across all 1,702 regions. These are different spectral/model specifications.

Actual retained ranks and correction denominators:

| trait | sample_size | retained_rank | correction_denominator | analysis |
| --- | --- | --- | --- | --- |
| GWAS1 | 49071.3 |    1200 | 0.975546 | gcorr_chr22 |
| GWAS2 | 40743.2 |    1200 | 0.970547 | gcorr_chr22 |
| GWAS1 | 49071.3 |   33971 | 0.307721 | gcorr_genome |
| GWAS2 | 40743.2 |   33971 | 0.166216 | gcorr_genome |

Reference LD spectral fraction retained by the genome cap (not a measured fraction of genetic signal):

| summary | retained_reference_eigenvalue_fraction |
| --- | --- |
| minimum | 0.24941 |
| median | 0.558405 |
| maximum | 0.991234 |

Estimate/uncertainty availability (raw points are preserved):

| analysis | rows | negative_regional_h2 | unavailable_points | unavailable_standard_errors | rg_outside_unit_interval |
| --- | --- | --- | --- | --- | --- |
| gcorr_chr22 |     100 |      26 |      19 |     100 |       4 |
| gcorr_genome |    6812 |    1353 |    1353 |    6812 |      91 |

Native reasons for unavailable uncertainty: Joint HESS uncertainty requires valid regional genetic and total residual parameters

All regional h2/covariance contributions sum to their corresponding joint totals. Correlations are ratios and do not sum. A finite point estimate is not evidence of calibrated uncertainty.

**The corrected CAD genome total is negative, so genome rg is undefined. All HESS/rho-HESS standard errors are unavailable under the native conditional checks. These results do not support interpretable CAD/T2D genetic-correlation inference under this configuration.**

The [focused HESS investigation](../hess-investigation/README.md) reproduces the
compact projections independently, explains the projected-energy shortfall and
the uncertainty gates, and records genomic-control/rank/reference diagnostics.
It preserves these numerical results and does not apply an inferred correction.

## Timing and memory

| analysis | threads | preparation_seconds | fit_seconds | analysis_seconds | serialization_seconds | process_wall_seconds | peak_working_set_GiB | sampled_peak_private_GiB |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gcorr_gnova |       4 | unavailable |   291.1 |  311.62 |    3.35 | 450.672 | 6.96782 | 16.8697 |
| gcorr_sumher |       4 | unavailable |     379 |  399.36 |    5.58 |  511.26 | 5.78675 | 14.8485 |
| gcorr_chr22 |       4 |  213.37 |    0.14 |  215.41 |    0.09 | 319.284 | 5.20144 | 7.11275 |
| gcorr_genome |       4 | 1048.58 |   58.16 | 1135.14 |   18.06 |  1275.5 | 6.43116 | 10.9847 |

![Corrected preparation time](gcorr-current-preparation-time.png)

![Corrected regional covariance](gcorr-current-regional-covariance.png)

Input/LDSC run: manuscript-20261007-003431
Score/chr22 run: manuscript-20261007-064651
Genome run: manuscript-20261007-070903

The subsequent [robust HESS adapter and saved-moment checks](../hess-robustness/README.md)
add automatic rank budgeting, explicit genomic-control provenance and inference
blocker diagnostics. The measured estimates and timings above are preserved;
the subsequent [signed-moment and automatic-rank confirmation](../hess-moment-uncertainty/README.md) records the new run separately.
