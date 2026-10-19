# BLR fitting from precomputed cross-products

One execution on 19 October 2026 using the existing gsuite native gbayes library and the existing installed-database gene scores. No gene-association computation or genotype analysis was run. The 22-fold validation timings remain separate.

| Stage | Wall seconds | Boundary |
| --- | ---: | --- |
| Load cached gene scores and annotation fixture | 0.11 | No gene testing |
| Construct sparse gene-feature matrix and standardize response | 0.46 | 17,955 genes; 1,324 Reactome features; 39,633 nonzero memberships |
| Prepare centered/scaled X'X, X'y and y'y | 0.20 | Sparse Matrix cross-products; outside sampling clock |
| Load serialized cross-products into native worker | 0.0153 | Outside sampling clock |
| Validate/bind summaries and set model controls | 0.0037 | Existing dense_cross_products / FixedGaussianStatistics API |
| BayesR native fit from prepared summaries | 14.2124 | One chain; 1,000 burn-in + 5,000 sampling sweeps; thinning 1; one thread; native validation/sampling/posterior assembly |
| Write posterior coefficients and parameter traces | 0.0424 | Outside sampling clock |

Compilation, process startup and the separate qgg input-validation check are excluded. These are single measurements, not repeated-run averages. Native-process peak working set is 22,618,112 bytes (21.57 MiB); it excludes the separate R preparation process and is not an end-to-end R workflow peak.

## Matching qgg

The installed qgg 1.1.6 R bodies for pops(), computeStat() and blr() were checked against the local source. The sparse summaries match installed qgg::computeStat() to floating-point precision: maximum XX difference 1.10e-11 and maximum Xy difference 2.28e-13. Gene evidence is standardized once, as in qgg PoPS; feature cross-products are centered and scaled without making the gene-by-feature matrix dense.

The actual qgg BayesR defaults are used: starting weights (0.95, 0.02, 0.02, 0.01), variance multipliers (0, 0.01, 0.1, 1), initial ve=0.5 and vb=0.5/(1324*0.001), inverse-chi-square df=4 with scales 0.25 for residual variance and 0.5*vb for effect variance, and Dirichlet(1,1,1,1) weight learning. Residual likelihood dimension is n=17,955 to match qgg's update, not n-1. All three parameter-learning flags are enabled. The seed is 32132745, obtained from the R seed draw qgg::blr() would use after set.seed(10).

For a dense XX, every qgg LD neighborhood contains all features: adjustE contributes zero. Consequently updateG changes a reported genetic-variance quantity but not the dense coefficient/residual updates. No nonzero reference-LD adjustment is substituted.

This matches the model scale, priors and sampling controls; it does not reproduce a bitwise qgg trajectory. Native marker visit order/RNG consumption and retained outputs differ. qgg retains all full coefficient draws; this native driver retains posterior summaries and parameter traces. No new qgg sampler execution or runtime comparison was performed.

## What exists and what is missing

The public native APIs gbayes::dense_cross_products(), FixedGaussianStatistics and fit_bayesr() already implement sufficient-statistic BLR. This measurement uses a thin benchmark executable linked to the existing installed native library. No new sampler was written and no native library was rebuilt.

The public R gbayes() wrapper currently expects a genotype/reference descriptor; it has no direct XX/Xy/yy/n entry point. The existing gsea R route also rebuilds its general covariance/design preparation. Therefore this is a native-engine benchmark, not a timing of a currently exposed R summary-BLR function.

The earlier 162.87-second gsea R call remains archived with its own two-chain settings and preparation boundary. It is not a sampler-only timing, has different priors/controls and must not be treated as a controlled speedup comparison with this one-chain native run.

## Proposed implementation

1. Add a serializable, checked prepared-BLR R object and a direct gbayes::blr() entry point accepting XX, Xy, yy, n and feature labels. Check dimensions/order, finite values, symmetry, normalization and prior/control settings; pass through to the existing native APIs. Do not fabricate SNP coordinates or LD files for gene features.
2. Add sparse preparation in the gene-feature workflow, preserving centering/scaling and gene/feature ordering. Cache the summaries for repeated priors/chains and return explicit preparation, sampling and output clocks. Keep the optional chromosome-held-out validation workflow separately labeled.
3. Audit the native exact-variance BayesR path: current source calls residual->rebuild(b) every sweep, which performs a dense matrix-vector multiplication. The sampler already maintains r = X'y - XX*b incrementally and can compute SSE from that residual. First profile this path and add safe zero-column fast paths: the dense matrix multiply currently visits columns even when the corresponding effect is zero. Preserve the existing exact refresh and numerical checks for that optimization. Only if needed, evaluate a separate explicit periodic/adaptive rebuild policy with drift monitoring and a final rebuild, documenting any changed numerical-refresh policy. This is a candidate performance cost, not a measured phase attribution or a promised speedup.
4. Verify sparse summaries against qgg, scalar conditional/variance-update formulas, fixed-seed prepared-path consistency and residual drift on small dense examples. Re-time the same frozen full workload only after those checks pass. Retain the existing 14.21-second result and identify the exact native DLL/profile hashes for any new run.

The implementation plan does not claim that native outputs are already qgg-identical or that Bayesian chains have converged. Sources, settings, input/output hashes, package session and the preliminary input-format diagnostic are retained alongside this report; raw matrices, binaries and fitted outputs remain Git ignored.
