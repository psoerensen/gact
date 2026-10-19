# Standalone counterparts after the original qgg website timings

Authorized by the user on 4 October 2026, after completion of the running qgg
gene examples. Do not launch while those timings are active. Preserve the old
simulated benchmarks; they are not matched comparisons for these workloads.

Use the same 503 EUR individuals, downloaded genotypes, checked CAD/T2D summaries
and gene/pathway memberships, at one thread. Gene LD is computed on the fly.
The primary clock includes required genotype/LD preparation, spectra and model
evaluation. No persistent gene-LD or eigenspectrum cache is introduced.

The current public glma/gmap interfaces require CSR LD. A benchmark adapter
copies the exact selected packed BED records into a temporary gene-sized PLINK
triplet; it does not simulate, recode or change the individuals. Public gbase
preparation validates that triplet and computes full within-gene LD with r2=0,
both windows disabled, and explicit allow_full_ld=TRUE. Public glma/gmap then
fit the gene. Only the fixed temporary triplet and five owned CSR files are
removed after a successful fit. Failed inputs remain diagnostic evidence.
Copying, public preparation/validation, and cleanup are inside the primary
operation clock. Component times only explain that combined cost.

This adapter's conversion/file overhead is part of the measured standalone
workflow; it is not an intrinsic native-solver timing. Temporary CSR persistence
uses float32 correlations whereas qgg calculates double-precision cor(B).
Compare quadratic statistics and report p-value/tail differences. Do not claim
bit-identical numerical results or an in-memory BED API that is not available.

## Workloads and controls

* VEGAS: every gene retained by the qgg multi-trait example, both studies,
  identical QC marker counts. qgg p floor 1e-10 is retained before forming Z.
  Public glma saddlepoint tails; relative PSD tolerance 1e-7; maximum_markers
  explicitly raised for the largest supplied set instead of dropping it.
  Gene spectra are computed once for both traits within each call, not reused
  between timed runs. Save per-gene component clocks to examine large outliers.
* Joint Bayesian pathways: identical cached 17,955-gene matrix and 99 KEGG sets,
  shared available-gene universe, BayesC, 200 burn-in and 1,000 sampling sweeps,
  one seed/chain. The original null/three-active pattern initialization is
  (.99,.01/3,.01/3,.01/3). New gsea learns effect covariance and pattern weights;
  residual trait covariance is explicitly fixed at .5 I, matching the original
  initialization. Original qgg learns its residual covariance; this policy is
  different, so this is a corresponding workload, not exact model equivalence.
* Fine-mapping: four original genes, their stricter QC marker sets, BayesR,
  100 burn-in and 1,000 sampling sweeps, seed 10, one chain. Use the same
  .995/.99/.95/.9 schedule with selection on convergence disabled, corresponding
  to the original unchecked first-threshold workflow. No convergence claim.
  Direct score/information inputs use the qgg dense working-summary ww/wy
  formulas, with recorded scale/initial variance and (.99,.006,.003,.001) weights.
  The original qgg eigen transform differs; this is a corresponding workflow,
  not a numerical reproduction of that transform. New gmap fixes residual
  variance because its learning contract requires in-sample LD/full rank;
  effect variance and weights are learned. Do not call gmap_stat's OLS helper
  on these binary meta-analysis summaries or assert in-sample LD.

Whole-genome fitting/scoring remain blocked by the original LD resource gate.
Do not substitute cached/simulated effects and label them as completed website
fits. A further standalone whole-genome workload needs its own resource and
control audit. No package interfaces or numerical implementations are changed
by this benchmark adapter.
