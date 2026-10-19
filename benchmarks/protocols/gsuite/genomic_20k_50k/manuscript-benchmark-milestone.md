# Minimal gact manuscript performance benchmark

## Dense regional fitting and familiar interfaces, 6 October

The [new comparison and verified archive](manuscript-results/manuscript-simulated-finemapping-dense-20261006/README.md)
preserve the frozen tutorial inputs and every earlier report. Three serial
one-thread repeats gave regional medians of 6.14 s for qgg and 5.48 s for the
standalone dense route. A fresh standalone eigen profile took 116.00 s on the
same current runtime: the dense route was about 21 times faster for this
four-region workload. Genome-wide medians were 9.96 s and 2.36 s, respectively,
using already prepared LD. Median process peaks were 0.84/0.56 GiB regionally
and 2.72/0.81 GiB genome-wide. Ranges and phase timings are in the report.

The qgg-style gmap/gbayes adapters and gbase marker/LD helpers are implemented.
Standalone `gmap::getCredibleSets()` now replaces the missing helper reported
below. Its groups and regional fit groups matched qgg exactly on the same
PIPs; full regional LD agreed within 3.6e-15. All regional PIPs matched the fresh
standalone eigen profile, and all four causal PIPs remained one in both
engines. The distinct likelihood/residual policies, effect units and
single-chain limits remain documented. Focused R/native tests passed. The gact
source example and local public documentation snapshot are updated; nothing
has been committed, pushed or publicly published by this milestone.

## Earlier completed comparisons, 5 October

The [complete regional VEGAS workflow](manuscript-results/manuscript-vegas-complete-20261005/README.md)
now has a successful standalone timing: 17,862 genes and two traits in 15.00
minutes, including full regional BED spectra and native evaluation. Its different
tail approximations and the original-qgg tail comparison are recorded separately.
This supersedes the failed initial standalone VEGAS attempt mentioned below.

The [published simulated-human fine-mapping comparison](manuscript-results/manuscript-simulated-finemapping-current-20261005/README.md)
uses frozen tutorial inputs, one-thread clocks and refreshed 0.1.0 standalone
packages. Four regional BayesC fits took 15.87 s in qgg and 122.66 s in gbase/gmap,
including regional LD in both. Genome BayesC fitting took 13.51 s and 7.78 s respectively,
with LD already prepared. Fresh sparse-LD construction took 190.53 s and 198.31 s.
All four causal markers have PIP 1 in both fits. The report retains different
likelihood/residual policies, single-chain diagnostic limits, the local BIM-map
correction and the missing standalone summary-fit posterior-set interface. The
common genome posterior-set comparison uses qgg postprocessing on both fits.
The earlier cached-package report remains intact. Current seeded standalone PIPs
and posterior means exactly match that earlier run; these single execution clocks
do not establish a package-refresh performance effect.

The consolidated table linked below records the earlier 30 executions; these
new comparisons have separate reports and never overwrite the original evidence.

Current manuscript evidence measures the published gact website examples in
installed qgg/gact. Database installation, EUR reference preparation, CAD/T2D
input QC, full-gene VEGAS, joint KEGG Bayesian MAGMA and four-gene fine-mapping
have completed. The gene-analysis stages requested four threads; preparation
and the additional website examples requested one. These are single descriptive
executions, not a uniform one-thread benchmark or a scaling study.

PoPS, the qgg Bayesian MAGMA examples, all four LDSC examples, and simulated
BayesC fine-mapping with credible-set reporting also completed on 5 October.
The current [Word overview](manuscript-results/20261005-structured/gact-benchmark-current-overview.docx)
records their scope and limitations. The [website workflow report](manuscript-results/20261004-qgg-website/README.md)
is an earlier execution snapshot; its statements about active stages are
historical. The [protocol](qgg-examples-protocol.md) records the original controls.
New genome-wide EUR LD preparation remains storage-gated; individual total and
pathway scoring and a matched standalone MAGMA comparison are incomplete.

The [current stage measurements](manuscript-results/20261005-structured/completed-stage-measurements-current.csv)
consolidate 30 successful executions with decimal-point numeric fields, relative
source paths and hashes of the original timing and process-measurement files.
Operation clocks, whole-process clocks, working-set memory, sampled private
memory and saved RDS sizes are distinct columns. Earlier tables remain intact.
This inventory records completed execution, not universal numerical or MCMC
qualification; failed and blocked stages remain in their original run records.

Corresponding standalone gsea joint-pathway and gmap four-gene runs have
completed, with differences in residual covariance/variance treatment and model
transformations. They are not exact reproductions or evidence of a matched
speedup. The initial standalone VEGAS attempt failed with an LD-resource reader
error; its successful replacement is linked above. One- and four-thread gene-preparation component
measurements are saved separately, including the flagged extraction outliers.
The [comparison protocol](qgg-package-comparison-protocol.md) requires on-the-fly
LD preparation inside primary clocks. Further integration into qgg is future
work; current standalone measurements must not replace original qgg manuscript
claims. Preserve all earlier evidence below.

## Earlier simulation study

The real-1000G-derived one-thread benchmark is complete. See the
[measured report](manuscript-results/20261004-real-1000g/README.md),
[exact controls](real-1000g-protocol.md) and [input audit](real-ld-input-audit.md).

The new standalone R packages are the primary implementations intended for
qgg/gact integration. Original gact/qgg models remain historical references;
no model equivalence is asserted. Simulated A/B traits illustrate the workflow
and must not be described as measured manuscript CAD/T2D analyses.

1. Input audit complete: cached official 1000G Phase 3 chr22, GRCh37 map and
   population panel verified; 503 EUR donors selected.
2. Real-reference preparation complete: public gsim generated 20,000 people
   (19,998 founders and two offspring), with 48,088 imported markers and 5,360
   retained at MAF >= .01. The cohort is preserved for reuse.
3. Median-size timing sets complete: 1,000 contiguous windows of 274 markers
   match the raw installed Genesplus median. They are timing sets, not mapped
   biological genes or a post-QC gene-size audit.
4. One-thread measurements complete: VEGAS 44.78 s; joint Bayesian pathways
   0.10 s; BayesR mcmc-eigen mapping 195.65 s; total scoring 1.03 s; three pathway
   scores 2.41 s. Memory, serialization, output/resource sizes and provenance
   are recorded. Dimensions and finite-output checks passed.
5. Database installation audit complete and reused: download/install 544.20 s,
   descriptor saving 12.48 s, peak R working set 2.677 GiB, installed files
   3,261,689,414 bytes. All seven component categories succeeded.

Mapping accepted only 3 of 40 region-trait fits; 37 are provisional. Its bounded
runtime includes convergence-triggered eigen retries and is not evidence of
inferential convergence. Float32 LD uses recorded rounding tolerances of 1e-7.
Scoring uses deterministic sparse simulation-effect weights, not posterior
weights: the auxiliary genome-wide BayesC fit overflowed on window-truncated
LD and was excluded. All diagnostic runs are preserved.

The previous artificial-LD one/four-thread measurements remain in their
[original report](manuscript-results/20261004/README.md), superseded for analysis
reporting. Their thread effects do not establish speedups on this real-LD panel.
No additional four-thread run was requested after selecting one-thread reporting.

A matched official standalone MAGMA gene-analysis comparison is specified in
the report for Reviewer 2 but has not been executed. Installed gbase source
attribution remains limited by four DLL hashes that differ from the older build
manifest; exact measured binary hashes are recorded. No commits or pushes made.
