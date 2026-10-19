# Published qgg/gact examples: EUR-only execution

User authorized direct website-example timings on 4 October 2026. This replaces
the proposed regional/shared-set design in qgg-website-protocol.md. Existing
measurements, the database and failed attempts are preserved. Use installed
qgg 1.1.6 and gact 0.0.0.1, one OpenMP/BLAS thread, fixed R seed 20261004.

Follow the published examples, changing local paths and documenting required
API compatibility adjustments. Do not substitute standalone implementations,
simulation-derived weights, regional gene samples or redesigned controls.

1. Process_1000G: downloadDB(what="1000G"), with the local descriptor's source
   list restricted to documented EUR DOI 10.5281/zenodo.14141438. The original
   descriptor/database is unchanged. Download to a fresh stage directory.
   writeBED filters by GAlist$rsids; gprep uses the published explicit BED/BIM/FAM
   arguments and study label "1000G EUR".
2. Compute_sparseLD_1000G: gprep(task="sparseld", msize=1000,
   rsids=unlist(Glist$rsids), overwrite=FALSE), one thread and a fresh LD path.
   Preserve LD instead of the tutorial's deletion so model/scoring examples
   can reuse it and sizes can be audited.
3. Gene_analysis_vegas: both CAD/T2D studies, all installed Genesplus sets.
   qgg checkStat only accepts data frames, whereas gact's multi-study accessor
   returns a list. Check studies separately and pass aligned p-values through
   qgg's public multi-trait p interface. QC and VEGAS defaults otherwise intact.
4. Gene_set_analysis_bayesian_magma: getVEGAS on GWAS1/GWAS2 and
   getFeatureSets(feature="KEGG", minsets=50); magma BayesC pi=.01,
   nit=1000, nburn=200. This uses the saved database gene evidence, as published.
5. Finemapping_bayesian_linear_regression_real_data, gene-region example:
   GWAS1, published QC, four named genes ENSG00000099204/ENSG00000118972/
   ENSG00000145996/ENSG00000140718; BayesR mcmc-eigen with thresholds
   .995/.99/.95/.9 and installed defaults for unspecified arguments. Capture
   Geweke diagnostics. checkConvergence defaults FALSE, so returned success
   does not establish convergence or trigger eigen retries.
6. Polygenic_scoring_bayesian_linear_regression: use the EUR Glist as the
   ancestry-specific replacement for the generic Glist filename. Replace the
   obsolete getMarkerStatDB name with its installed getMarkerStat equivalent.
   Fit BayesC nit=1000, then total score. Fit BayesR nit=5000/nburn=1000,
   then total and KEGG pathway scores; reuse that identical BayesR fit rather
   than time its repeated identical invocation twice. All EUR reference people
   are score recipients; fitting and scoring have separate clocks.

Plotting/printing examples are excluded from analysis clocks. verbose=FALSE
avoids rendering diagnostic plots; numerical diagnostics are saved. Primary
analysis wall time excludes loading, previously completed QC and serialization;
the VEGAS operation includes constructing its aligned two-trait p matrix. Process wall
time and peak R working set cover the whole stage. Fresh processes per stage,
100-ms private-memory sampling, exact output/resource sizes and failure logs.

Sources: https://psoerensen.github.io/gact/Document/Process_1000G.html,
Compute_sparseLD_1000G.html, Gene_analysis_vegas.html,
Gene_set_analysis_bayesian_magma.html,
Finemapping_bayesian_linear_regression_real_data.html,
Polygenic_scoring_bayesian_linear_regression.html.

## Resource gate observed during execution

The published filter retained 9,985,562 markers in 503 EUR people. The
msize=1000 LD representation requires 79,924,438,248 bytes (74.44 GiB),
whereas C: had 47,847,862,272 bytes free at the stage check. Computation was
blocked before writing LD; the failed-stage wall time is not an LD benchmark.
The successful Glist and all prior results remain available. Gene-level VEGAS
and gmap calculate LD directly from BED, and joint MAGMA uses database gene
statistics, so these independent examples can proceed using the prepared
Glist. This does not qualify whole-genome fitting or scoring as completed.

The EUR download operation completed with validated files, but its enclosing
process exited nonzero after a worker-footer parse error caused by editing the
worker during that run. Preserve and label that diagnostic run; the download
operation timing is available, but it is not a clean full-process timing.

## Authorized follow-up

After the qgg analysis timings complete, the user requested the corresponding
new standalone R packages in gsuite. Reuse the same EUR reference, GWAS data,
QC-retained marker identities, gene/pathway memberships and one-thread policy.
Match iteration/burn-in settings where interfaces support them, and document
method/LD/tail/control differences explicitly. Do not reuse the earlier
simulated 20K/50K figures as a matched comparator. Preparation/LD resource costs
must be separate from model clocks. Audit the public interfaces and resource
requirements before launching: standalone VEGAS consumes a prepared sparse-LD
Glist, whereas original qgg VEGAS calculates full within-gene correlations from
BED, so equal membership alone does not establish equal LD inputs. Preserve
the original qgg run, including any failed/blocked stages.
The user clarified that gene LD is normally calculated on the fly. Primary
gene-analysis clocks must therefore include computing within-gene LD for both
implementations. Do not introduce a persistent gene-LD cache or present
evaluation of cached spectra as the primary matched workload. If the current
standalone public interface requires temporary CSR resources, create them
within that timed workflow and remove only their explicitly owned files after
use, preserving timing/provenance/output evidence. Keep secondary component
clocks solely to explain the combined wall time.
