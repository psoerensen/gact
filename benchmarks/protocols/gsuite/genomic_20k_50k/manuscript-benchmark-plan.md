# Proposed minimal gact manuscript performance benchmark

Historical audit and initial proposal. The current scope is defined in
[the benchmark milestone](manuscript-benchmark-milestone.md): standalone R
packages are primary and analyses use one and four threads. The original-qgg
priority and one-thread proposal below are superseded.

Status: superseded initial proposal; the milestone and measured report contain
the completed standalone analysis and database-installation results.
Audit date: 4 October 2026.

## Existing evidence and reuse

Preserve `results/`, the existing drivers, and
`build/benchmarks/genomic-20k-50k/` (including its qgg library and frozen data).
The existing 20,000-individual / 50,000-marker simulated benchmark supplies
the stage logging, isolated-library and provenance conventions, not manuscript
data. Its timings are one descriptive execution. Its qgg mapping is BayesC
`mcmc-eigen`; its standalone mapping is multi-effect. Neither is the requested
BayesR workload. Existing physical RAM provenance is not peak process memory.
Do not rerun the existing drivers into their current result filenames.

The local gact tutorials reference
`C:/Users/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds`; that file was not present.
No matching manuscript/database file was found in the targeted local gact,
qgdata, human, aipublish, examples or Documents/Data searches. This is a bounded
search, not proof that the manuscript data are absent from this computer.

## Workloads to agree

The user selected the public download/install example for database timing:
`gact::gact(version="hsa.0.0.1", dbdir=<fresh benchmark directory>, task="download")`,
followed by saving GAlist. Measure this as **database download and installation**,
including network transfer, extraction and metadata loading, rather than database
construction from raw source data. Package installation remains outside timing.
Record total time and, where observable without changing the workflow, component
times and transferred bytes. Do not infer local computation performance from
network-dependent total time.

The local installer requests annotation, marker, gsea, gstat, gbayes, atc and dgi
components; optional raw reference downloads are not enabled by this call.
Inspect the installed manifest to establish whether usable LD/BED resources
are included before planning analysis stages. The installer catches some
download errors and returns a logged GAlist, so successful function return
alone is insufficient: validate download status and required files.

The user selected the existing simulated 20K/50K workloads. Reuse the frozen
PLINK files, reference (18,000 individuals), summaries (6,000 individuals per
trait) and held-out score recipients (2,000 individuals). Keep labels A/B and
describe these as representative synthetic timings, not measured CAD/T2D data.

### Minimal synthetic protocol proposed for final agreement

| Timed stage | Frozen inputs and controls |
| --- | --- |
| VEGAS | qgg implementation; traits A/B; all 1,000 existing artificial genes; existing qgg sparse LD |
| Joint BLR-MAGMA | qgg `magma` on two-column A/B VEGAS Z scores; three existing overlapping pathway sets; BayesC, pi=.01, nit=1000, nburn=200; set.seed(20260932) |
| BayesR fine-mapping | qgg `gmap`, A/B separately, first 20 existing 100-marker regions per trait; mcmc-eigen; thresholds .995/.99/.95; nit=1200, nburn=500, nthin=1, nrun=1, seed=10; h2=.3; updateG/updateE/updatePi=TRUE; other current qgg defaults captured explicitly in run manifest |
| Total scoring | qgg `gscore`; frozen first-chain qgg-current BayesC posterior weights; 2,000 held-out individuals; scaleMarker=TRUE, scaleGRS=FALSE, msize=128, ncores=1 |
| Pathway scoring | Same frozen weights/recipients and controls; marker unions of the three existing pathway sets, deduplicated within each pathway and permitting overlap across pathways |

Use one thread for the original qgg analysis table and one descriptive execution.
Do not regenerate simulation, rerun GWAS/LDSC/PoPS, or fit another genome-wide
scoring model. Loading the frozen weights is outside score application timing;
the result must disclose that the weights are BayesC, not BayesR mapping output.
The three pathways and 1,000 genes are artificial, not KEGG manuscript annotations.
This updates the existing scorer's 20 gene-set workload to actual synthetic
pathway membership; report that scope change rather than reuse its old timings.

Standalone results remain separate existing evidence unless an additional
matched standalone run is explicitly selected. The optional official MAGMA
comparison would use trait A and the same artificial reference, gene definitions
and three pathways; it remains pending agreement and executable availability.

| Stage | Proposed original gact/qgg workload | Required decisions |
| --- | --- | --- |
| Database preparation | Prepare a fresh output database from frozen local CAD/T2D summaries, EUR reference and required annotations | Database/version/build, study identities, input files, QC, annotation scope; separate local preparation from downloads and copying an existing database |
| VEGAS | CAD and T2D over the same frozen Genesplus marker sets and EUR LD | Gene universe, 40 kb upstream / 10 kb downstream mapping, QC, shared preparation and per-trait evaluation boundaries |
| Joint BLR-MAGMA | One two-column CAD/T2D fit using qgg `magma(method="bayesC", pi=.01, nit=1000, nburn=200)` and KEGG `minsets=50` | Confirm actual manuscript controls; tutorial uses GWAS1/GWAS2, which must be verified against database metadata |
| Fine-mapping | qgg `gmap(method="bayesR", algorithm="mcmc-eigen")` | Choose the four tutorial genes or chromosome-17 LD regions; fix nit, nburn, nrun, seeds, eigen thresholds, hyperparameter learning and credible-set controls explicitly |
| Scoring | Separate total and KEGG pathway `qgg::gscore` application using one frozen posterior weight table and recipient BED dataset | Recipient count/ancestry, markers, pathway universe, scaling and whether weights come from the measured mapping or an existing genome-wide fit |

The tutorial alternatives are not one frozen manuscript protocol. The fine-mapping
four-gene example uses thresholds .995/.99/.95/.9 with implicit sampler defaults;
the chromosome-17 example uses 5000 iterations, 1000 burn-in and
.99/.97/.95/.9. Scoring tutorial weights come from a separate genome-wide
BayesR fit. These differences must be resolved before timing.

## Measurement contract

Use one descriptive execution per agreed workload and one fixed thread budget;
no scaling campaign or repeated timing study. Launch each stage in a fresh R
process with explicit isolated libraries. Record parent wall time and analysis
wall time separately so startup/loading, computation and serialization boundaries
are clear. Reused inputs must be labelled; resumed results are not new replicates.

On Windows record the process peak working set (bytes), with peak private bytes
as a distinct metric if available. Report the measurement API and whether child
processes exist; a parent-only high-water mark is not aggregate process-tree
memory. Do not substitute R allocation counts or machine RAM for process memory.

Record input bytes, newly written output bytes, final database bytes and per-stage
database growth using file lengths; identify shared files so totals do not double
count them. State whether compressed files and logs are included. Preserve input
hashes, retained dimensions, warnings/failures and lightweight output checks.

Record CPU model, physical/logical cores, RAM, OS, R/sessionInfo, installed package
versions and library paths, source revisions and dirty status, compiler/native
library provenance where available, BLAS/OpenMP settings, thread limits, seeds,
full arguments, timestamps and storage location/type. Distinguish installed
binary provenance from the current source checkout revision.

Store new large intermediates in a separate run directory under the existing
`build/benchmarks/work/` tree, and compact review tables in a new run directory under
this benchmark. Refuse existing output destinations. Never modify the source
database or existing results. No measured values are available yet.

## Standalone packages and Reviewer 2

The primary manuscript table should measure the original gact/qgg workflow.
Any standalone measurements must be separately labelled: gbase preparation,
glma VEGAS (re-exported by gsea), gsea Bayesian pathways, gmap/gbayes mapping
and gscore application. Shared purpose does not establish matching priors,
likelihoods, tail evaluation, defaults or inference. No replacement model
should silently enter the original table.

Consider one official MAGMA executable comparison for a single agreed trait:
same frozen SNP intersection, EUR reference individuals, genome build,
gene boundaries, pathway universe and covariates. Time annotation, SNP-wise
mean gene analysis and competitive pathway analysis separately, plus their
total. Pin executable version/hash and exact commands. This is a workflow
cost comparison; official competitive MAGMA is not the joint Bayesian
CAD/T2D BLR-MAGMA model. Confirm Reviewer 2's requested comparison and whether
an official executable and compatible local inputs are already available.

## Gate before execution

Obtain manuscript/reviewer text and data location; inspect study metadata and
file manifest read-only; agree the above exact workload dimensions and controls.
Only then prepare concrete stage drivers and execute the agreed benchmark.
