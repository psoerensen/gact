# Standalone workflow performance measurements

**Input limitation identified during review:** the frozen benchmark driver
generates genotypes using binomial draws and within-block copying. It uses gsim
to generate phenotypes on that BED, not gsim's real-reference founder workflow.
These results therefore measure artificial LD and are provisional for the
manuscript's requested real-LD benchmark. A size-matched real-reference workload
has now been located, authorized and started under the
[real-1000G protocol](../../real-1000g-protocol.md). No artificial chromosome remapping is accepted.
The [completed real-1000G report](../20261004-real-1000g/README.md) supersedes
the analysis table below; the database evidence here remains applicable.
The attempted extension in `manuscript-20261004-152758` failed during preparation
and produced no valid timing; its diagnostic files are preserved and its driver
extension has been withdrawn.

Measured on 4 October 2026 using the existing frozen 20K-individual / 50K-marker
simulation. These are representative synthetic timings, not measured CAD/T2D
analyses or a method-ranking study. One execution per thread setting is reported.
The original benchmark results were preserved. At the user's request, the
one-thread measurements are the primary manuscript results. Four-thread
measurements are retained as supporting evidence, not as a second primary table.

## Primary manuscript table: one thread

| Workload | Wall time (seconds) | Peak R working set | Serialized output / database size |
| --- | ---: | ---: | ---: |
| Database download/install | 544.20 | 2.677 GiB | 3.262 GB database files |
| VEGAS, A/B and 1,000 genes | 6.80 | 203.3 MiB | 351.2 KiB output |
| Joint Bayesian pathways, A/B and 3 pathways | 0.11 | 147.5 MiB | 12.3 KiB output |
| BayesR eigen mapping, 20 regions per trait | 24.00 | 189.2 MiB | 1037.7 KiB output |
| Total scoring, 2,000 recipients | 8.03 | 192.3 MiB | 748.9 KiB output |
| Three pathway scores, 2,000 recipients | 6.30 | 192.8 MiB | 840.8 KiB output |

Analysis wall times exclude loading/preparation and serialization, as detailed
below. Database descriptor saving adds 12.48 seconds and 52,318,792 bytes;
download/install plus descriptor saving totals 556.68 seconds. Peak memory covers
the entire R process lifetime. GB is decimal; GiB/MiB/KiB use powers of 1024.
These synthetic workloads are small and fast. More traits or larger analyses may
offer more parallel work, but their scaling has not been measured here. The
supporting runs show a clear benefit for scoring at this scale.

## Supporting thread comparison

| Stage | Analysis seconds, 1 thread | Analysis seconds, 4 threads | Peak R working set MiB, 1 / 4 | Serialized output KiB, 1 / 4 |
| --- | ---: | ---: | ---: | ---: |
| VEGAS, A/B and 1,000 genes | 6.80 | 8.97 | 203.3 / 206.5 | 351.2 / 351.2 |
| Joint Bayesian pathways, A/B and 3 pathways | 0.11 | 0.17 | 147.5 / 149.3 | 12.3 / 12.4 |
| BayesR eigen mapping, 20 regions per trait | 24.00 | 24.08 | 189.2 / 191.7 | 1037.7 / 1037.7 |
| Total scoring, 2,000 recipients | 8.03 | 3.23 | 192.3 / 192.9 | 748.9 / 748.9 |
| Three pathway scores, 2,000 recipients | 6.30 | 2.66 | 192.8 / 192.7 | 840.8 / 840.8 |

The mapping columns vary the process OpenMP limit; gmap exposes no thread
control and loops over regions/traits serially. They are not evidence of
four-thread mapping. VEGAS reports four preparation threads; evaluation remains
serial. The joint sampler reports four threads across four identical chains.
Scoring exposes a native thread setting but does not report an observed runtime
thread count. BLAS is limited to one thread for both settings.

All substantive output comparisons passed at relative tolerance 1e-10:
VEGAS statistics, pathway PIPs/results, mapping estimates/parameters/sets and
score matrices. Scoring outputs are finite and contain 2,000 recipients; joint
pathway output contains six pathway/trait rows. This is an execution check,
not statistical calibration or convergence qualification.

## Exact scope and model distinctions

### Gene-set size audit

The 1,000 synthetic genes partition the 50,000 markers without overlap. Set
sizes range from 20 to 80 markers, with mean/median 50 and quartiles 35/65.
They are a deliberately bounded computational fixture, not a realistic
genome-wide distribution of human gene-marker sets.

The installed gact `Genesplus` annotation contains 18,279 sets before
GWAS/reference intersection or QC: median 274 markers, quartiles 183/438,
90th percentile 767, 99th percentile 2,465.32 and maximum 20,376. Only 2.22%
fall between 20 and 80 markers. These raw annotation counts are not the
post-alignment VEGAS dimensions; matching a particular GWAS and LD reference
will change them. `gene-set-size-audit.csv` records this read-only audit.
Consequently the reported VEGAS timing does not establish genome-wide human
performance or the cost of large gene spectra. The simulated sets also omit
overlapping gene windows and derive their boundaries from artificial block splits.

VEGAS includes gene-spectrum preparation and evaluation of both traits, using
the saddlepoint tail approximation and saved sparse LD. Joint analysis reuses
the frozen gene evidence. It uses gsea/gbayes BayesC with independent gene
sampling covariance, identity trait sampling and effect covariance, scaled
responses/features, and four chains with seeds 11/29/47/71. Each chain has
200 burn-in and 800 sampling sweeps. Pattern weights are .985 null, .005 A-only,
.005 B-only and .005 both, giving marginal inclusion probability .01 per trait.
Covariance and pattern probabilities remain fixed. This is not a claim of
identity to original qgg joint BLR-MAGMA.

Mapping uses gmap/gbayes BayesR with 500 burn-in and 700 sampling sweeps per
chain, the same four seeds, residual variance 1, effect variance .04, mixture
weights .7/.1/.1/.1 and variance multipliers 0/.01/.1/1. Variances/weights are
fixed. Eigen retention uses .995/.99/.95 with convergence-based selection;
inspect the saved diagnostics for attempts and retained ranks. Regions contain
100 markers each. This replaces the historical multi-effect workload with an
explicit standalone BayesR model, not a qgg `mcmc-eigen` equivalence claim.

Scoring uses the frozen standalone BayesC posterior weights, standardized marker
coding and unscaled scores. Pathways use deduplicated unions of the existing
synthetic gene memberships, preserving overlap between pathways. No new
weight fitting is timed. Preparation/loading is outside analysis time.

## Measurement boundaries and reproducibility

`measurements.csv` includes analysis, serialization and full process wall time,
peak R working-set bytes, sampled maximum private bytes, output bytes and thread
evidence. Process peaks cover the entire R lifetime, including input loading and
serialization. They exclude the launcher and helper processes. The memory API
is System.Diagnostics.Process.PeakWorkingSet64 on the PID emitted by Sys.getpid;
private bytes are sampled every 100 ms and may miss transient peaks. Serialized
RDS sizes include fit metadata and retained weights, not just result matrices.
Input files are inventoried once with byte lengths and MD5 hashes.

Hardware: Intel Core i7-1365U, 10 physical / 12 logical cores, approximately
15.6 GiB RAM, Windows 11 Enterprise, local C: NTFS storage. R 4.4.1 is used.
Detailed installed package/native versions and binary hashes, source revisions,
library paths and session information accompany the CSV files. Current source
revisions do not prove installed binary identity; binary hashes are recorded
separately. Build-cache records identify Release builds with Rtools44 GCC 13.2.0
and OpenMP. The archived compatibility manifests match 35 installed DLL hashes;
four gbase DLLs differ from the earlier manifest. Their exact installed hashes
are retained, but their source revision/build attribution remains unverified.
These are private validation packages, not source-release qualification.
Locale startup warnings occurred in both settings.

From the gsuite repository root:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File benchmarks/genomic_20k_50k/run_manuscript.ps1
powershell.exe -NoProfile -ExecutionPolicy Bypass -File benchmarks/genomic_20k_50k/run_manuscript.ps1 -Database
```

The driver creates fresh timestamped directories under `build/benchmarks/`.
The reported corrected run is `manuscript-20261004-142856`. The earlier
`manuscript-20261004-142639` run is retained as diagnostic evidence and excluded
from this table because its memory readings measured the launcher. The first
preflight directory `manuscript-20261004-142625` contains no benchmark results.

## Database download and installation

The separate `manuscript-20261004-143135` run used the published
`gact(version="hsa.0.0.1", task="download")` workflow in a fresh directory.
The installed gact version was 0.0.0.1. Download/install took **544.20 seconds**;
saving the GAlist descriptor took **12.48 seconds**. Full process wall time was
**575.41 seconds**, including startup, validation and inventory writing.
Peak R working set was **2.677 GiB**; sampled peak private bytes were
3,089,489,920. Helpers and the launcher are excluded from memory measurement.

Installed database files total **3,261,689,414 bytes**, including the retained
440,170,912-byte annotation ZIP and compressed summary resources. The separately
saved descriptor is **52,318,792 bytes**. The whole stage directory, including
descriptor, logs and inventories, totals 3,315,628,347 bytes. These are logical
file lengths, not allocated disk blocks or network traffic counts. The five
core Zenodo records advertise 2,757,229,929 bytes of downloadable resources;
that excludes auxiliary HTTP downloads and retries/protocol overhead.

All seven requested component categories (annotation, marker, gsea, gstat,
gbayes, atc and dgi) report success. Zenodo file checksums passed. Final checks
verified the annotation map, marker identities, nonempty installed files and
study-file paths. The study manifest identifies GWAS1 as T2DM
(PMID:30297969) and GWAS2 as CAD (PMID:26343387). The benchmark analyses still
use synthetic A/B inputs. This installation reuses published prepared resources;
it does not measure raw database construction or optional 1000G BED/LD preparation.
Network conditions contribute to the installation time. An encoding translation
warning was recorded; no installation component failed.

`database-measurements.csv`, `database-files.csv`, `database-downloads.csv` and
`database-studies.csv` retain the compact installation evidence. Detailed logs
and the usable database/descriptor remain in the ignored benchmark cache.

The optional matched official MAGMA comparison has not been executed; no local executable
was found in the targeted gsea search. It remains outside this primary table.
