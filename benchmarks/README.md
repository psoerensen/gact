# GACT benchmark evidence

Existing timing evidence for the GACT website and manuscript revisions is collected here. Real CAD/T2D tutorial workloads and gsim-based synthetic workloads have different scopes and remain separate. The original collection was organized from saved evidence. A separate dated collection adds authorized gsuite ecosystem runs measured on 19 October 2026.

The [consolidated timing table and missing evidence](revision/consolidated-timings.md) select existing measurements for website and revision use.

## Evidence collections

| Collection | Location | Intended use |
| --- | --- | --- |
| Classical gsuite and official MAGMA reference | [20 October comparison](results/20261020-classical-magma-comparison/README.md) | Shared saved gene statistics/correlations; separate preparation/fitting, process clocks and numerical agreement |
| Current direct-cross-product reruns | [Rerun report](results/20261019-direct-cross-products/README.md) | Separate preparation, full-data MAGMA/PoPS models and 22-fold PoPS validation |
| Earlier revision measurements | [19 October report](results/20261019-gsuite-revision/README.md) | Current Reactome gene-set/PoPS, 40-cell scoring grid and additional module checks |
| Website timing source snapshot | [Tutorial timings](scripts/website-snapshot/Tutorial_timings_and_memory.Rmd) | Published timing definitions and reported workloads as imported |
| Simulated 20K by 50K workflow | [One thread](results/gsuite/genomic_20k_50k/results/summary.csv), [four threads](results/gsuite/genomic_20k_50k/results/summary-4cores.csv) | Existing gsim workflow including score application to 2,000 held-out recipients |
| Complete CAD/T2D VEGAS | [Report](results/gsuite/genomic_20k_50k/manuscript-results/manuscript-vegas-complete-20261005/README.md) | Preparation, first fit, reuse costs and historical qgg qualifications |
| Simulated fine mapping | [Dense route report](results/gsuite/genomic_20k_50k/manuscript-results/manuscript-simulated-finemapping-dense-20261006/README.md) | Three serial repeats, model settings and memory |
| Current correlation timings | [Report](results/gsuite/genomic_20k_50k/manuscript-results/gcorr-cad-t2d-20261007/current-timings/README.md) | Current tutorial timing and memory definitions |
| Historical qgg and standalone runs | [Results](results/gsuite/genomic_20k_50k/manuscript-results) | Earlier evidence, including superseded and blocked runs; consult each report |

## Organization

- `scripts/gsuite/genomic_20k_50k/` preserves existing driver and summarizer files unchanged.
- `scripts/website-snapshot/` preserves selected tutorial sources. Live tutorial sources remain in `../Document/`; these are dated evidence snapshots.
- `protocols/` preserves original measurement plans and limitations, including explicitly superseded plans.
- `results/` holds existing small tables, reports and provenance records. Directory names and filenames retain original run identifiers.
- `revision/` identifies reusable evidence and remaining reporting gaps for this manuscript.
- `provenance.csv` records original absolute paths, modification times, byte sizes and SHA-256 checksums for every imported file.
- `retained-at-source.csv` inventories omitted large and binary artifacts at their original locations.
- `import-metadata.json` records import-time repository state, distinct from measured-run versions.

## Execution and preservation

Drivers are archival copies, not newly configured GACT launchers. They retain original gsuite-relative paths, libraries, data dependencies and output locations. Do not invoke them from the GACT root. Read the original protocol and prepare a separate writable run directory before any future authorized execution. No installation or analysis is required to inspect these tables.

Original files and large evidence archives remain in gsuite. Relative links inside copied reports retain their original layout assumptions; consult the provenance manifest for source locations when an asset is not included. A source archive checksum is not verification of an omitted archive's current contents.

Store future genotype panels, individual data, simulations, LD caches and large binary outputs outside Git. Preserve old measurements; add a dated run directory for new executions. Capture measured code/package versions, controls, hardware, workload counts, thread policy, timing boundaries, memory definition and completion status with each run.

Website and revision tables should select explicit evidence snapshots. Do not merge numerical results from different estimators or label a tutorial timing as a manuscript reproduction. The import does not establish convergence, statistical calibration, or matched performance across different methods.

## Local installed database

The already installed hsa.0.0.1 database is copied to `benchmarks/data/hsa.0.0.1/`, with a relocated `GAlist_local.rds` descriptor. All 96 original files were verified byte-for-byte using SHA-256. The original installation remains available. The data directory, descriptors and raw benchmark outputs are ignored by Git and are not website assets. Copy provenance is retained locally in `benchmarks/data/installation-provenance.json`; the dated evidence report records its checksum and storage totals.

The [prepared-summary BLR report and implementation plan](results/20261019-gsuite-revision/summary-blr-qgg-settings/README.md) separate sparse cross-product preparation from the native sampler and explain the missing direct R interface.
