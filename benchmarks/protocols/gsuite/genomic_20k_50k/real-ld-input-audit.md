# Existing gsim real-reference inputs

Read-only input audit, 4 October 2026. No simulation or analysis benchmark ran
during this audit. The local gsim documentation confirms earlier real 1000G
reference-based experiments. No new download is needed to recover their inputs.

## Instructions and earlier measurements

- `../gsim/docs/README.md` identifies the real-reference production benchmark.
- `../gsim/tools/benchmark/benchmark_1000G_chr22.R` separates VCF import,
  reference opening, founder generation and pedigree stages. It requires explicit
  external data/output directories; downloads are opt-in.
- `../gsim/docs/qualification/production_benchmark_chr22.md` documents 500
  reference donors, region 22:20000000-22200000 (GRCh37), 48,088 retained markers,
  and the official recombination map. Its original founder run used 1,000 founders.
- `../gsim/docs/qualification/parallel_founder_benchmark_chr22.md` documents
  10,000 founders on the same real reference, with identical one/eight-thread output.
- `gsim_import_vcf()` followed by `gsim_simulate_founders()` is the supported
  real-reference path. Phenotype simulation alone does not generate real-reference LD.

## Located local inputs and fresh checks

| Input | Existing path | Verified SHA-256 |
| --- | --- | --- |
| Official phased chr22 VCF | `../qgdata/cache/1000g_phase3/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` | `a90c16c4ff2b3196476d506ae13cb3047fae8670163c7c932c4b0239aef3daf5` |
| Official GRCh37 chr22 map, gzip | `../qgdata/cache/1000g_phase3/genetic_map_GRCh37_chr22.txt.gz` | Uncompressed content: `fe688fb1722c0f2aaeb691d813ad5d759fd033e345369e50694f5cb04b12d8a8` |
| Existing 5,000-individual / 50,000-marker EUR mosaic BED | `../qgdata/human_1000g_eur/human_1000g_eur.bed` | `39e8ff00f11ebafe8d39e3d602b08df21dfd38f5972af41176d39261bdef6dce` |

The VCF and uncompressed map hashes match the gsim production report exactly.
The map was decompressed without changing bytes into the gsuite ignored
`build/benchmarks/gsim-reference-audit/` cache; original inputs were unchanged.
The EUR BED hash matches its existing README. That panel was made by qgdata's
base-R copying mosaic, not gsim; it must not be labelled a saved gsim output.
Its BIM/FAM hashes have not been freshly verified in this audit.

The documented large gsim output paths are environment-selected external
directories, not committed artifacts. The targeted searches did not locate
those prepared outputs, but the original input VCF and map are present and verified.

## Revised workload direction

Use the cached official reference/map through gsim's documented founder path,
retaining physical positions and genetic distances. Do not remap artificial
100-marker chromosomes. Preserve all existing benchmark results as artificial-LD
evidence and the completed database-installation evidence separately.

For the revised performance run, settle the exact founder count, ancestry/donor
selection and marker extent before execution: the documented real-reference
fixture has 48,088 markers, whereas the previous artificial fixture has 50,000.
The gsim production report's first 500 VCF samples are not automatically an EUR
selection. The population panel is available in the same qgdata cache.

Use 1,000 gene-sized sets with median marker count matched to the relevant
real input after alignment, or an explicitly size-controlled 274-marker target
if the user chooses the raw Genesplus median. Preserve real LD in either case.
Gene-count extrapolation should be labelled an estimate at fixed set size and
reference size; the large-set tail and preparation overhead are not captured by
matching a median alone.
