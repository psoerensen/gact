# Classical gsuite and official MAGMA: saved comparison records

These completed 20 October runs use **17,862 saved genes and 1,324 Reactome pathways per trait**, with the same saved MAGMA gene statistics, gene properties and gene correlations in both classical pathway implementations. The underlying matched inputs used 3,192,403 SNPs and a 503-person EUR reference. This is a **shared-input classical pathway check**, separate from the website's Bayesian models, which use a different 17,955-gene score collection. **No gene statistics were recomputed for these gsuite fits or this website update.**

### Gene-based analysis: saved gsuite VEGAS and MAGMA timings

| Trait scope / operation | gsuite VEGAS analysis clock (s) | Official MAGMA gene-process clock (s) |
| --- | ---: | ---: |
| GWAS1 initial analysis | Not timed separately | 4333.469 |
| GWAS2 initial analysis | Not timed separately | 4029.865 |
| Both traits: initial analysis | **899.94**, shared preparation and initial gene analysis | **8363.334**, sum of the two processes above |
| Both traits: alternative prepared refit | **71.44**, prepared-object loading excluded | Not timed |

The saved gsuite VEGAS measurement covers 17,862 genes and a 503-person EUR reference. Its **899.94 s is an analysis clock, not a whole-process or fitting-only clock**: it includes shared regional spectra, prepared-object conversion, alignment, initial gene analysis and intermediate checkpoint/cleanup overhead. Input-file loading (**110.45 s**) and neighbouring-gene correlation preparation are excluded. Preparation uses four decoding threads while regional products/eigensolvers are serial. Recorded components include **643.36 s** for spectra, **79.52 s** for conversion and **97.69 s** for the initial two-trait VEGAS call; these overlap the complete analysis and must not be added to it.

The **71.44-second refit is an alternative reuse measurement**, not an additional initial-analysis stage. It excludes a separately measured **17.26 s** to load the prepared object. Per-trait gsuite gene times are unavailable; the two-trait total has not been split or duplicated. The archived reference contains 9,985,562 physical markers, whereas 3,192,403 SNPs were selected for the later matched MAGMA input; physical reference size and selected GWAS SNP count are different quantities.

**The timing boundaries are different.** MAGMA's gene-process clocks include startup, loading, gene analysis, gene-correlation work and output. Standard gsuite and MAGMA also retain different gene-tail algorithms. The earlier VEGAS clock is not a matched identical-work speedup comparison with MAGMA. These gene-stage records and the pathway-stage records below come from separate saved runs and must not be added into a purported measured end-to-end workflow. No analyses were rerun to add this table.

See the [saved VEGAS report](../gsuite/genomic_20k_50k/manuscript-results/manuscript-vegas-complete-20261005/README.md) and [component clocks](../gsuite/genomic_20k_50k/manuscript-results/manuscript-vegas-complete-20261005/phase-summary.csv).

### Classical gene-set analysis: separate preparation and fitting

| Trait | gsuite preparation (s) | gsuite fitting (s) | gsuite process wall (s) | MAGMA gene-set process wall (s) |
| --- | ---: | ---: | ---: | ---: |
| GWAS1 | 32.691 | 0.0226 | 50.299 | 29.451 |
| GWAS2 | 31.123 | 0.0195 | 32.395 | 79.332 |

Preparation includes set mapping, gene-response/property transformations, powered regional correlation eigensystems/generalized inverses, internal-covariate adjustment, and cached shared/member-index cross-products. Fitting uses those sufficient statistics for the small pathway regressions and inference. The prepared weights are reused across all pathways: the marginal preparation uses 33 shared block applications and 223,548 within-set precision entries. The two traits were run as separate processes, so regional preparation was measured twice; no cross-trait reuse is claimed here.

These are **native C++ gsea/gsolve Release/MKL measurements with a one-thread limit**. They do not measure R interface conversion overhead. Input loading is separate: 0.849 s for GWAS1 and 0.850 s for GWAS2. The process clock also includes startup, output and cleanup; it is not the sum of the two displayed phase columns. An additional **16.736 s (GWAS1) and 0.403 s (GWAS2)** falls outside all internal loading/preparation/fitting clocks and has not been attributed to individual causes. These excess times are retained in the process totals. Peak memory was not measured for these classical runs.

The official MAGMA 1.10 gene-set clocks include loading, correlation preparation, regression and output. **Separate MAGMA preparation and fitting clocks are unavailable**; its whole-process time must not be compared with gsuite fitting alone. The saved MAGMA gene-analysis clocks are a separate upstream stage, include gene-correlation work, and are not added to the gene-set clocks. The earlier 899.94-second gsuite two-trait gene analysis excludes file loading and gene correlations; it is not an identical-work speed comparison with these MAGMA gene-process clocks.

### Numerical agreement and interpretation

| Trait | Maximum absolute beta difference | Maximum absolute SE difference | Maximum absolute p difference |
| --- | ---: | ---: | ---: |
| GWAS1 | 3.731e-05 | 4.994e-06 | 5.003e-06 |
| GWAS2 | 3.806e-05 | 4.997e-06 | 5.005e-06 |

**All 2,648 pathway fits pass the unchanged 2e-4 absolute beta/SE/p tolerance** against the rounded MAGMA text outputs; p-value Spearman agreement exceeds 0.99999998 in both traits. Gene backgrounds and pathway IDs match, with zero pathway-size mismatches. The working correlation power is 1.1 and the absolute eigenvalue cutoff is 0.1, retaining rank 13,814 and marginal residual degrees of freedom 13,812. All gene rows remain in the background; pruning removes correlation directions.

The direct-set refactor differs from the previous correlated gsuite route by less than 2e-15 in beta and p (less than 3e-16 in SE). Its fitting phase replaces repeated dense products with cached cross-products; the shared work is now included in preparation. Previous preparation/fitting times were 27.558/88.300 s and 27.527/86.257 s. Whole-process times changed from 117.195 to 50.299 s and from 115.025 to 32.395 s. These are single runs, not estimates of repeated-run uncertainty or a universal speedup.

This agreement isolates the classical pathway model on **shared saved MAGMA statistics**. It does not establish equivalence between the different gene-tail algorithms used by standard gsuite and MAGMA, or validate Bayesian MAGMA/PoPS against the external executable. No calibration or convergence claim follows from these timings.

## Evidence and reproduction

This collection imports completed measurements; no analysis was executed to update the website.

- [Exact timing records](measurements.csv) and [agreement summary](agreement.json).
- [Input fingerprints](input-provenance.json), [measured code fingerprints](code-hashes.json), [collection provenance](provenance.json) and [completion status](status.json).
- [Existing native comparison driver](classical_driver.cpp) and [CMake target](CMakeLists.txt) use the installed public gsea interfaces. They are archival source, not instructions to launch an analysis during a documentation update. A compatible native gsea with the direct-set refactor, installed gsolve and the original frozen inputs are required for any separately authorized reproduction.

The imported source run is `matched-magma-gsuite-20261019/classical-direct-sets-20261020`. Original gene-statistic/correlation files, per-pathway fits and binary outputs remain in the local run directory. Their hashes are retained here, not their contents. The previous correlated gsuite times are historical outputs, not newly rerun baselines. No official MAGMA source or executable is redistributed.

Preparation/fitting were measured in the native driver, not the R interface. The updated R interface passed focused numerical checks separately. Five native CTest checks passed; the R correlation oracle and dispatcher checks passed (14 and 40 expectations).
