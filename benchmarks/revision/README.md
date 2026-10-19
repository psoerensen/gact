# Timing evidence for the manuscript revision

Use the updated website workflows for representative computational requirements, with existing gsim-based score application as synthetic evidence. Original manuscript-result provenance is a separate question.

| Operation | Existing evidence | Qualification |
| --- | --- | --- |
| Database download and setup | 544.20 seconds and 2.68 GiB in the website snapshot | Network transfer and installation of an existing database, not reconstruction from raw sources |
| CAD/T2D VEGAS | 899.94 seconds, 17,862 genes, 5.88 GiB | Complete initial analysis; prepared refit has separate clocks; different tail evaluator from historical qgg |
| Correlation workflows | Current timing report and website-data CSVs | Use per-row measurement and estimator definitions |
| Simulated regional and genome fits | Dense simulated fine-mapping report | Short fixed settings and one chain in timed runs; not manuscript BayesR fits |
| Simulated score application | 14.878789 seconds at one thread and 5.309708 seconds at four | 2,000 recipients, 50,000 markers, total plus 20 artificial sets; memory is not supplied for this combined run; separate scoring workloads with memory are available in the consolidated table |
| Gene-set fitting | Historical qgg and standalone joint results | Verify implementation, KEGG versus Reactome, model controls and version before selecting a revision row |

The current benchmark summary should state hardware/software, workload, threads, preparation versus fitting, memory definition and output sizes where available. Preserve missing fields rather than substituting machine RAM or unrelated run peaks. Historical comparisons and older failed/blocked runs are retained as supporting evidence, not the headline current table.

No new timings, fitted models or synthetic datasets were generated during organization. No UK Biobank scoring benchmark is required for this collection. The dedicated benchmark overview and timing-page navigation have been rendered locally from existing evidence; publication remains pending.

The [updated consolidation](consolidated-timings.md) includes recovered separate scoring time/memory and the narrower remaining gaps. The [website overview](../../Document/Benchmark_overview.html) presents these records with workload qualifications.
