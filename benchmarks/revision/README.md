# Timing evidence for the manuscript revision

The [consolidated timing table](consolidated-timings.md) combines qualified earlier tutorial measurements with the [19 October gsuite measurements](../results/20261019-gsuite-revision/README.md). The new collection measures saved-input loading, current Reactome Bayesian gene-set fits, a single full-data BayesR BLR fit and separately labeled chromosome-held-out ridge/BayesC/BayesR validation, and a controlled 40-cell simulated scoring grid with separate process-memory measurements. The dated report records the completion status of the two-chain mapping and chromosome-component examples.

All new gene-set and PoPS runs reuse the already installed hsa.0.0.1 database. Scoring reuses frozen gsim-generated effects and simulated genotypes; it requires no UK Biobank data. Earlier one-chain mapping repeats, original gene-score workflows and historical qgg evidence remain distinct. These measurements describe representative updated workflows rather than reproducing original manuscript fits.

Use the workload counts, package versions, requested threads, operation boundaries, output sizes and memory definitions with every selected timing. Execution does not establish Bayesian convergence or external MAGMA/FUMA concordance. Remaining measurement gaps are listed in the consolidated table and [website overview](../../Document/Benchmark_overview.html).

Current Bayesian MAGMA/PoPS table rows come from the [direct-cross-product reruns](../results/20261019-direct-cross-products/README.md). Preparation, single full-data fits and 22-fold validation are separate; upstream gene-statistic computation is reused and excluded. Earlier records remain available as historical evidence.
