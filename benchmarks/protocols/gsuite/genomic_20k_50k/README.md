# Representative 20K-by-50K genomic workflow timing

This benchmark records one descriptive execution of the connected genomic
workflow on 20,000 simulated individuals, 50,000 markers and three continuous
traits. It is intended for a reviewer-facing runtime table, not for method
ranking or formal performance claims.

The run generates one packed PLINK dataset, prepares one 18,000-individual LD
reference, computes sparse LD once and reuses it throughout. The remaining
2,000 individuals are held out for score application. Three disjoint groups of
6,000 individuals supply the association summaries.

Run from the repository root after building the qualified task-package
ecosystem:

```powershell
& 'C:\Program Files\R\R-4.4.1\bin\Rscript.exe' --vanilla benchmarks/genomic_20k_50k/run.R
```

The OpenMP-capable Windows binaries can be exercised on the same frozen data
with an explicit thread limit:

```powershell
& 'C:\Program Files\R\R-4.4.1\bin\Rscript.exe' --vanilla benchmarks/genomic_20k_50k/run.R --threads=4
```

This writes separate `*-4cores.csv` results and a separate intermediate cache.
It does not assume that every stage benefits from OpenMP; the measured thread
limit and per-stage timings must be reported together.

The qgg comparison reuses that frozen simulation. Install the qgg source being
measured in the isolated `build/benchmarks/genomic-20k-50k/qgg-library`
directory, then run:

```powershell
& 'C:\Program Files\R\R-4.4.1\bin\Rscript.exe' --vanilla benchmarks/genomic_20k_50k/run_qgg.R
```

Package compilation and installation are outside the timed stages. The qgg
driver records its source commit and uses separate Glist and LD files. Its
one-thread run evaluates four sequential BayesC chains to match the chain count.

Large intermediate data are written below
`build/benchmarks/genomic-20k-50k/`. Small tables intended for review and later
website publication are written to `benchmarks/genomic_20k_50k/results/`.
Passing `--resume` reuses completed stages; it is for recovering an interrupted
run and must not be described as a new timing replicate.

The reported stages are data simulation, Glist preparation, sparse-LD and
LD-score preparation, three linear GWAS scans plus `checkStat()`, ordinary and
partitioned LDSC, VEGAS for 1,000 artificial genes, MAGMA-style pathway
analysis and PoPS-style ranking, one BayesC fit, multi-effect fine-mapping of 20
regions across two traits, and total plus marker-set score application.

The columns describe comparable workflow goals, not a claim that every method
is identical. In particular, qgg's partitioned LDSC route estimates
heritability rather than genetic correlation, qgg PoPS consumes set-valued
features, and qgg fine-mapping uses BayesC with `mcmc-eigen` while the measured
gsuite route uses multi-effect eigen mapping. These distinctions are retained
in the downloadable table and must accompany any performance interpretation.

## Tutorial verification against the source-release packages

After the 20K/50K simulation has been generated and the seven 0.1.0 source
packages have been installed in `build/source-repository/library`, run:

```powershell
& 'C:\Program Files\R\R-4.4.1\bin\Rscript.exe' --vanilla tools/validation/verify_package_tutorials.R
& 'C:\Program Files\R\R-4.4.1\bin\Rscript.exe' --vanilla tools/validation/verify_package_tutorials.R --full-clumping
```

The first command re-executes the tutorial's main analysis routes and focused
conditional, ORA, colocalisation and C+T examples. The second creates a
window-complete LD reference for the full 50,000-marker C+T example. It uses
the saved simulated PLINK files but does not reuse saved analysis fits, build
or install packages, or count as another timing replicate. Results are saved
in the ignored four-core cache as `tutorial-verification.csv` and
`tutorial-verification-clumping.csv`. This is a Windows release-package
walkthrough, not cross-platform source-install qualification or a calibration
study.
