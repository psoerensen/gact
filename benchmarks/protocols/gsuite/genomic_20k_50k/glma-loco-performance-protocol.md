# Standalone LOCO performance checks

Use `run_manuscript.ps1 -GlmaComparison` with an explicit frozen fixture,
`-OnlyStages glma_loco` or `glma_prior`, a thread count and a timeout. Every
fit runs in a fresh R process. The harness records analysis time, process CPU
time and peak resident memory. Timeouts preserve evidence and are incomplete
fits; they are not successful timings.

Current `glma_loco` measurements use the replacement Gaussian-only joint
engine (`infinitesimal_loco`, native glma 0.1.1). It uses fixed marker power -1,
one ridge prior, blocked RHE initialization and a shared three-point Monte
Carlo heritability grid. Global variance estimates stay fixed for chromosome
refits. The default grid uses ten probes on this fixture. Old variance-search
controls and projection algorithms are rejected by the standalone interface.
Keep automatic structure scaling enabled for the full frozen measurement.
Record the estimator and loaded DLL versions to distinguish these runs from
the historical `glma_loco` results below. The replacement uses a different
estimator; old timings and association agreement do not establish equivalence.

The website fixture has 5,000 individuals, 37,991 retained markers and 22
chromosomes. Its phenotype and marker/sample order were replayed and checked
against the frozen qgg example. Reuse the fixture from
`build/benchmarks/manuscript-20261006-123111/glma_fixture-1/output.rds`.
Do not regenerate the phenotype for a timing comparison.

For a completed before/after comparison, the smaller native diagnostic input
has 1,000 individuals, 10,000 markers and three chromosomes. Its phenotype and
additional fixed covariate must both be retained when preparing the R fixture.
The preserved baseline library includes both glma and gbase and their native
DLL closures. Select it through `GLMA_BENCH_LIBRARY`; an empty value selects
the updated isolated standalone installation.

The historical variance-search measurements retained all background markers, 15 simulation pairs,
the existing Monte Carlo seeds, variance-ratio convergence criteria,
calibration settings and the default prior grid. Their scan block remained 256;
updated variance-search fits used ridge chunks no wider than 128 so their Gram
cache covered this input. Those one/four-thread measurements are preserved.

Inspect numerical agreement outside the timing process. Compare association
statistics, calibration and chromosome variance estimates for the small
before/after fits, and selected priors and associations across thread counts
for prior-selected fits. The independent native tests also check ridge
solutions and association calibration against dense reference calculations.

The earlier user-stopped full-data infinitesimal run consumed about 88 minutes
of CPU without returning a result. It is diagnostic evidence of slowness,
not a completed baseline or a denominator for an exact speedup claim.
