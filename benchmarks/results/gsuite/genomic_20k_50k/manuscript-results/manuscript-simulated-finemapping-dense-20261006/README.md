# Simulated-human fine-mapping: dense regional route and familiar interfaces

Fresh regional eigen profile: **116.00 s**. New dense regional median:
**5.48 s**, a **21.17×** reduction on this
fixed four-region workload. Historical pre-milestone measurements (73.77 s and
122.66 s) remain in their separate 20261005 reports; they used the older solver
runtime. The fresh profile uses the same current runtime as the dense repeats.

Three serial repeats, one thread, unchanged frozen inputs: 5,000 actual people,
50,000 physical markers, 37,991 QC-retained markers. Regional sizes are 850,
785, 1091 and 1399. Statistics n=4998 is the regression residual degrees of
freedom. Regional runs use 1000 post-burn sweeps and 100 burn-in sweeps; genome
runs use 100 sweeps with no burn-in. Download, simulation, QC and whole-genome
sparse LD preparation are excluded from these fit timings and reused unchanged.
Serialization and process startup are recorded separately.

| Stage | Median operation | Range | Median process peak RSS |
| --- | ---: | ---: | ---: |
| cmp_qgg_regional | 6.14 s | 5.98–16.91 s | 0.84 GiB |
| cmp_qgg_genome | 9.96 s | 9.80–10.56 s | 2.72 GiB |
| cmp_qgg_credible | 3.78 s | 3.43–3.98 s | 0.66 GiB |
| cmp_suite_regional | 5.48 s | 5.22–6.28 s | 0.56 GiB |
| cmp_suite_genome | 2.36 s | 2.03–2.38 s | 0.81 GiB |
| cmp_suite_credible_oracle | 3.95 s | 3.83–4.16 s | 0.52 GiB |

The legacy stage name `cmp_suite_credible_oracle` now measures standalone
`gmap::getCredibleSets()`, not qgg postprocessing. An untimed audit applied qgg's
helper to the **same standalone PIPs** and verified exact group membership in
all repeats. Full regional LD was independently compared to qgg at tolerance
1e-12. Original BED/BIM/FAM MD5 checks passed. Thirteen invalid genetic-map
entries were corrected only in the separate owned BIM metadata copy; physical
coordinates, alleles and source files remain unchanged.

| Region | Markers | qgg causal PIP | Dense causal PIP | PIP correlation |
| --- | ---: | ---: | ---: | ---: |
| 19:55551662 | 850 | 1.000 | 1.000 | 0.99899 |
| 11:935343 | 785 | 1.000 | 1.000 | 0.99808 |
| 13:92838400 | 1091 | 1.000 | 1.000 | 0.99779 |
| 11:92296987 | 1399 | 1.000 | 1.000 | 0.99635 |

All four genome-wide causal PIPs also remained one in each engine. These are
fixed-input checks: all regional PIPs matched the fresh standalone eigen profile
exactly in all three dense repeats. The matrix representation and grouping
implementation changed; this is not an assertion of identical effect draws.
These are
single-chain demonstrations, not convergence, calibration or posterior
equivalence evidence. Dense sampling preserves complete SPD LD with no jitter,
repair or truncation. Native gbayes sampling replaces the eigen operator;
floating-point behavior and RNG consumption can change. qgg's independent vg
draw is not recreated. Genome-wide standalone fitting retains its standardized
summary likelihood and reference-LD residual adjustment .9; qgg uses supplied
marker-specific ww/wy and neighborhood adjustment. Native genome bm remains in
standardized units. Marginal-PIP group mass is expected active-marker count,
not component credible coverage. Timing improvements apply to this workload.

The regional phase CSV separates decoding, LD product, preparation (SPD check
or eigen decomposition plus verification), binding, native sampling and set
grouping. For the eigen baseline its LD total includes separate regional sparse
file preparation; decoding/product subphases were not instrumented there.
Exact priors, effective controls, diagnostics, package/native version
information and session records are retained in output RDS files and archives.
Native gmap is 0.3.0 ABI2; the refreshed gsolve runtime is 0.26.0 ABI4.

Validation: all focused gbase/gmap/gbayes package tests and clean binary loads;
three focused native gmap tests; full-LD and same-PIP grouping audit. The
canonical package build stopped on an oracle-test input error after creating
the unchanged gbayes binary; corrected tests were rerun successfully against
that binary, followed by a fresh clean-path load. Source
checks were explicitly skipped by the selected canonical package workflow;
this is not complete R CMD check or suite-wide qualification.

Reproduce with `benchmarks/genomic_20k_50k/run_manuscript.ps1 -SimComparison
-QggInterfaces -QggCompareRun build/benchmarks/manuscript-20261005-231349
-QggReuseRun build/benchmarks/manuscript-20261005-231155 -OnlyStages
cmp_qgg_regional,cmp_qgg_genome,cmp_qgg_credible,cmp_suite_regional,cmp_suite_genome,cmp_suite_credible_oracle,cmp_compare`.
For a fresh eigen profile omit QggInterfaces and select cmp_suite_regional.

`source-and-evidence.zip` includes exact installed R/DLL copies, native gmap
sources, interface sources, benchmark scripts, all new fits and phase records,
the frozen statistics and corrected BIM copy. Original genotype files and
large sparse-LD caches stay in place and are not duplicated into the ZIP.
Source commits plus diffs record uncommitted development state. Historical
archives are preserved. No commit, push or public publication was performed.
