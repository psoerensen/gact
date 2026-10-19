# Regional eigen preparation: focused development evidence

Date: 2026-10-05. Working changes based on gsuite `4d2bf0aa3641a4fb4d36e827f78f61dc08fdda80` and
gmat `4e260c9429adba58c2c0116af15736e6092744ba`. The archive includes exact changed source files (including
new files), tracked patches, comparison script, spectra and validation logs.
This evidence is separate from the earlier download and sparse-LD timings.

Existing chromosome 22 reference: 503 samples, 141,121 physical markers.
Both prepared descriptors use the identical BED/BIM/FAM fixture in
`build/benchmarks/manuscript-chr22-ld-20261005-200555/reference_fixture-1`.
Regions select the first variable markers using gbase's stored allele frequencies;
the batch uses ten distinct consecutive groups of 1,000 such markers.
No data download or new sparse-LD computation was timed.

## Scope and observations

`timings.csv` reports one development run per case, without confidence intervals.
The gbase measurements include descriptor validation, panel construction, regional
genotype decoding, dense product and eigenvalues-only decomposition. qgg includes
`getG(scale=TRUE)`, `cor()` and either values-only eigen or its current full eigen
calculation. Package loading and descriptor RDS loading are excluded. No external
thread or BLAS setting was changed.

- Ten 1,000-marker regions in one gbase call, four decoder threads: **2.34 s**.
- Corresponding current qgg computation: **6.18 s**.
- Maximum absolute eigenvalue difference over all ten regions: **1.1653e-12**.
- A standalone 1,000-marker gbase call: **1.58 s**, versus **0.69 s** for current qgg.
- The 100-marker case is likewise slower as a standalone gbase call.
- One- and four-thread gbase spectra were identical in the individual-region checks.

Pass all region sets in one `gprep(Glist, task="eigen", sets=regionSets)` call
to share descriptor validation and panel construction. These small measurements
support the batched workflow, not a universal speed claim or a paper-wide timing
replacement. The decoder is threaded; the current portable gsolve product and
Eigen eigensolver are serial. No OpenBLAS dependency was added.

When marker count exceeds sample count, gbase uses a 503-by-503 sample-space
Gram matrix and pads the spectrum with exact null zeros. It keeps the full
marker-length eigenvalue vector, without clipping or a sparse-LD cutoff.
The descriptor stores selected samples, marker memberships and diagnostics.
Current glma VEGAS does not yet consume this new resource.

## Validation

- Focused gmat BED-panel test passed, including regional normalization and errors.
- All gbase tests passed (122 expectations, including 32 new eigen expectations).
- Source R CMD check: no errors or warnings; two NOTES for bundled DLL size and
  native-library termination symbols. This is gbase validation, not complete
  gsuite ecosystem qualification.
- Allowlisted documentation build: 7,413 references checked in 92 HTML pages,
  zero failures; existing local gtools snapshot refreshed without committing.

Archive SHA256: `22fe610d1e0f9c5807a8ae4dae9c97b26844d03d923fc0f2b1522b0d280f1cc5`.
