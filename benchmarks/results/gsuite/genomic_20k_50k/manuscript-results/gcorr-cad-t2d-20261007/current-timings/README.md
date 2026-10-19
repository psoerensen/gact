# Current CAD/T2D tutorial timings

The four score-based stages use the current measurements with gcorr R 0.1.5
and native 0.1.4. Those numerical paths did not change in the HESS reader update.
The three HESS stages were rerun once with gcorr R 0.1.6, native gcorr 0.1.4
and the updated installed gmat physical-row reader. The previous HESS timing
rows have been replaced. The earlier evidence is preserved in the local build
backup; no earlier measurements appear in the tutorial.

All fits use the 7,221,923-marker website availability mask and extracted
stat$n. Stored LD scores are reused. Native score-based fits receive numeric
arrays and integer jackknife membership, without per-marker IDs or weights.
Chromosome sets use integer membership. Saved GAlist is never loaded by these
timing workers.

HESS resolves marker names to physical BED rows once, retains one component's
BIM metadata during descriptor validation and releases preparation buffers
before numeric decoding. The reader does not reconstruct a BedPanel or BIM/FAM
name maps. The worker releases the redundant cached input matrices before its
baseline. It recomputes 1,702 regions and 18,696 retained components with the
automatic cap (11), four threads and no saved LD or eigenvector files.
The two subsequent pair fits load these newly computed compact moments.

All regional counts, ranks, eigenvalues and quadratic/cross moments agree
with the established full-data results within 1e-10. Both pair fits' estimates
and standard errors agree within 1e-10; maximum moment change was
0. The selected native reader and
projection tests passed; the canonical R package builder also ran its existing
installed test files. No repeated numerical test campaign was run.

Each worker records a baseline after loading its required inputs and garbage
collection. Process resident/private memory is sampled every 200 ms during
computation. Loading and output serialization are excluded from times and
peaks. Peaks cover the whole process, including inputs and runtime overhead.
Resident growth is peak minus baseline, floored at zero. Private allocation
includes committed paged memory; neither measure is solver-only storage.
HESS validation took 106.47 seconds; its phase peaks are recorded in
preparation-memory-phases.csv. Native panel setup is now zero. Detailed raw
memory samples and the R heap diagnostic remain in the local benchmark run.
No build or other numerical analysis ran concurrently with the timed work.

source.sha256 records the current reviewed source snapshot. The score timings
predate the HESS-only source changes described above. The public page presents
the current figures without the internal evidence inventory.
