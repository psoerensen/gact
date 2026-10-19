# Published gact/qgg website workflow, EUR reference

Primary implementation: installed qgg 1.1.6 and gact 0.0.0.1, one thread.
These measurements are separate from the earlier standalone-package synthetic
and 1000G-derived simulation studies, whose results are preserved unchanged.
The execution follows [the agreed protocol](../../qgg-examples-protocol.md).

## Completed preparation

| Operation | Operation wall time | Peak process working set | Output size |
| --- | ---: | ---: | ---: |
| Previously downloaded gact database | 544.20 s | 2.677 GiB | 3,261,689,414 bytes of database files |
| EUR reference download/extraction | 134.63 s | 1.112 GiB | 511,626,945-byte archive; 2,855,798,067-byte BED |
| Filter reference and prepare Glist | 1,399.86 s | 5.430 GiB | 1,258,180,815-byte filtered BED |

Database descriptor serialization added 12.48 s and 52,318,792 bytes.
Glist serialization added 148.05 s; its complete process took 1,588.99 s.
Download files passed the downloader validation, but that process had a footer
parse error after saving its outputs/timing. Its operation timing is retained
with this qualification; it is not a successful full-process execution.

The EUR reference contains 503 individuals and 9,985,562 retained markers.
The archive MD5 is 1919cb5c79bbe7871aed71ae4abe6217 (Zenodo DOI
10.5281/zenodo.14141438). Filtered BIM/FAM sizes are 294,619,969/13,078 bytes.

## Incomplete stages and resource limits

The published LD window, msize=1000, requires 79,924,438,248 bytes
(74.44 GiB). C: had 47,847,862,272 bytes available at the stage check.
LD computation did not start. Its failed process duration must not be reported
as the time to calculate LD. Whole-genome fitting and total/pathway scoring
remain pending that prerequisite. No regional or cached-weight substitute was
introduced.

Installed `gbayes` calls `getSparseLD` with its default r2=0, which reads all
reference LD rows rather than only the QC-retained GWAS rows. If nearly all
window entries are nonzero, 9,985,562 x 2,001 entries require approximately
298 GiB for eight-byte values and eight-byte indices, before list/object
overhead. This is a source-based capacity estimate, not a measured allocation;
actual sparsity has not been measured because LD computation was blocked.
The machine has 16 GB RAM, so additional disk alone does not establish that
the unmodified whole-genome fitting example can run here.
The installed onebased=FALSE branch subtracts a numeric 1 from the integer
sequence, promoting indices to R doubles; the payload estimate accounts for
that promotion.

The independent VEGAS, joint CAD/T2D Bayesian MAGMA and four-gene mcmc-eigen
examples are continuing from the successful prepared Glist. Their measurements
will be added after completion and output validation. No standalone MAGMA
comparison has been executed.

## Gene-size audit during the long VEGAS run

The supplied 18,279 gene sets have raw marker-count median 274, 95th percentile
1,142, and maximum 20,376 (CSMD1, ENSG00000183117). The next largest is 17,739
(RBFOX1, ENSG00000078328), followed by PTPRD (12,772). Symbols are from the
installed database's ensg2sym.rds. These are before GWAS/reference QC; effective counts will
come from the completed VEGAS output. Installed VEGAS calls eigen(cor(B)) for
each retained gene, so a small number of very large genes can dominate runtime.
This is a supported explanation, not a measured attribution to a particular
gene: the run has no per-gene clock/progress record. Raw sizes are saved in
raw-gene-size-audit.csv. A brief separate size/source audit ran concurrently;
these are descriptive single-run measurements on the user's workstation.

## Measurement definition and provenance

Operation clocks exclude previously loaded inputs and final serialization.
Whole-process clocks include startup, loading, computation and saving. Peak
working set is the Windows process high-water mark; sampled private allocation
is recorded separately. These are single-run measurements, not medians or a
scaling study. Output bytes refer to actual retained files.

Hardware: Intel Core i7-1365U, 10 physical/12 logical cores, 16,794,288,128 bytes
RAM; Windows 11 Enterprise 10.0.26200. Software: R 4.4.1 UCRT x64, qgg 1.1.6,
gact 0.0.0.1. R session reports and installed argument defaults are preserved
per stage. OMP, OpenBLAS and MKL thread limits were set to one.

Raw runs: build/benchmarks/manuscript-20261004-143135 (database),
manuscript-20261004-171036 (EUR download diagnostic),
manuscript-20261004-171426 (Glist and storage gate), and
manuscript-20261004-174351 (independent gene analyses).
