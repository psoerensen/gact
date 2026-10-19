# Real-1000G-derived performance workload

Authorized on 4 October 2026. Previous artificial-data measurements and the
completed database installation are preserved separately.

Use the verified cached official 1000G Phase 3 chr22 VCF, official GRCh37 map
and population panel. Select EUR donors by `super_pop`, import interval
22:20000000-22200000, and report actual marker count without padding to 50,000.
EUR selection differs from the earlier gsim report's first 500 VCF samples.

Generate a 20,000-person reference-conditioned cohort using public gsim:
EUR weight 1, N=donor count, Ne=10000, rho=.02, mutation age 1e9,
seed=20260930, batches of 4096, one thread. The public pedigree interface requires
two generations: the all-new-founder setting produces 19,998 founders and two
offspring. Convert these haplotypes to BED through that public interface. Preserve
physical coordinates and official map distances. These are synthetic individuals
derived from real reference LD, not original 1000G participants.

Retain markers with cohort MAF >= .01. Simulate two continuous A/B phenotypes
with fixed sparse effects (causal probability .02, h2 .3/.4). They are not actual
CAD/T2D phenotypes. Use 6,000 disjoint GWAS samples per trait, 18,000 LD samples
and 2,000 held-out score recipients. Prepare signed LD without r2 pruning and
with a 273-variant window: complete pairwise coverage for each 274-marker set.
The completed cohort retains 5,360 of 48,088 imported markers at this threshold;
the analysis workload is not a post-QC 50K-marker panel.

Measure one-thread VEGAS on 1,000 distinct contiguous 274-marker windows.
Size 274 matches the raw Genesplus median before reference/GWAS filtering.
These overlapping windows are size-controlled timing sets, not biological
annotations or a matched post-QC gene universe. Joint Bayesian pathways consume
this measured evidence and the same three synthetic membership patterns, using
four fixed chains with 200 burn-in and 800 sampling sweeps. BayesR mapping uses
20 distributed 100-marker regions per trait, four fixed chains, 500 burn-in,
700 sampling sweeps, and the previous fixed prior/eigen controls. Mapping region
size is separate from the VEGAS gene-size requirement.
Float32 LD requires VEGAS relative PSD tolerance 1e-7 and mapping eigen
verification tolerance 1e-7. Declare GRCh37 explicitly on the mapping Glist.

Construct fixed sparse simulation-effect weights, normalized to unit L2 norm,
then measure total and pathway scoring on 2,000 recipients with standardized
markers and unscaled scores. These are not posterior weights. The attempted
auxiliary genome-wide BayesC fit failed with finite-range arithmetic overflow
on window-truncated LD; preserve it as excluded diagnostic evidence. Old
artificial marker weights cannot be reused. Record import, founder generation,
BED conversion, phenotype/GWAS/LD preparation and weight construction separately
from primary analysis timings. Serialization is separate;
peak memory covers the whole stage process lifetime.

Run the existing measurement launcher with `-RealLD`. Fresh outputs live under
gsuite `build/benchmarks/`; gsim and qgdata inputs remain read-only. Stop on stage
failure and preserve evidence. Record hashes, installed software/native
provenance, controls, dimensions, finite-output checks, mapping diagnostics and
sizes. Gene-count extrapolations remain estimates at fixed marker/reference
size; matching a median does not measure the large-gene tail.
