# qgg/gact website benchmark protocol

Authorized 4 October 2026. Use installed qgg 1.1.6 and gact 0.0.0.1 as the
primary manuscript implementations. Preserve prior standalone evidence and the
completed database installation audit. Do not modify source databases or inputs.

The downloaded database contains summary statistics, gene/pathway annotations,
saved VEGAS results and fitted BLR results, but its Glist directory is empty.
Prepare the missing bounded reference from the cached official 1000G Phase 3
chr22 VCF with PLINK2: 503 EUR samples, GRCh37 20,000,000–22,200,000, biallelic
ACGT SNPs, MAF .05 and one thread. qgg computes complete per-gene LD directly
from this BED in VEGAS and fine-mapping; no window-truncated LD is substituted.

Use actual installed GRCh37 Genesplus memberships for genes fully contained
in that interval. Intersect with the same post-QC CAD/T2D marker universe for
both methods. Website fine-mapping QC: MAF .05, frequency difference .05,
INFO .8, ambiguous/indel/duplicate exclusion, missingness .05 and HWE 1e-12.
Record actual counts and gene-size distribution. The website's four named
fine-mapping genes are replaced by the complete available regional gene list;
do not call this a literal full-genome reproduction.

Run at one thread, one successful execution per stage, with fresh R processes:

- VEGAS: both real studies, shared gene sets, public qgg defaults (including
  p-value floor 1e-10 and eigenvalue floor 1e-7).
- Joint BLR-MAGMA: real whole-genome CAD/T2D saved VEGAS statistics and KEGG
  sets from gact::getFeatureSets(minsets=50), BayesC pi=.01, nit=1000,
  nburn=200, matching the website's joint example. Cached VEGAS input preparation
  is excluded; this stage does not consume the separately timed regional VEGAS.
- BayesR mcmc-eigen: both real studies, exactly the VEGAS gene sets,
  eigen thresholds .995/.99/.95/.9, nit=1000, nburn=100, nrun=1, seed20261004.
  Other controls remain installed defaults, including checkConvergence=FALSE,
  as in the website example. Preserve Geweke diagnostics; the returned fit is
  not proof of convergence and default-disabled checking suppresses retries.
- Total and KEGG pathway scores: cached fitted BLR1/T2D effects restricted to
  available regional reference markers, applied to 2,000 existing simulated
  recipients. Remap recipient BIM labels by GRCh37 position and allele pair
  without altering BED dosages. Retain nonempty pathway unions (>=2 markers).
  Public qgg scaling defaults TRUE/TRUE; fitting cost is reused, not measured.

Report analysis wall time separately from loading, preparation and RDS saving.
Peak process memory covers the R process lifetime. PLINK is an auxiliary child
whose peak memory is not included in the R-only preparation memory metric.
Record input/output sizes, exact membership, libraries, controls, warnings,
diagnostics, binary hashes and hardware. No standalone MAGMA comparison is part
of this qgg run; keep the Reviewer 2 proposal separate.
