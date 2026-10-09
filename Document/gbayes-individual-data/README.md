# Individual-level gbayes teaching run

Recorded 9 October 2026. See the parent tutorial source for the complete recipe.
Eight fits compare BayesC and BayesR, full and scheduled updates, for four- and
500-causal-marker traits. Two chains per fit use 200 burn-in and 500 sampling
sweeps, thinning effects by five. Training-only QC retains 37,994 markers.

These are short teaching runs with incomplete mixing, not convergence
qualification or evidence for ranking models. Scheduling is approximate;
polygenic BayesR visited every marker with these settings. The independent
simulated panel does not validate realistic population-LD behavior.

The summary tables and figures contain no individual predictions or genotypes.
Input SHA-256 values identify the reused qgdata files and corrected gact BIM.
The original panel generation seed and genome build are unrecorded.
Variance identities and all-post-burn-in means were checked for all eight fits.
Focused BED and continuation adapter tests passed; no full R CMD check is claimed.

Inputs: https://github.com/psoerensen/qgdata/tree/main/simulated_human_data
Corrected map: https://psoerensen.github.io/gact/Document/glma-simulated-data/

The simple-posterior table was recalculated from the original saved chain
draws, using every post-burn-in BED variance sweep and type-8 equal-tail
quantiles. Existing VB/VE rows retain their native effect-thinned summaries.
The original runs predate the gbayes 0.1.3 posterior-table adapter update.

The mixing investigation preserves the original priors and prediction results.
Four original fits were continued for two 500-sweep batches. A controlled
polygenic scheduled BayesR check starts from the same original saved states
using native correction c198232 and R package 0.1.4. Complete scheduled marker
visits now learn mixture weights on every sweep: 500 updates instead of 50
per batch. Prediction accuracy remains similar, but the polygenic posterior
is still not qualified. Parameter groups and per-chain update counts are
reported separately; the screening thresholds are not convergence proof.
See mixing-investigation-record.json for source revisions, DLL hashes and
focused validation scope. Original eight-fit evidence remains unchanged.
