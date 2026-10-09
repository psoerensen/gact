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
