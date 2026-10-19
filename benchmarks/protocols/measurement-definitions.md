# Measurement definitions

Retain the definitions in each original report. Operation wall time may exclude startup, input loading, preparation and serialization. Whole-process time includes more of those costs; never add overlapping clocks. Separate reference preparation, fitting, score application and reuse.

Memory may be the whole-process Windows peak working set or a sampled analysis-only peak. Neither is solver-only incremental memory, total machine RAM, or summed parent-and-child memory. Missing per-stage peaks remain missing.

Keep real and simulated inputs distinct. Record individuals, variants, genes or regions, sets, traits, iterations, burn-in, chains and threads. State whether measurements are one descriptive run or repeated executions. A thread limit does not imply every operation uses that many threads.

For scoring, report recipient and marker counts, weight source, coding/scaling and marker-set memberships. Existing 20K/50K scores use 2,000 held-out recipients and total plus 20 artificial gene-marker sets; these are not Reactome pathways or UK Biobank timings. Existing timings do not constitute a systematic scaling study.

Import checksums establish fidelity to the local files at import time. They do not identify the exact measured code unless the original run provenance supplies that identification.
