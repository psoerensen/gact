from pathlib import Path
import json,re,hashlib
repo=Path('C:/Users/au223366/Documents/GitHub/gact');archive=repo/'benchmarks/results/20261019-direct-cross-products'
rows=json.loads((archive/'measurements.json').read_text());m={r['stage']:r for r in rows}
def f(v):return f'{v:.2f}'
main=['| Full-data model | Preparation (s) | BLR fit (s) | Whole call (s) | Process peak (GiB) |','| --- | ---: | ---: | ---: | ---: |']
for key in ['magma_single','magma_joint','pops_single_c','pops_single_r']:
 r=m[key];main.append(f"| {r['label']} | {f(r['preparation_seconds'])} | {f(r['blr_fit_seconds'])} | {f(r['operation_seconds'])} | {f(r['peak_working_set_gib'])} |")
validation=['| Validation workflow | Preparation (s) | Sum of BLR fits (s) | Prediction (s) | Whole call (s) | Process peak (GiB) |','| --- | ---: | ---: | ---: | ---: | ---: |']
for key in ['pops_validation_c','pops_validation_r']:
 r=m[key];validation.append(f"| {r['label']} | {f(r['preparation_seconds'])} | {f(r['blr_fit_seconds'])} | {f(r['prediction_seconds'])} | {f(r['operation_seconds'])} | {f(r['peak_working_set_gib'])} |")
def section(link):
 return '''## Current Bayesian MAGMA and single full-data PoPS models

These reruns use **17,955 saved gene scores and 1,324 Reactome sets** from the already installed hsa.0.0.1 database, gsea 0.1.2 and gbayes 0.1.6. **Gene-association statistics were not recomputed**: their upstream cost is excluded, not measured as zero. Cached fixture loading is separate from the calls below; no new database acquisition, genotype processing or simulation was performed.

'''+ '\n'.join(main)+'''

Preparation includes input checks, set mapping and model-specific summaries. Scalar independent-gene set models compute cross-products directly without forming X. The joint-trait model retains its general projected preparation. Native BLR time includes checks, sampling, posterior accumulation and original-unit coefficient back-transformation. Whole-call time also includes posterior summaries, prediction where relevant and remaining interface overhead; serialization and process startup are excluded. Memory is the complete R worker's lifetime peak, including loading and saving, and cannot be assigned to an individual phase.

All rows use two chains, seeds 10/11 and one requested thread. Single-trait MAGMA uses BayesC with pi=.01 and 5,000 sampling plus 1,000 burn-in sweeps; joint MAGMA uses 1,000 plus 200. MAGMA retains its default fixed hyperparameters; PoPS learns slab variance, residual variance and mixture weights. The two single full-data PoPS models use their current PoPS priors and learning flags, pi=.001, 5,000 sampling plus 1,000 burn-in sweeps. For BayesR that scalar pi is divided among three active components; it is not the .95/.02/.02/.01 mixture from the earlier qgg-oriented probe. Pathway residual normalization remains n-1 for an intercept-only fit. These are representative current workflows, not original manuscript model reruns or convergence evidence.

## Separate 22-fold chromosome-held-out validation

'''+ '\n'.join(validation)+'''

These are **22 training models**, each with two chains and 5,000 sampling plus 1,000 burn-in sweeps. They are validation costs and are excluded from the single-fit table. Saved gene scores are reused in every fold; preparation uses only training genes. The native fitting column is the sum over folds, with preparation and held-out prediction recorded separately. Each setting has one execution; there is no repeated-run timing uncertainty estimate.

The updated independent set path needs no dense gene-by-feature matrix. The 1,536 MiB preparation limit is retained for consistency; it is a matrix-payload guard, not measured process memory. The joint model still uses projected preparation. Earlier preparation-included runs and the separate one-chain 14.21-second native probe remain historical evidence, not current rows or matched speedup baselines. A direct sufficient-statistic R interface is now available as `gsea::blr()`, backed by native gbayes. The qgg R/native vb/vg positional mismatch is not reproduced; these timings make no external qgg or MAGMA/FUMA comparison.

See the [rerun report and phase records]('''+link+''') for settings, per-fold timing, binary/input hashes, memory definitions and the excluded duplicate-membership harness diagnostic.

'''
# Preserve existing evidence and unrelated table rows; replace only the obsolete
# current-model and validation blocks.
for rel,start,end,link in [
 ('Document/Benchmark_overview.Rmd','## Current gsuite gene set and single-model BLR measurements','## Controlled simulated scoring grid','../benchmarks/results/20261019-direct-cross-products/README.md'),
 ('benchmarks/revision/consolidated-timings.md','## Current gsuite measurements from 19 October 2026','## Historical supplementary evidence','../results/20261019-direct-cross-products/README.md')]:
 p=repo/rel;s=p.read_text(encoding='utf-8');assert start in s and end in s
 a=s.index(start);b=s.index(end,a);s=s[:a]+section(link)+s[b:];p.write_text(s,encoding='utf-8')
# Keep the tutorial example synchronized, without executing any chunks.
p=repo/'Document/Gene_ranking_bayesian_pops.Rmd';s=p.read_text(encoding='utf-8');a=s.index('## Recorded computational workload')
s=s[:a]+'''## Recorded computational workload

The [benchmark overview](Benchmark_overview.html) records preparation and one full-data BayesC/BayesR fit separately from the 22-fold validation calls above. All reruns use 17,955 saved gene scores and 1,324 Reactome sets. Gene-association statistics are reused, not recomputed in each fold.

With independent genes, no extra covariates and no screening, set-based Bayesian PoPS now computes training cross-products and held-out predictions directly from memberships. Adjusted and screened models retain general preparation. `fit$timing$folds` separates preparation, native BLR, posterior summaries and prediction for each fold. `fit$timing$gene_statistics$seconds` is NA because upstream gene-statistic computation is not measured by this call.

The full-data comparison uses two chains and the same PoPS prior and learning settings as validation, with 1,000 burn-in plus 5,000 sampling sweeps. The [dated rerun report](../benchmarks/results/20261019-direct-cross-products/README.md) contains exact phase clocks and provenance. A direct R interface is also available as `gsea::blr()`; the earlier one-chain qgg-oriented native probe uses different controls and is retained only as historical evidence. Timings do not establish convergence or prediction accuracy.
'''
p.write_text(s,encoding='utf-8')
p=repo/'Document/Gene_set_analysis_bayesian_magma.Rmd';s=p.read_text(encoding='utf-8');a=s.index('## Recorded computational workload')
s=s[:a]+'''## Recorded computational workload

The [benchmark overview](Benchmark_overview.html) records rerun single-trait and joint Bayesian models using 17,955 saved gene scores and 1,324 Reactome sets from the existing hsa.0.0.1 installation. Gene statistics are reused; KEGG extraction shown above is not timed. Scalar independent-gene set preparation now uses direct cross-products without forming X; the joint model retains projected preparation.

`fit$timing` separates preparation, native BLR and posterior summaries. Upstream gene-statistic computation is unmeasured (NA), not zero. The [dated rerun report](../benchmarks/results/20261019-direct-cross-products/README.md) records both two-chain settings and whole-R-process memory. Successful execution is not evidence of convergence or equivalence to the external MAGMA executable.
'''
p.write_text(s,encoding='utf-8')
p=repo/'Document/Tutorial_timings_and_memory.Rmd';s=p.read_text(encoding='utf-8');s=s.replace('current Reactome gene-set measurements, a native BLR fit from prepared summaries and separate chromosome-held-out PoPS validation','updated Reactome Bayesian MAGMA and full-data PoPS models with separate preparation/BLR clocks, and separately timed 22-fold PoPS validation');p.write_text(s,encoding='utf-8')
p=repo/'benchmarks/README.md';s=p.read_text(encoding='utf-8');needle='| Current revision measurements |';at=s.index(needle);s=s.replace('| Current revision measurements |','| Earlier revision measurements |');at=s.index('| Earlier revision measurements |');s=s[:at]+'| Current direct-cross-product reruns | [Rerun report](results/20261019-direct-cross-products/README.md) | Separate preparation, full-data MAGMA/PoPS models and 22-fold PoPS validation |\n'+s[at:];p.write_text(s,encoding='utf-8')
# Earlier evidence directories remain byte-preserved.
p=repo/'benchmarks/revision/README.md';s=p.read_text(encoding='utf-8');s+='\nCurrent Bayesian MAGMA/PoPS table rows come from the [direct-cross-product reruns](../results/20261019-direct-cross-products/README.md). Preparation, single full-data fits and 22-fold validation are separate; upstream gene-statistic computation is reused and excluded. Earlier records remain available as historical evidence.\n';p.write_text(s,encoding='utf-8')
# Refresh hashes for every selected source, preserving the existing inventory.
p=repo/'benchmarks/revision/consolidated-sources.json';data=json.loads(p.read_text(encoding='utf-8'))
for path in ['results/20261019-direct-cross-products/measurements.json','results/20261019-direct-cross-products/measurements.csv','results/20261019-direct-cross-products/README.md']:
 if not any(r['path']==path for r in data['sources']):data['sources'].append({'path':path})
for item in data['sources']:item['sha256']=hashlib.sha256((repo/'benchmarks'/item['path']).read_bytes()).hexdigest()
p.write_text(json.dumps(data,indent=2),encoding='utf-8')
p=repo/'Document/Benchmark_overview.Rmd';s=p.read_text(encoding='utf-8');s=s.replace('- Separate per-phase VEGAS memory remains unavailable; complete-worker memory is already measured.','- Per-phase memory for VEGAS, gene-set preparation and BLR is not separately measured; whole-worker peaks are available for the listed workflows.');p.write_text(s,encoding='utf-8')
print('Updated overview, consolidation and two tutorial workload notes.')
