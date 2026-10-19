from pathlib import Path
import json,csv,hashlib,shutil,datetime
base=Path.cwd();repo=Path('C:/Users/au223366/Documents/GitHub/gact');run=base/'work/direct-workflow-rerun-20261019'
archive=repo/'benchmarks/results/20261019-direct-cross-products'
archive.mkdir(parents=True,exist_ok=True)
rows=[]
labels={'magma_single':'Bayesian MAGMA: single-trait BayesC','magma_joint':'Bayesian MAGMA: joint GWAS1/GWAS2 BayesC','pops_single_c':'PoPS: one full-data BayesC model','pops_single_r':'PoPS: one full-data BayesR model','pops_validation_c':'PoPS BayesC: 22 held-out chromosome folds','pops_validation_r':'PoPS BayesR: 22 held-out chromosome folds'}
for stage,label in labels.items():
 d=run/stage
 timing=json.loads((d/'timing.json').read_text(encoding='utf-8-sig'));phases=json.loads((d/'stages.json').read_text(encoding='utf-8-sig'));process=json.loads((d/'process.json').read_text(encoding='utf-8-sig'))
 assert timing['status']=='completed' and process['exit_code']==0
 r=dict(stage=stage,label=label,scope='22-fold validation' if 'validation' in stage else 'one full-data model',
  preparation_seconds=phases['preparation_seconds'],blr_fit_seconds=phases['blr_fit_seconds'],
  posterior_summary_seconds=phases['posterior_summary_seconds'],prediction_seconds=phases.get('prediction_seconds',0),
  operation_seconds=timing['operation_seconds'],input_loading_seconds=timing['input_loading_seconds'],
  peak_working_set_gib=process['peak_working_set_bytes']/1024**3,chains=2,threads=1,
  sampling_sweeps=1000 if stage=='magma_joint' else 5000,burnin=200 if stage=='magma_joint' else 1000,
  gene_statistics_seconds=None)
 rows.append(r)
 dest=archive/stage;dest.mkdir(exist_ok=True)
 for name in ['timing.json','stages.json','process.json','settings.json','effective-settings.json','session.txt','fold-timings.csv','fold-effective-settings.json','stdout.log','stderr.log']:
  if (d/name).exists():shutil.copy2(d/name,dest/name)
with (archive/'measurements.csv').open('w',encoding='utf-8',newline='')as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(archive/'measurements.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
scripts=archive/'scripts';scripts.mkdir(exist_ok=True)
for name in ['direct-workflow-benchmark-worker.R','direct-workflow-benchmark-worker-v2.R','run-direct-workflow-benchmarks.ps1','run-direct-validation-benchmarks.ps1']:
 shutil.copy2(base/'work'/name,scripts/name)
for name in ['diagnostic-pops-single-c-duplicate-memberships']:
 if (run/name).exists():
  dest=archive/'diagnostics'/name;dest.mkdir(parents=True,exist_ok=True)
  for f in (run/name).iterdir():
   if f.suffix in ['.json','.txt','.log']:shutil.copy2(f,dest/f.name)
raw=repo/'benchmarks/work/direct-workflow-rerun-20261019';raw.mkdir(parents=True,exist_ok=True)
retained=[]
for stage in labels:
 dest=raw/stage;dest.mkdir(exist_ok=True)
 for name in ['output.rds','predictions.csv']:
  src=run/stage/name
  if src.exists():shutil.copy2(src,dest/name);retained.append(dict(original_path=str(src),local_ignored_copy=str(dest/name),bytes=src.stat().st_size,sha256=hashlib.sha256(src.read_bytes()).hexdigest()))
(archive/'retained-artifacts.json').write_text(json.dumps(retained,indent=2),encoding='utf-8')
sources=[base/'work/revision-benchmark/input-preparation/output.rds']
for rel in ['packages/gsea/R/pops.R','packages/gsea/R/pops-bayesian.R','packages/gsea/R/pathways.R','packages/gsea/R/qgg-interface.R','packages/gsea/src/pathways.cpp','build/task-packages/prefix/bin/libgsea.dll','build/task-packages/prefix/bin/libgsea_bayes.dll','build/task-packages/prefix/bin/libgbayes.dll','build/task-packages/library/gsea/libs/x64/gsea.dll']:
 sources.append(Path('C:/Users/au223366/Documents/GitHub/gsuite')/rel)
for rel in ['src/set_cross_products.cpp','src/bayesian_pathways.cpp','include/gsea/bayesian_pathways.hpp']:
 sources.append(Path('C:/Users/au223366/Documents/GitHub/gsea')/rel)
manifest=[dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())for p in sources]
(archive/'provenance.json').write_text(json.dumps(dict(input='Frozen local installed hsa.0.0.1 saved gene statistics and Reactome sets',no_gene_testing=True,loading_separate=True,single_runs=True,files=manifest),indent=2),encoding='utf-8')
shutil.copy2(run/'validation-checks.json',archive/'validation-checks.json')
print('Archived',len(rows),'successful stages to',archive)
