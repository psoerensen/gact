from pathlib import Path
import csv,json,hashlib,numpy as np,math
from collections import Counter
import sys
if len(sys.argv)!=3:raise SystemExit('usage: compare_saved_gene_results.py SAVED_RUN_DIRECTORY OUTPUT_DIRECTORY')
base=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
prior=json.loads((base/'standard-gene-agreement.json').read_text());report={'scope':'Descriptive verification from saved results only; no gene or gene-set model rerun','alignment':'Gene IDs, not row positions','gene_tail_methods':{'gsuite':'standard saddlepoint with recorded moment-matched chi-square fallbacks','MAGMA':'standard integration/simulation'},'traits':{},'source_files':{}}
def ranks(x):
 _,inverse,count=np.unique(x,return_inverse=True,return_counts=True);end=np.cumsum(count)
 return ((end-count+end-1)/2)[inverse]
for study in ['GWAS1','GWAS2']:
 suite_file=base/f'{study}-standard-gsuite-genes.csv';magma_file=base/f'standard-{study}.genes.out.txt'
 with suite_file.open() as f: sr=list(csv.DictReader(f))
 lines=[x.split() for x in magma_file.read_text().splitlines() if x.strip() and not x.startswith('#')]
 mr=[dict(zip(lines[0],row)) for row in lines[1:]]
 a={x['GENE']:x for x in mr};b={x['set']:x for x in sr}
 assert len(a)==len(mr)==len(b)==len(sr)==17862 and set(a)==set(b)
 ids=list(a);pa=np.array([float(a[g]['P']) for g in ids]);pb=np.array([float(b[g]['p_value']) for g in ids])
 valid=np.isfinite(pa)&np.isfinite(pb)&(pa>0)&(pb>0)&(pa<=1)&(pb<=1)
 delta=np.abs(pa[valid]-pb[valid]);logdiff=np.abs(-np.log10(pa[valid])+np.log10(pb[valid]))
 gene_cutoff=.05/len(ids);ma=pa<=gene_cutoff;sb=pb<=gene_cutoff
 result={'genes':len(ids),'identical_gene_background':True,'per_gene_snp_count_mismatches':sum(int(a[g]['NSNPS'])!=int(b[g]['marker_count']) for g in ids),
 'effective_parameter_count_mismatches':sum(float(a[g]['NPARAM'])!=float(b[g]['effective_parameters']) for g in ids),
 'valid_p_value_pairs':int(valid.sum()),'invalid_or_unavailable_pairs':int((~valid).sum()),
 'median_absolute_p_difference':float(np.median(delta)),'maximum_absolute_p_difference':float(delta.max()),
 'spearman_negative_log10_p':float(np.corrcoef(ranks(-np.log10(pa[valid])),ranks(-np.log10(pb[valid])))[0,1]),
 'median_absolute_log10_p_difference':float(np.median(logdiff)),'maximum_absolute_log10_p_difference':float(logdiff.max()),
 'descriptive_bonferroni':{'alpha':.05,'denominator':len(ids),'p_threshold':gene_cutoff,'MAGMA_significant':int(ma.sum()),'gsuite_significant':int(sb.sum()),'both_significant':int((ma&sb).sum()),'MAGMA_only':int((ma&~sb).sum()),'gsuite_only':int((sb&~ma).sum()),'same_threshold_decision':int((ma==sb).sum())},
 'gsuite_tail_method_counts':dict(Counter(row['tail_method'] for row in sr)),
 'gsuite_recorded_moment_fallback_count':sum(row['moment_fallback']=='TRUE' for row in sr),
 'qualification':'Strong rank concordance across different standard gene-tail evaluators, not numerical equivalence. Bonferroni counts are a descriptive summary at a stated conventional threshold, not a prespecified equivalence gate.'}
 for old,new in [('median_absolute_p_difference','median_absolute_p_difference'),('maximum_absolute_p_difference','maximum_absolute_p_difference'),('spearman_negative_log10_p','spearman_negative_log10_p')]:
  assert abs(result[new]-prior[study][old])<1e-12,(study,old)
 assert result['per_gene_snp_count_mismatches']==len(prior[study]['marker_count_mismatch'])
 report['traits'][study]=result
 for file in [suite_file,magma_file]:report['source_files'][file.name]={'sha256':hashlib.sha256(file.read_bytes()).hexdigest(),'bytes':file.stat().st_size}
report['existing_summary_verified']=True
(out/'gene-agreement-summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report['traits'],indent=2))
