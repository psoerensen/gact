"""Compare installed outputs to pinned original HESS variance functions.

Usage: python compare_hess_software_reference.py PINNED_ESTIMATION_PY MOMENTS_CSV OUTPUT_DIR
The source is caller-supplied, never a runtime/package dependency. Only its
second-stage function definitions are executed. Projection/load adapters use
saved moments and native point estimates; no LD/eigenvectors are recomputed.
"""
import ast
import hashlib
import logging
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import scipy.stats

source, moments_path, output = map(Path,sys.argv[1:])
text=source.read_text()
expected='d53aa10941746d7af8c164404aabf4bd3ea8e9a614078fb1e5d6bafd6b9cb5ee'
assert hashlib.sha256(text.encode()).hexdigest()==expected, 'Wrong pinned source'
np.int=int
np.float=float
tree=ast.parse(text)
tree.body=[node for node in tree.body if isinstance(node,ast.FunctionDef)
  and node.name in {'local_hsqg_var','local_rhog_step2','estimate_gw_rg'}]
namespace=dict(np=np,pd=pd,scipy=scipy,logging=logging,sys=sys,xrange=range,eps=1e-8)
exec(compile(tree,str(source),'exec'),namespace)
m=pd.read_csv(moments_path)
info=[pd.DataFrame(dict(N=m[f'n{j}'],rank=m.numerical_rank,nsnp=m.markers,
  CHR=1,start=np.arange(len(m)),stop=np.arange(len(m))+1)) for j in (1,2)]
rows=[]
for scenario,factors in [('raw',[1.,1.]),('partial_gc',[1.,1.18])]:
  native=pd.read_csv(output/f'{scenario}-estimates.csv')
  region=native[native.component=='region']
  fits=[]
  for j in (0,1):
    values=region[(region.quantity=='h2') & (region.trait1==f'GWAS{j+1}')]
    assert list(values.region)==list(m.region)
    h=pd.Series(values.estimate.to_numpy())
    v=namespace['local_hsqg_var'](info[j],h).to_numpy()
    current=values.standard_error.to_numpy()**2
    good=(v>=0) & np.isfinite(v)
    assert np.all(np.isnan(current[~good]))
    difference=float(np.max(np.abs(current[good]-v[good])))
    assert difference<1e-12
    rows.append(dict(scenario=scenario,quantity='h2',trait=f'GWAS{j+1}',
      max_variance_difference=difference,negative_variances=int((v<0).sum())))
    fits.append(pd.DataFrame({'local_h2g':h,'var':v}))
  namespace['load_local_hsqg_est']=lambda _:fits
  namespace['load_local_rhog_step1']=lambda _:(info[0],info[1],None,None)
  namespace['local_bilin_form']=lambda eig,prjprod,max_k,min_eigval,reinflate:pd.DataFrame(
    {'sum':m.q12*np.sqrt(np.prod(reinflate)),'k':m.k})
  saved=[]
  original_csv=pd.DataFrame.to_csv
  pd.DataFrame.to_csv=lambda frame,*a,**k:saved.append(frame.copy())
  try:
    namespace['local_rhog_step2'](None,None,11,1,factors,None,0.,0.,'unused')
  finally:
    pd.DataFrame.to_csv=original_csv
  values=region[region.quantity=='cov']
  v=saved[0]['var'].to_numpy()
  difference=float(np.max(np.abs(values.standard_error.to_numpy()**2-v)))
  assert difference<1e-12
  rows.append(dict(scenario=scenario,quantity='cov',trait='GWAS1/GWAS2',
    max_variance_difference=difference,negative_variances=int((v<0).sum())))
pd.DataFrame(rows).to_csv(output/'upstream-variance-parity.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
