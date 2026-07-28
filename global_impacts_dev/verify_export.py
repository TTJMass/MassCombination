import numpy as np, sys, time
sys.path.insert(0,'.')
from pyconvino.combiner import Combiner
from pyconvino.result import to_dict, export_npz
t0=time.time()
c=Combiner.from_config('pyconvino/ConvinoSetups/CMSOnly/rho_config.txt', prefix='/tmp/sewuchte/vexp', compute_impacts=False)
r=c.combine()
dd=to_dict(r)
print("combine %.0fs; new keys present: x_sys_cov=%s sys_group=%s"%(time.time()-t0, 'x_sys_cov' in dd, 'sys_group' in dd),flush=True)
Sxa=np.asarray(dd['x_sys_cov']); sg=list(dd['sys_group']); names=list(dd['sys_names'])
print("x_sys_cov shape=%s  n sys_group=%d  ungrouped=%d/%d"%(Sxa.shape, len(sg), sg.count('ungrouped'), len(sg)),flush=True)
# closure from EXPORTED arrays only: per-group quadrature == total == w Sxa Sxa^T w
rng=np.random.default_rng(11); w=rng.normal(size=Sxa.shape[0])
dm = w@Sxa                       # per-nuisance global impact (len nsys)
from collections import defaultdict
grp=defaultdict(float)
for i,g in enumerate(sg): grp[g]+=dm[i]**2
syst_from_groups=np.sqrt(sum(grp.values()))
syst_direct=np.sqrt(float(dm@dm))
syst_matrix=np.sqrt(float(w@(Sxa@Sxa.T)@w))
print("closure: sqrt(sum_g q_g)=%.6g  sqrt(sum_i dm^2)=%.6g  sqrt(wSSw)=%.6g"%(syst_from_groups,syst_direct,syst_matrix),flush=True)
print("ratios grp/direct=%.6f  direct/matrix=%.6f"%(syst_from_groups/syst_direct, syst_direct/syst_matrix),flush=True)
# export_npz round-trip to confirm keys land in the .npz file
export_npz(r,'/tmp/sewuchte/vexp_result.npz')
z=np.load('/tmp/sewuchte/vexp_result.npz', allow_pickle=True)
print("in .npz file: x_sys_cov=%s sys_group=%s  x_sys_cov shape=%s"%('x_sys_cov' in z, 'sys_group' in z, z['x_sys_cov'].shape),flush=True)
print("DONE",flush=True)
