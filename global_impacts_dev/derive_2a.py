import numpy as np, time, sys
sys.path.insert(0,'.')
from pyconvino.combiner import Combiner
cfg='pyconvino/ConvinoSetups/CMSOnly/rho_config.txt'
t0=time.time()
c=Combiner.from_config(cfg, prefix='/tmp/sewuchte/derive2a', compute_impacts=False)
r=c.combine()
print("combine(no-impacts) %.1fs converged=%s"%(time.time()-t0,r.converged),flush=True)
sys_names=list(r.sys_names); nsys=len(sys_names); nest=len(r.combined_names)
cov_full=np.asarray(r.cov_full,float); A=np.asarray(r.impact_matrix,float)
Sig=cov_full[:nsys,:nsys]
Vsyst = A@Sig@A.T                       # verified = freeze total_syst to 4.6%
# group membership from config
groups=c.config.impact_groups           # dict label -> members
name2idx={n:i for i,n in enumerate(sys_names)}
print("n groups=%d nsys=%d nest=%d"%(len(groups),nsys,nest),flush=True)
# per-group within-block contribution D_g = A_g Sig[g,g] A_g^T
den=np.max(np.abs(Vsyst))
Dsum=np.zeros_like(Vsyst); covered=set()
Dg_list=[]
for label,members in groups.items():
    idx=[name2idx[m] for m in members if m in name2idx]
    if not idx: continue
    covered.update(idx)
    Ag=A[:,idx]; Sgg=Sig[np.ix_(idx,idx)]
    Dg=Ag@Sgg@Ag.T
    Dg_list.append((label,Dg,idx))
    Dsum+=Dg
missing=[i for i in range(nsys) if i not in covered]
print("nuisances not in any group: %d"%len(missing),flush=True)
print("max|sum_g D_g - Vsyst|/max|Vsyst| = %.3e  (cross-group residual)"%(np.max(np.abs(Dsum-Vsyst))/den),flush=True)
print("trace(sum_g D_g)/trace(Vsyst) = %.4f"%(np.trace(Dsum)/np.trace(Vsyst)),flush=True)
# mass-level test with random w (matrix fact, w-independent conclusions)
rng=np.random.default_rng(2); 
for trial in range(3):
    w=rng.normal(size=nest)
    tot=float(w@Vsyst@w)
    per=np.array([float(w@Dg@w) for _,Dg,_ in Dg_list])
    print("  w-trial: sum_g (wDgw)=%.4g  wVsw=%.4g  ratio=%.4f  min_group=%.3g"%(per.sum(),tot,per.sum()/tot,per.min()),flush=True)
np.savez('/tmp/sewuchte/derive2a.npz', cov_full=cov_full, impact_matrix=A, sys_names=np.array(sys_names),
         group_labels=np.array([l for l,_,_ in Dg_list], dtype=object),
         Vsyst=Vsyst)
print("SAVED /tmp/sewuchte/derive2a.npz",flush=True)
