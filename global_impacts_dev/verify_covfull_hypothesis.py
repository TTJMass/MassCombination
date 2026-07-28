import numpy as np, time, sys
sys.path.insert(0,'.')
from pyconvino.combiner import Combiner
from pyconvino.result import to_dict
cfg='pyconvino/ConvinoSetups/CMSOnly/rho_config.txt'
t0=time.time()
c = Combiner.from_config(cfg, prefix='/tmp/sewuchte/covtest_cmsonly', compute_impacts=True)
r = c.combine()
print("combine time %.1fs converged=%s" % (time.time()-t0, r.converged), flush=True)
cov_full=np.asarray(r.cov_full,float); nsys=len(r.sys_names); nest=len(r.combined_names)
print("nsys=%d nest=%d cov_full=%s"%(nsys,nest,cov_full.shape), flush=True)
Sigma_post=cov_full[:nsys,:nsys]
A=np.asarray(r.impact_matrix,float)
dd=to_dict(r)
Vsyst=np.asarray(dd['total_syst_covariance'],float)
constraints=np.asarray(r.constraints,float)
den=np.max(np.abs(Vsyst))
recon=A@Sigma_post@A.T
print("\nHYP1 Vsyst = A Sigma_post A^T : max rel err = %.3e  trace ratio=%.6f"%(np.max(np.abs(recon-Vsyst))/den, np.trace(recon)/np.trace(Vsyst)), flush=True)
recon_d=(A*constraints[None,:]**2)@A.T
print("HYP2 A diag(c^2) A^T (diag approx): max rel err = %.3e"%(np.max(np.abs(recon_d-Vsyst))/den), flush=True)
off=Sigma_post-np.diag(np.diag(Sigma_post))
print("Sigma_post offdiag/diag = %.3f ; diag vs c^2 rel = %.2e"%(np.max(np.abs(off))/np.max(np.abs(np.diag(Sigma_post))), np.max(np.abs(np.diag(Sigma_post)-constraints**2))/np.max(constraints**2)), flush=True)
np.savez('/tmp/sewuchte/cmsonly_test.npz', cov_full=cov_full, impact_matrix=A, total_syst_covariance=Vsyst, stat_only_covariance=np.asarray(dd['stat_only_covariance'],float), combined_covariance=np.asarray(dd['combined_covariance'],float), constraints=constraints, sys_names=np.array(r.sys_names))
print("SAVED /tmp/sewuchte/cmsonly_test.npz", flush=True)
