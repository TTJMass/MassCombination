#!/usr/bin/env python3
"""Toy check of the linearised (J C J^T) covariance of the normalised combined cross sections."""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import uncertainties as unc
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'mtpole-ttj-pyconvino'))
from exp_xsec import exp_xsec_npz

p = argparse.ArgumentParser()
p.add_argument('npz')
p.add_argument('--out', default='normToyCheck')
p.add_argument('--ntoys', type=int, default=1_000_000)
p.add_argument('--seed', type=int, default=1)
a = p.parse_args()
os.makedirs(a.out, exist_ok=True)

e = exp_xsec_npz(a.npz)
x, C = e.xsec.astype(float), e.cov_matrix
yhat = np.array([v.n for v in e.u_xsec])
Clin = np.array(unc.covariance_matrix(e.u_xsec))
slices = {m: slice(e.index_meas[m][0], e.index_meas[m][1] + 1) for m in e.input_meas}


def normalise(X):
    # exact (non-linearised) normalisation, same convention as exp_xsec.normaliseXsec
    Y = X.copy()
    for m, s in slices.items():
        b = np.diff(e.binning[m])
        if m.split('_')[4] == 'abs':
            Y[:, s] /= b
        Y[:, s] /= (Y[:, s] @ b)[:, None]
    return Y


Y = normalise(np.random.default_rng(a.seed).multivariate_normal(x, C, size=a.ntoys))
Ctoy = np.cov(Y.T)
sig_lin, sig_toy = np.sqrt(np.diag(Clin)), np.sqrt(np.diag(Ctoy))
corr_lin, corr_toy = Clin / np.outer(sig_lin, sig_lin), Ctoy / np.outer(sig_toy, sig_toy)
n = len(x)
labels = {m: m.split('_')[0] + ' ' + m.split('_')[1] for m in e.input_meas}
edges = [s.start - 0.5 for s in slices.values()][1:]

# 1) width ratio and bias per bin
fig, ax = plt.subplots(2, 1, figsize=(7, 5), sharex=True)
eps = 1 / np.sqrt(2 * a.ntoys)
ax[0].errorbar(range(n), sig_toy / sig_lin - 1, yerr=eps, fmt='o', color='k', ms=4)
ax[0].set_ylabel(r'$\sigma_{\rm toy}/\sigma_{\rm lin} - 1$')
ax[0].margins(y=0.25)
ax[1].errorbar(range(n), (Y.mean(0) - yhat) / sig_lin, yerr=1 / np.sqrt(a.ntoys), fmt='o', color='k', ms=4)
ax[1].set_ylabel(r'$(\langle y\rangle_{\rm toy} - \hat y)/\sigma_{\rm lin}$')
ax[1].set_xlabel('bin index')
for axi in ax:
    axi.axhline(0, color='grey', lw=1)
    for ed in edges:
        axi.axvline(ed, color='grey', ls=':', lw=1)
for m, s in slices.items():
    ax[0].text((s.start + s.stop - 1) / 2, 0.97, labels[m], transform=ax[0].get_xaxis_transform(),
               ha='center', va='top', fontsize=8)
fig.tight_layout()
fig.savefig(f'{a.out}/width_bias.pdf'); fig.savefig(f'{a.out}/width_bias.png', dpi=150)

# 2) correlation matrices
fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
dmax = np.max(np.abs(corr_toy - corr_lin))
for axi, M, t, vm in [(ax[0], corr_lin, r'linear $J C J^T$', 1), (ax[1], corr_toy, f'toys ({a.ntoys:.0e})', 1),
                      (ax[2], corr_toy - corr_lin, 'toys - linear', dmax)]:
    im = axi.imshow(M, cmap='RdBu_r', vmin=-vm, vmax=vm)
    axi.set_title(t)
    for ed in edges:
        axi.axvline(ed, color='k', lw=0.8); axi.axhline(ed, color='k', lw=0.8)
    fig.colorbar(im, ax=axi, fraction=0.046)
fig.tight_layout()
fig.savefig(f'{a.out}/correlations.pdf'); fig.savefig(f'{a.out}/correlations.png', dpi=150)

# 3) chi2 of each toy vs the fit convention (first bin per measurement dropped), and bin-drop invariance
def chi2(drop):
    keep = np.setdiff1d(np.arange(n), drop)
    R = (Y - yhat)[:, keep]
    return np.einsum('ij,jk,ik->i', R, np.linalg.inv(Clin[np.ix_(keep, keep)]), R)

ndf = n - len(slices)
c_first = chi2([s.start for s in slices.values()])
c_last = chi2([s.stop - 1 for s in slices.values()])
fig, ax = plt.subplots(figsize=(6, 4))
bins = np.linspace(0, 50, 101)
ax.hist(c_first, bins=bins, density=True, histtype='step', color='k', label='toys')
ax.plot(bins, stats.chi2.pdf(bins, ndf), color='C3', label=rf'$\chi^2$ pdf, ndf = {ndf}')
ax.set_xlabel(r'$\chi^2$ (first bin per measurement dropped)')
ax.set_ylabel('density')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(f'{a.out}/chi2.pdf'); fig.savefig(f'{a.out}/chi2.png', dpi=150)

pvals = stats.chi2.sf(c_first, ndf)
print(f'toys: {a.ntoys}   MC precision on sigma ratio: {eps:.1e}')
print(f'max |sigma_toy/sigma_lin - 1| = {np.max(np.abs(sig_toy / sig_lin - 1)):.4f}')
print(f'max |corr_toy - corr_lin|     = {dmax:.4f}')
print(f'max |<y>_toy - yhat|/sigma    = {np.max(np.abs(Y.mean(0) - yhat) / sig_lin):.4f}')
print(f'chi2: mean {c_first.mean():.3f} (expect {ndf}), var {c_first.var():.3f} (expect {2 * ndf}); '
      f'frac p<0.05 = {np.mean(pvals < 0.05):.4f}, p<0.01 = {np.mean(pvals < 0.01):.4f}')
print(f'bin-drop invariance: max |chi2(first dropped) - chi2(last dropped)| = {np.max(np.abs(c_first - c_last)):.2e}')
