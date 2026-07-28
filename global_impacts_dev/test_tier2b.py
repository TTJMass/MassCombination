"""Stub test for Tier 2b consumer math (doGlobalImpacts in fit_object.py), ahead of the
pyconvino Tier-2a export existing.

Contract under test (see coordinator's Tier-2b revision, 2026-07-23):
    x_sys_cov : (n_total, nsys) cross-covariance Cov(combined_x, nuisance_i), absolute units.
    sys_group : (nsys,) string label per nuisance ('ungrouped' if none).
    dm_i      = w_abs @ x_sys_cov[:, i]                      (per-nuisance global impact)
    impact_g  = sqrt(sum_{i in group} dm_i^2)                (per-group global impact)
Claim to verify: sqrt(sum_g impact_g^2) == sqrt(sum_i dm_i^2) == sqrt(w_abs . Sxa Sxa^T . w_abs)
exactly (up to floating point), i.e. partitioning nuisances into groups and re-combining in
quadrature recovers the total exactly -- this is pure algebra (Pythagorean sum over a disjoint
partition), so it should hold to ~1e-12 relative, not just "a few %".

x_sys_cov itself does not exist in the NPZ yet; here we stub it as Sxa = cov_full[nsys:, :nsys]
from the saved cmsonly_test.npz (nuisance-nuisance block cov_full[:nsys,:nsys] was already used
for the Tier-1-era hypothesis checks; the cross block cov_full[nsys:,:nsys] is the (n_total, nsys)
piece that plays the role of x_sys_cov here).
"""
import numpy as np

d = np.load('global_impacts_dev/cmsonly_test.npz', allow_pickle=True)
sys_names = [str(s) for s in d['sys_names']]
nsys = len(sys_names)
cov_full = np.asarray(d['cov_full'], dtype=float)
n_total = cov_full.shape[0] - nsys
assert cov_full.shape == (nsys + n_total, nsys + n_total)

Sxa = cov_full[nsys:, :nsys]   # (n_total, nsys) stub for x_sys_cov
assert Sxa.shape == (n_total, nsys)

rng = np.random.default_rng(0)
w_abs = rng.normal(size=n_total)

# Arbitrary group assignment + an 'ungrouped' bucket, mirroring doGlobalImpacts.
sys_group = []
for i, name in enumerate(sys_names):
    if i % 7 == 0:
        sys_group.append('ungrouped')
    else:
        sys_group.append('group{}'.format(i % 5))

dm = w_abs @ Sxa   # (nsys,) per-nuisance dm_i

groups = {}
for name, grp, v in zip(sys_names, sys_group, dm):
    groups.setdefault(grp, []).append(v)
impact_g = {g: float(np.sqrt(np.sum(np.square(vals)))) for g, vals in groups.items()}

n_ungrouped = len(groups.get('ungrouped', []))
print('n_total={}  nsys={}  n_groups={} (incl. ungrouped, {} members)'.format(
    n_total, nsys, len(groups), n_ungrouped))

quad_sum_groups = float(np.sqrt(sum(v ** 2 for v in impact_g.values())))
quad_sum_nuisances = float(np.sqrt(np.sum(dm ** 2)))
direct = float(np.sqrt(w_abs @ Sxa @ Sxa.T @ w_abs))

ratio_group_vs_nuis = quad_sum_groups / quad_sum_nuisances
ratio_nuis_vs_direct = quad_sum_nuisances / direct

print('sqrt(sum_g impact_g^2)   = {:.10f}'.format(quad_sum_groups))
print('sqrt(sum_i dm_i^2)       = {:.10f}'.format(quad_sum_nuisances))
print('sqrt(w_abs Sxa Sxa^T w_abs) = {:.10f}'.format(direct))
print('closure ratio (group/nuisance) = {:.12f}'.format(ratio_group_vs_nuis))
print('closure ratio (nuisance/direct) = {:.12f}'.format(ratio_nuis_vs_direct))

assert abs(ratio_group_vs_nuis - 1.0) < 1e-9
assert abs(ratio_nuis_vs_direct - 1.0) < 1e-9
print('PASS: group partition and direct quadratic form agree to < 1e-9')
