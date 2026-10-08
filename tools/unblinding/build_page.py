"""Build the click-to-reveal unblinding page from the frozen fit outputs.

usage: build_page.py blind|real frozen|lep075 OUTDIR
  frozen: frozen 2 Oct matrix (lepton ATLAS 8-13 corr 0%)
  lep075: codecheck_20261008/fits_new (lepton corr 0.75), collected in matrix_lep075.csv
  blind: salt = 0, page shows the BLINDED values (safe for visual checks)
  real : salt read from file, page shows the true values. Never prints a value.
"""
import sys, os, json, glob, base64, subprocess, hashlib, html
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

MODE, CONFIG, OUT = sys.argv[1], sys.argv[2], sys.argv[3]; assert CONFIG in ('frozen', 'lep075'); os.makedirs(OUT, exist_ok=True)
H = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.expanduser('~/Claude_MassComb/ReviewTracking/unblinding_2026-10-08')  # matrix_lep075.csv, scanrun*/
SALT_PATH = os.path.expanduser('~/Claude_MassComb/.masscomb_blind_salt')
FP = '4de0c068e866b4fb'
assert hashlib.sha256(open(SALT_PATH, 'rb').read()).hexdigest()[:16] == FP
SALT = 0.0 if MODE == 'blind' else float(open(SALT_PATH).read())
R = os.path.expanduser('~/Claude_MassComb/local_runs_20261002')
MATRIX = R + '/fitmatrix_20261005_expratio/matrix.csv' if CONFIG == 'frozen' else DATA + '/matrix_lep075.csv'
SCAN = DATA + ('/scanrun' if CONFIG == 'frozen' else '/scanrun_lep075') + '/scan_capture.json'
RUN1, RUN1_U = 172.52, 0.33

df = pd.read_csv(MATRIX); df['PDF'] = df['PDF'].astype(str)
L = lambda v: json.loads(v) if isinstance(v, str) and v.startswith('[') else [float(v)]


def row(exp, order='NNLO', pdf='CT18NNLO', poi='single', tnp=True, poly=2, src='stripper_json'):
    s = df[(df.exp == exp) & (df.order == order) & (df.PDF == pdf) & (df.poi_config == poi) & (df.tnp == tnp)
           & (df.poly_order == poly) & (df.expScale == False) & (df.theory_source == src)]
    assert len(s) == 1, (exp, order, pdf, poi, len(s))
    return s.iloc[0]


def unblind(m, blinded, fp):
    """Remove the salt from a blinded mass; single-dataset fits pass unchanged."""
    if blinded:
        assert fp == FP, fp
        return m - SALT
    return m


def masses(r):
    b = str(r['is_blinded']) == 'True'
    return [unblind(m, b, r['blind_salt_fingerprint']) for m in L(r['mass'])], L(r['mass_unc_total'])


def qrun(name):
    j = json.load(open(glob.glob(f'{R}/q_runs/{name}/output/*/results.json')[0]))
    assert j['is_blinded'] and j['blind_salt_fingerprint'] == FP
    return [m - SALT for m in j['mass']], j['mass_unc_total'], j['mass_corr'], j['mass_labels']


f2 = lambda x: f'{x:.2f}'
V = {}  # key -> dict(val, unc)

# main
main = row('ATLAS_813TeV_CMS_13TeV_npz'); (m0,), (s0,) = masses(main)
assert abs(s0 - 0.91) < 0.02
V['main'] = dict(val=f2(m0) + ' ± ' + f2(s0))
d = m0 - RUN1; sd = np.hypot(s0, RUN1_U)
V['run1'] = dict(val=f'{d:+.2f} ± {sd:.2f}', unc=f'{abs(d) / sd:.1f} σ')
for pdf, key in (('NNPDF31_nnlo_as_0118_hessian', 'nnpdf'), ('ABMP16als118_5_nnlo', 'abmp'),
                 ('MSHT20nnlo_as118', 'msht'), ('PDF4LHC21_40', 'pdf4lhc')):
    (m,), (s,) = masses(row('ATLAS_813TeV_CMS_13TeV_npz', pdf=pdf)); V[key] = dict(val=f2(m), unc='± ' + f2(s))
# legacy NLO: total = MINOS (+) externalised scale, decision D2
for pdf, key in (('14400', 'legct18'), ('93300', 'legpdf4lhc')):
    r = row('ATLAS_813TeV_CMS_13TeV_npz', order='NLO', pdf=pdf, tnp=False, src='legacy_root'); (m,), _ = masses(r)
    up = np.hypot(L(r['mass_unc_total_minos_up'])[0], L(r['mass_unc_scale_ext_up'])[0])
    dn = np.hypot(L(r['mass_unc_total_minos_down'])[0], L(r['mass_unc_scale_ext_down'])[0])
    V[key] = dict(val=f2(m), unc=f'+{up:.2f} −{dn:.2f}')
# 2 POI CMS vs ATLAS
r = row('ATLAS_813TeV_CMS_13TeV_npz', poi='split_exp'); ms, ss = masses(r)
lab = json.loads(r['mass_labels'].replace("'", '"')); corr = np.array(json.loads(r['mass_corr']))
ic, ia = lab.index('CMS'), lab.index('ATLAS'); mc, ma, sc, sa, rho = ms[ic], ms[ia], ss[ic], ss[ia], corr[ic, ia]
V['p2cms'] = dict(val=f2(mc), unc='± ' + f2(sc)); V['p2atlas'] = dict(val=f2(ma), unc='± ' + f2(sa))
rat = ma / mc; rs = rat * np.sqrt((sa / ma) ** 2 + (sc / mc) ** 2 - 2 * rho * sa * sc / (ma * mc))
V['p2ratio'] = dict(val=f'{rat:.3f}', unc=f'± {rs:.3f}')
SPLIT_EXP = (ma, mc, sa, sc, rho)
# 8 vs 13 TeV
if CONFIG == 'frozen':
    ms, ss, cr, lab = qrun('nnlo_splitEcm')
else:
    r = row('ATLAS_813TeV_CMS_13TeV_npz', poi='split_ecm'); ms, ss = masses(r)
    cr = json.loads(r['mass_corr']); lab = json.loads(r['mass_labels'].replace("'", '"'))
i8, i13 = lab.index('8TeV'), lab.index('13TeV')
V['e8'] = dict(val=f2(ms[i8]), unc='± ' + f2(ss[i8])); V['e13'] = dict(val=f2(ms[i13]), unc='± ' + f2(ss[i13]))
SPLIT_ECM = (ms[i8], ms[i13], ss[i8], ss[i13], cr[i8][i13])
# 3 POI
r = row('ATLAS_813TeV_CMS_13TeV_npz', poi='split_meas'); ms, ss = masses(r)
lab = json.loads(r['mass_labels'].replace("'", '"'))
for k, key in (('ATLAS_8TeV', 'p3a8'), ('ATLAS_13TeV', 'p3a13'), ('CMS_13TeV', 'p3cms')):
    i = [n.startswith(k) for n in lab].index(True); V[key] = dict(val=f2(ms[i]), unc='± ' + f2(ss[i]))
# pairwise
for exp, key in (('ATLAS_8TeV_CMS_13TeV_npz', 'pw8c'), ('ATLAS_13TeV_CMS_13TeV_npz', 'pw13c'), ('ATLAS_813TeV_npz', 'pwa')):
    (m,), (s,) = masses(row(exp)); V[key] = dict(val=f2(m), unc='± ' + f2(s))
# single datasets (never blinded)
SINGLE = []
for exp, lbl in (('ATLAS_8TeV_npz', 'ATLAS 8 TeV'), ('ATLAS_13TeV_npz', 'ATLAS 13 TeV'), ('CMS_13TeV_npz', 'CMS 13 TeV')):
    r = row(exp); assert str(r['is_blinded']) == 'False'
    (m,), (s,) = masses(r); SINGLE.append((lbl, f2(m) + ' ± ' + f2(s)))

# ---------------- plots
sys.path.insert(0, os.path.expanduser('~/Claude_MassComb/ReviewTracking/review_2026-10-04_prl/r2'))
import prlstyle as S
MT = r'$m_t^{\mathrm{pole}}$'
# chi2 scan (points from the blinded rerun of the frozen nominal fit, x shifted by the salt)
cap = json.load(open(SCAN))
sc_ = [c for c in cap if 'scan' in str(c['label'])][0]
x, y = np.array(sc_['x']) - SALT, np.array(sc_['y'])
near = np.abs(x - x[np.argmin(y)]) < 1.0
c = np.polyfit(x[near], y[near], 2); xm = -c[1] / (2 * c[0]); ymin = np.polyval(c, xm); sig = 1 / np.sqrt(c[0])
assert abs(xm - m0) < 0.02, 'scan minimum does not close with the frozen fit'
S.use(); fig, ax = plt.subplots(figsize=(S.COL_W, 2.5))
ax.plot(x, y - ymin, 'o', ms=3, color=S.OI['black'], label=r'$\chi^2$ scan')
xx = np.linspace(xm - 1.1 * sig, xm + 1.1 * sig, 100)
ax.plot(xx, np.polyval(c, xx) - ymin, color=S.OI['vermillion'], lw=1.2, label='Parabola near the minimum')
ax.axhline(1, color='0.5', ls='--', lw=0.7)
for s_ in (xm - sig, xm + sig): ax.axvline(s_, color='0.5', ls=':', lw=0.7)
ax.set_xlim(xm - 2, xm + 2); ax.set_ylim(0, 6.2)
ax.set_xlabel(MT + ' [GeV]', loc='right'); ax.set_ylabel(r'$\Delta\chi^2$', loc='top')
ax.text(0.97, 0.95, rf'{MT} = {xm:.2f} $\pm$ {sig:.2f} GeV' + '\n' + rf'$\chi^2_{{\min}}/\mathrm{{ndf}}$ = {ymin:.1f}/15',
        transform=ax.transAxes, ha='right', va='top', fontsize=6.5)
ax.legend(loc='upper left', frameon=False); S.label(ax)
fig.savefig(OUT + '/chi2_scan.png'); plt.close(fig)


def contour(p, nx, ny, name):
    px, py, sx, sy, rho = p
    cv = np.array([[sx * sx, rho * sx * sy], [rho * sx * sy, sy * sy]]); w, Vv = np.linalg.eigh(cv)
    ang = np.degrees(np.arctan2(Vv[1, 1], Vv[0, 1]))
    S.use(); fig, ax = plt.subplots(figsize=(S.COL_W, S.COL_W * 0.92))
    for k, a in ((5.99, 0.18), (2.30, 0.45)):
        ax.add_patch(Ellipse((px, py), 2 * np.sqrt(k * w[1]), 2 * np.sqrt(k * w[0]), angle=ang, color=S.OI['blue'], alpha=a, lw=0))
    lo, hi = m0 - 4.5, m0 + 4.5
    ax.plot(px, py, 'o', ms=4, color='k'); ax.plot([lo, hi], [lo, hi], 'k--', lw=0.7); ax.plot(m0, m0, '*', ms=9, color=S.OI['orange'])
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect('equal')
    ax.set_xlabel(nx + ' [GeV]', loc='right'); ax.set_ylabel(ny + ' [GeV]', loc='top')
    ax.text(0.97, 0.05, rf'$\rho$ = {100 * rho:.0f}%', transform=ax.transAxes, ha='right', fontsize=7)
    ax.legend(handles=[plt.Line2D([], [], color='k', marker='o', ls='', ms=4, label='Best fit'),
                       plt.Line2D([], [], color=S.OI['blue'], alpha=.45, lw=6, label='68% CL'),
                       plt.Line2D([], [], color=S.OI['blue'], alpha=.18, lw=6, label='95% CL'),
                       plt.Line2D([], [], color=S.OI['orange'], marker='*', ls='', ms=8, label='One-parameter fit'),
                       plt.Line2D([], [], color='k', ls='--', lw=0.7, label='Equal masses')], loc='upper left', frameon=False)
    S.label(ax); fig.savefig(f'{OUT}/{name}.png'); plt.close(fig)


contour(SPLIT_EXP, r'$m_t^{\mathrm{ATLAS}}$', r'$m_t^{\mathrm{CMS}}$', 'twopoi_exp')
contour(SPLIT_ECM, r'$m_t^{8\,\mathrm{TeV}}$', r'$m_t^{13\,\mathrm{TeV}}$', 'twopoi_ecm')
# summary plot (paper Fig. 2 layout): prl_summary_options4.py removes the salt itself (dry run in blind mode)
subprocess.run([sys.executable, 'prl_summary_options4.py', MATRIX, OUT],
               cwd=os.path.expanduser('~/Claude_MassComb/MassCombination/MtopSummaryPlot'), check=True,
               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
               env=dict(os.environ, MASSCOMB_UNBLIND_DRYRUN='1' if MODE == 'blind' else '0'))
img = lambda n: 'data:image/png;base64,' + base64.b64encode(open(f'{OUT}/{n}.png', 'rb').read()).decode()

# ---------------- page
tpl = open(H + '/template.html').read()
for k, v in V.items():
    tpl = tpl.replace('{{%s}}' % k, html.escape(v['val'])).replace('{{%s_u}}' % k, html.escape(v.get('unc', '')))
tpl = tpl.replace('{{source}}', html.escape(MATRIX.replace(os.path.expanduser('~/Claude_MassComb/'), '')) + (' (lepton ATLAS 8–13 correlation 0.75)' if CONFIG == 'lep075' else ' (frozen, lepton correlation 0%)'))
tpl = tpl.replace('{{single}}', ''.join(f'<tr><td>{a}</td><td>{b} GeV</td></tr>' for a, b in SINGLE))
for n in ('prl_summary_round', 'chi2_scan', 'twopoi_exp', 'twopoi_ecm'):
    tpl = tpl.replace('{{img_%s}}' % n, img(n))
tpl = tpl.replace('{{banner}}', '' if MODE == 'real' else
                  '<div class="mock">BLINDED CHECK COPY. Combined values still include the salt.</div>')
assert '{{' not in tpl
open(f'{OUT}/unblinding.html', 'w').write(tpl)
print('wrote', f'{OUT}/unblinding.html', '(mode', MODE + ')')
