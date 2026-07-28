# Uncertainty breakdown of the combined top-quark mass: three impact definitions

*Math note for the two-stage combination (`pyconvino` cross-section combination →
`mtpole-ttj-pyconvino` χ² mass fit). Covers (A) the classic freeze/conditional
impact, (B) the marginal "global impact" of arXiv:2307.04007, and (C) the exact
additive variance-share generalisation used for correlated priors.*

All numbers quoted are **uncertainties only** (the fit is blinded); they come from the
3-experiment `Combination_ATLAS813CMS13_corrV2` fit (ATLAS 8 TeV + ATLAS 13 TeV +
CMS 13 TeV, NNLO / CT18NNLO, `polyOrder 1`).

---

## 0. Setup and notation

**Stage 1 — combination.** `pyconvino` combines the input differential measurements
into a vector of combined cross sections $x\in\mathbb{R}^{n}$ (here $n=19$ bins),
profiling $N$ experimental nuisance parameters $\theta\in\mathbb{R}^{N}$ (here
$N=614$). Each nuisance carries a Gaussian prior

$$\theta \sim \mathcal{N}\!\big(0,\ C\big),\qquad C^{-1}\equiv \texttt{inv\_C}\ (\text{prior precision}).$$

For the usual "unit, independent" priors $C=\mathbb{1}$; the `corrV2` setup deliberately
makes $C$ (hence $\texttt{inv\_C}$) **non-diagonal**, correlating common systematics
*across experiments*. The combination returns the joint post-fit covariance

$$\text{cov\_full}=
\begin{pmatrix}\Sigma_{xx} & \Sigma_{x\theta}\\[2pt] \Sigma_{\theta x} & \Sigma_{\theta\theta}\end{pmatrix},
\qquad
\Sigma_{xx}\equiv C_{\rm tot}=\texttt{combined\_covariance},\quad
\Sigma_{x\theta}=\texttt{x\_sys\_cov}.$$

**Stage 2 — mass fit.** In the Gaussian/linear regime the mass estimator is a linear
functional of the combined cross sections,

$$\boxed{\ \hat m - m_0 = w^\top\,(x-x_0)\ }\qquad
w=\text{Fisher weights (\texttt{\_fisher\_weights})}.$$

(The code applies a normalisation Jacobian $J$ so that $w_{\rm abs}=J^\top w$ acts in
absolute-cross-section space; that is an implementation detail and is suppressed below.)

**The one vector that drives everything.** Define the covariance of the mass with each
nuisance,

$$\boxed{\ g_i \equiv \mathrm{Cov}(\hat m,\theta_i) = \big(w^\top \Sigma_{x\theta}\big)_i\ }
\qquad (\text{code: } \texttt{dm\_i} = w_{\rm abs}\!\cdot\!\texttt{x\_sys\_cov}[:,i]).$$

Also write $\sigma_{{\rm post},i}=\sqrt{\Sigma_{\theta\theta,ii}}$ (post-fit width of
$\theta_i$) and $\sigma_{{\rm pre},i}=\sqrt{C_{ii}}$ (prior width). The total experimental
mass variance is $\sigma_m^2 = w^\top C_{\rm tot}\,w$.

Everything below is a different way of attributing $\sigma_m^2$ (or its systematic part)
to the individual $\theta_i$.

---

## A. Classic impact — freeze / conditional method

**Recipe.** Fix one nuisance to its best fit $\pm$ its **post-fit** uncertainty,
$\theta_i=\hat\theta_i\pm\sigma_{{\rm post},i}$, re-minimise, read the shift in $\hat m$.

**Closed form.** For a Gaussian joint posterior of $(\hat m,\theta)$ the conditional mean
gives

$$\boxed{\ I^{\rm freeze}_i=\frac{\mathrm{Cov}(\hat m,\theta_i)}{\sigma_{{\rm post},i}}
=\frac{g_i}{\sigma_{{\rm post},i}}=\sigma_m\,\rho(\hat m,\theta_i)\ }$$

i.e. the familiar "impact $=\sigma_m\times$ correlation".

**Stat/syst split (Tier-1).** The systematic is defined by *freezing all* nuisances:

$$\sigma_{\rm stat}^{\rm freeze}=\big(\text{mass unc. with all }\theta\text{ fixed}\big),\qquad
\big(\sigma_{\rm syst}^{\rm freeze}\big)^2=\sigma_m^2-\big(\sigma_{\rm stat}^{\rm freeze}\big)^2
= w^\top V_{\rm syst}\,w,$$

with $V_{\rm syst}=C_{\rm tot}-V_{\rm stat}$ (code: `total_syst_covariance`). Equivalently
$f=\big(w^\top V_{\rm syst} w\big)/\big(w^\top C_{\rm tot} w\big)$ and
$\sigma_{\rm syst}=\sigma_m\sqrt{f}$, $\sigma_{\rm stat}=\sigma_m\sqrt{1-f}$.

**Key property / caveat.** The freeze *stat/syst split closes* ($\sigma_{\rm stat}\oplus
\sigma_{\rm syst}=\sigma_m$), **but the per-nuisance $I^{\rm freeze}_i$ do NOT add up** to
$\sigma_{\rm syst}^{\rm freeze}$ — the classic "impacts don't sum" fact, because the
$\theta_i$ are correlated in the post-fit and each $I^{\rm freeze}_i$ double-counts shared
constraint.

*Live numbers:* $\sigma_{\rm stat}^{\rm freeze}=424.7$ MeV, $\sigma_{\rm syst}^{\rm freeze}=733.1$ MeV
($\oplus=847.2$ MeV experimental total).

---

## B. Global impact — marginal (arXiv:2307.04007)

**Recipe.** Shift the *external/auxiliary observable* $y_i$ that constrains $\theta_i$ by
its **prefit** width $\sigma_{{\rm pre},i}$, let all $\theta$ **float**, refit, read $\Delta\hat m$.

**Relation to the freeze impact.** The two differ by the *post-fit nuisance-parameter
covariance* factor $\sigma_{{\rm post},i}/\sigma_{{\rm pre},i}$:

$$\boxed{\ I^{\rm glob}_i = I^{\rm freeze}_i\times\frac{\sigma_{{\rm post},i}}{\sigma_{{\rm pre},i}}
=\frac{g_i}{\sigma_{{\rm pre},i}}\ }\ \xrightarrow{\ \text{unit priors }\sigma_{\rm pre}=1\ }\
I^{\rm glob}_i = g_i = \mathrm{Cov}(\hat m,\theta_i).$$

This $I^{\rm glob}_i=g_i$ is exactly the code's marginal `dm_i` and the marginal forest bars.

**Additive property (the reason the method exists).** *When the auxiliary observables are
independent* ($C$, hence $\texttt{inv\_C}$, diagonal),

$$\big(\sigma_{\rm syst}^{\rm glob}\big)^2=\sum_i \big(I^{\rm glob}_i\big)^2=\sum_i g_i^2 .$$

**Stat/syst interpretation (paper slide 8).** The global split books the *in-situ
statistical constraint of the nuisances* as **statistical**, so
$\sigma_{\rm syst}^{\rm glob}=$ fluctuations of the external measurements only, and
$\sigma_{\rm stat}^{\rm glob}=\sqrt{\sigma_m^2-(\sigma_{\rm syst}^{\rm glob})^2}$ absorbs
that in-situ constraint. Same total, different stat/syst line than the freeze method —
neither is "more correct".

**Caveat under correlated priors.** With $\texttt{inv\_C}$ non-diagonal (our `corrV2`
case) the plain quadrature $\sum_i g_i^2$ is **not** $\sigma_{\rm syst}^2$ — see §C.

*Live numbers:* $\sqrt{\sum_i g_i^2}=954.7$ MeV, which even **exceeds** the experimental
total (852.1 MeV) — a clear symptom that quadrature-summing the marginals is wrong here.

---

## C. Global impact — exact additive variance share (correlated priors)

**The correct systematic variance.** Propagate fluctuations of the auxiliary data $y$
(covariance = prior covariance $C=\texttt{inv\_C}^{-1}$) through the profiled fit. Using
$\partial \hat x/\partial y=\Sigma_{x\theta}\,\texttt{inv\_C}$,

$$\frac{\partial\hat m}{\partial y}=w^\top\Sigma_{x\theta}\,\texttt{inv\_C}=g^\top\texttt{inv\_C},$$

$$\boxed{\ \sigma_{\rm syst}^2=\frac{\partial\hat m}{\partial y}\,C\,\frac{\partial\hat m}{\partial y}^{\!\top}
=\big(g^\top\texttt{inv\_C}\big)\,\texttt{inv\_C}^{-1}\,\big(\texttt{inv\_C}\,g\big)
= g^\top\,\texttt{inv\_C}\,g\ }$$

This is positive ($\texttt{inv\_C}\succ0$) and provably $\le\sigma_m^2$; it reduces to
$\sum_i g_i^2$ only when $\texttt{inv\_C}=\mathbb{1}$.

**Exact additive attribution.** Define $a\equiv\texttt{inv\_C}\,g$. Then

$$\sigma_{\rm syst}^2=g^\top a=\sum_i \underbrace{g_i\,a_i}_{\textstyle s_i},\qquad
s_g=\sum_{i\in g} g_i\,a_i,\qquad
\boxed{\ \sum_i s_i=\sum_g s_g=\sigma_{\rm syst}^2\ (\text{exactly})\ }$$

Plot the signed root $\operatorname{sign}(s)\sqrt{|s|}$ so the (signed) squares still add to
$\sigma_{\rm syst}^2$.

**Why this is the paper's method, generalised.** With independent priors
$\texttt{inv\_C}=\mathbb{1}\Rightarrow a_i=g_i\Rightarrow s_i=g_i^2=\big(I^{\rm glob}_i\big)^2$
— the two agree term by term. Restoring the correlation metric $\texttt{inv\_C}$ is the
*only* change; it turns the oblique (correlated) axes into the right inner product.

**Geometric picture — why the marginals (§B) don't sum but these do.**
$\sigma_{\rm syst}$ is the *length* of $g$ in the $\texttt{inv\_C}$ metric. Cholesky
$\texttt{inv\_C}=LL^\top$ and $h=L^\top g$ give an **orthogonal** (whitened) basis where
Pythagoras holds exactly, $\sum_k h_k^2=g^\top\texttt{inv\_C}\,g$ (verified to $10^{-13}$).
The named nuisances are the *oblique* axes; summing the squares of their projections
$g_i^2$ ignores the cross terms
$\sum_{i\ne j} g_i(\texttt{inv\_C})_{ij}g_j$. The variance share $s_i=g_i a_i$ folds each
row's full cross term back in, restoring closure.

**Why a share can be negative.** $s_i<0\iff\operatorname{sign}(g_i)\ne\operatorname{sign}(a_i)$:
the source's own mass pull $g_i$ opposes the correlation-weighted pull
$a_i=(\texttt{inv\_C}\,g)_i$, i.e. it is *anti-correlated through the prior* with the
dominant systematics, so including it slightly **reduces** the total variance — the exact
analogue of a negative covariance term in an ordinary correlated error budget.

*Live numbers:* $\sqrt{g^\top\texttt{inv\_C}\,g}=830.3$ MeV (linear) $\le 852.1$ total;
per-group shares close to $\sigma_{\rm syst}=825.5$ MeV (HESSE-anchored) with **1 of 49**
groups negative (`ATLAS_13TeV_modelling_singletop`, $-1.9$ MeV: $g=-1.7$, $a=+2.2$);
per-nuisance shares close to the same, 23 of 614 negative.

---

## D. Side-by-side

| | **A. Freeze / classic** | **B. Global — marginal** | **C. Global — additive share** |
|---|---|---|---|
| Per-source formula | $I^{\rm freeze}_i=g_i/\sigma_{{\rm post},i}$ | $I^{\rm glob}_i=g_i/\sigma_{{\rm pre},i}$ | $s_i=g_i\,(\texttt{inv\_C}\,g)_i$ |
| Reads as | "shift if $\theta_i$ fixed to $\pm\sigma_{\rm post}$" | "shift if aux obs moves $\pm\sigma_{\rm pre}$" | "contribution of $\theta_i$ to $\sigma_{\rm syst}^2$" |
| $\sum$ over sources $=\sigma_{\rm syst}^2$? | ✗ (never) | ✓ only if priors independent | ✓ **always (exact)** |
| Sign | $\ge 0$ | $\ge 0$ | signed (can be $<0$) |
| Books in-situ NP constraint as | systematic | statistical | (attribution of the global $\sigma_{\rm syst}$) |
| Total exp. syst | $733.1$ MeV | (marginals overshoot: $954.7$) | $825.5$ MeV |
| Code output | `doExpStatSystBreakdown` (Tier-1) | `global_impacts_{groups,nuisances}` | `global_impacts_{groups,nuisances}_additive` |

**Relations:** $I^{\rm glob}_i = I^{\rm freeze}_i\cdot\sigma_{{\rm post},i}/\sigma_{{\rm pre},i}$;
and $s_i \xrightarrow{\ \texttt{inv\_C}\to\mathbb{1}\ } (I^{\rm glob}_i)^2$. Both stat/syst
splits close to the same experimental total ($847.2$ MeV here); they only draw the
stat/syst line in different places.

## E. Practical guidance

- **Headline stat vs syst:** either split is defensible; quote which convention (freeze
  vs global) is used. Freeze $=$ traditional; global $=$ "external-only systematic".
- **"Which sources dominate?"** — read the **marginal** bars (§B): honest "size of each
  source on its own". Do **not** quadrature-sum them under correlated priors.
- **"Budget that adds up to $\sigma_{\rm syst}$"** — use the **additive shares** (§C).
  Expect a few small negative entries (anti-correlated sources); that is correct.
- Everything here is the *experimental* systematic (in the combined-cross-section
  covariance). Theory (PDF, scale, interpolation) is added separately:
  $\sigma_{\rm total}^2=\sigma_{\rm exp}^2+\sigma_{\rm PDF}^2+\sigma_{\rm interp}^2+\sigma_{\rm scale}^2$.

---

*Reference:* A. Freitas et al., *"Uncertainty breakdown via global impacts"*,
arXiv:2307.04007 (CMS-week slides by A. Gilbert, slides 4–9). §A ≡ conditional method
(slide 3); §B ≡ global impact / marginal additive form (slides 5–6, independent-aux case);
§C ≡ its exact quadratic-form generalisation for correlated auxiliary observables.
