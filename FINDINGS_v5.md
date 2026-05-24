# FINDINGS v5 — Non-stationary Hawkes process: the SSN-driven background

**Defensive publication. Not an operational risk model. Not financial, engineering, or insurance advice.**
**Diatom Sky R&D — open methodology, open data, open prior art.**

> v4 showed that the G4+ event series is self-exciting (η = 0.284, exponential
> kernel, 1/β = 1.74 d). But v4's residuals retained two warning signs: KS p = 2.8e−3
> against Exp(1) and a lag-1 autocorrelation of r = +0.23 (p = 3.3e−4). Both
> indicated that *something else* on a long timescale was still left in the data.
> v5 adds the obvious candidate — the solar cycle itself — as a time-varying
> background, and the residuals collapse to white Exp(1) noise.

## Model

A univariate exponential-kernel Hawkes process with a **non-stationary background**
modulated by the 13-month-smoothed monthly mean sunspot number $S(t)$ from SILSO:

$$
\lambda(t) \;=\; \mu(t) \;+\; \sum_{t_i < t} \alpha\, e^{-\beta(t - t_i)},
\qquad
\mu(t) \;=\; \mu_0\,\bigl(\tfrac{S(t)}{\bar S}\bigr)^{\gamma}
$$

where $\bar S = 93.50$ is the mean smoothed SSN over the 1932–2025 window. Three
nested models are compared:

| Model | Parameters | Notes |
|---|---|---|
| v0 — Poisson | $\lambda$ | constant rate |
| v4 — stationary Hawkes | $(\mu,\alpha,\beta)$ | self-exciting, constant background |
| v5 — non-stationary Hawkes | $(\mu_0,\gamma,\alpha,\beta)$ | self-exciting + SSN-modulated background |

v5 nests v4 at $\gamma = 0$.

## MLE result

Ogata-recursive log-likelihood, daily-grid trapezoidal integration of $\mu(t)$,
multi-start Nelder–Mead from 8 distinct initializations. **All starts converged
to the same global optimum to 5 significant figures.**

| Parameter | Value | Interpretation |
|---|---|---|
| $\mu_0$ | 0.00548 /day | 2.00 events/yr at $S = \bar S$ (mean-cycle background) |
| $\gamma$ | **0.846** | sub-linear SSN response |
| $\alpha$ | 0.1707 | excitation amplitude |
| $\beta$ | 0.6413 /day | 1/β = **1.56 d** (excitation decay) |
| $\eta = \alpha/\beta$ | **0.266** | branching ratio (subcritical) |

| Diagnostic | v5 | v4 | Poisson |
|---|---|---|---|
| Log-likelihood | **−1295.11** | −1335.61 | −1460.88 |
| AIC | **2598.21** | 2677.23 | 2923.76 |
| ΔAIC vs v4 | — | — | — |
| ΔAIC (v5 − v4) | **−79.02** | — | — |
| ΔAIC (v5 − Poisson) | **−325.55** | — | — |
| LR (v5 vs v4) | **81.02**, χ²(1) p ≈ 0 | — | — |
| Time-rescaling KS p vs Exp(1) | **0.438** | 2.81e−3 | 4.5e−16 |
| Lag-1 autocorr of τ | **r = +0.019**, p = 0.77 | r = +0.228, p = 3.3e−4 | — |
| Background fraction (declust) | 73.4% | 71.6% | — |

**v5 wins decisively.** AIC drops by 79 units, the likelihood-ratio test
rejects v4 at any conventional level, and the rescaled inter-arrival residuals
become statistically indistinguishable from Exp(1) — both in marginal
distribution (KS p = 0.44, up from 2.8e−3) and in serial independence (lag-1
autocorrelation collapses from r = +0.23 to r = +0.02). The two warning signs
that v4 left on the table are both gone.

## What γ = 0.846 means

If the background were proportional to SSN, $\gamma$ would be 1. If it were
proportional to SSN² (a CME-flux-like scaling), $\gamma$ would be 2. The fitted
value 0.846 ± (informally O(0.1) by likelihood curvature) sits just *below* linear,
consistent with the picture that G4+ days are not caused by SSN itself but by
high-speed streams and CMEs whose *occurrence* tracks SSN strongly but whose
*geo-effectiveness* depends on additional factors (Russell–McPherron alignment,
IMF orientation, ICME shock structure) that are not in the smoothed SSN regressor.
A genuinely linear response would also be surprising because the smoothed SSN
saturates around solar max while G4+ counts continue to depend on whether the
specific CMEs that month happen to be Earth-directed.

The numerical reduction of η from 0.284 (v4) to 0.266 (v5) is the expected
consequence of moving long-timescale variance out of the kernel and into the
background: some of what v4 attributed to "delayed offspring" is now correctly
attributed to "the cycle was simply more active that year."

## Goodness-of-fit summary

The time-rescaling theorem says that if the model is correctly specified, then
$\tau_k = \Lambda(t_k) - \Lambda(t_{k-1})$ should be i.i.d. Exp(1) where
$\Lambda(t) = \int_0^t \lambda(s)\,ds$.

- **v0 Poisson**: KS p = 4.5e−16 — catastrophically wrong (this is the v3/v4 result).
- **v4 stationary Hawkes**: KS p = 2.8e−3 — better, but the tail of $\tau$
  drifts above y = x past quantile ~3, and the lag-1 autocorrelation of +0.23
  indicates a smooth periodic component (the solar cycle) still in the residuals.
- **v5 non-stationary Hawkes**: KS p = 0.44, lag-1 autocorrelation of +0.02 — **the
  residuals look like white Exp(1) noise**. See `figures/13_v5_residuals.png`.

Empirical $\tau$ variance is 1.211 (theoretical 1.000). This is the only
remaining mild departure: a ~20% over-dispersion in residuals, much smaller than
v4's residual variance of 2.1, and could be addressed in a future v6 by a power-law
kernel or marked Hawkes with G4 vs G5 magnitudes.

## Decadal hazard, re-computed

5,000-trial Ogata-thinning Monte Carlo with the fitted parameters. Two scenarios
because the answer now *depends on where in the solar cycle the decade lives*:

| Scenario | Mean G4+/decade | Var/Mean | 95% CI | P(≥4 G4+ in any 7-day window) |
|---|---|---|---|---|
| Poisson (v0) | 26.3 | 1.0 | — | — |
| v4 stationary Hawkes | 26.1 | 2.02 | — | — |
| **v5 random historical decade** | **26.6** | **3.17** | [10, 45] | **45.1%** |
| **v5 SC25-like decade (2016–2025)** | **17.6** | 1.90 | [8, 30] | **33.2%** |

Two things to notice:

1. **The mean count agrees** across all three Hawkes formulations (~26/decade),
   as it must — the integrated background plus offspring has to recover the
   observed long-run rate. The disagreement is entirely in the variance and tail
   structure.

2. **The Scenario A Var/Mean of 3.17 is dramatically larger than v4's 2.02.**
   This is because random-decade selection now adds **between-decade variance**
   from the SSN level itself, on top of the within-decade self-excitation. Put
   another way: v4 treated every decade as equally storm-prone; v5 says some
   decades are inherently 2–4× more storm-prone than others, and the difference
   is predictable in advance from the smoothed sunspot number.

3. **SC25 looks below average.** Cycle 25 has peaked near 2024 with a smoothed
   SSN of ~160 — close to but slightly under the long-run mean of 93.5 weighted
   by a full cycle's worth of SSN. Scenario B integrates a real (recent) SC25
   profile and returns mean = 17.6 G4+/decade, against the all-history scenario
   mean of 26.6. This is consistent with SC25 being a *moderate* cycle, not a
   record-breaker. **This is not a forecast for SC26**, which we cannot yet
   model — see caveats.

## Stochastic declustering

For each observed G4+ event, the posterior probability of being a background
event (immigrant) under v5 is saved in
`data/derived_hawkes_v5_declustering.csv`:

- 180.5 of 246 events (**73.4%**) attributed to background
- 65.5 of 246 events (**26.6%**) attributed to offspring of an earlier event

This is the closest the model comes to a falsifiable physical picture: ~3 in 4
G4+ days are "fresh" from solar activity that month, ~1 in 4 are aftershocks of
the preceding G4+ within the 1.56-day decay envelope.

## Caveats — important

This is an exploratory, defensive-publication analysis by an independent
researcher. None of it is operational risk advice.

1. **Univariate, magnitude-blind.** A G4 day and a G5 day are treated as identical.
   A marked Hawkes with G4/G5 differentiation is a natural v6.
2. **One kernel, one regressor.** The exponential kernel might be the wrong shape;
   Nurhan et al. (2021) suggested a power-law inter-arrival tail. SSN is the
   simplest regressor; F10.7, ap-monthly, or CME-rate proxies could perform
   better.
3. **Smoothed SSN is causally backward at the edges.** The 13-month centered smoother
   uses 6 future months; for forecasting one must use a one-sided smoother or
   a forecast SSN profile. Scenario B used the actual (now known) SC25 trace and is
   therefore a *hindcast under the realized cycle*, not a forecast for SC26.
4. **The framework does not address grid resilience.** A G4+ in 2026 with a
   modernized GIC-blocking transformer fleet is a very different event from
   a G4+ in 1989.
5. **No claim of novelty.** Non-stationary Hawkes / ETAS models are standard in
   seismology and have been applied to space weather in other forms; the specific
   combination of Ogata MLE with SSN-power-law background for *G4+ Kp-defined
   geomagnetic storms* is, to my knowledge, not in the published literature, but
   I am a hobbyist and may simply not have found it. Diatom Sky R&D publishes
   methodology in case it is useful as prior art.

## Outputs

| File | Description |
|---|---|
| `scripts/analyze_hawkes_nonstationary.py` | full v5 pipeline (MLE + GOF + sim) |
| `data/hawkes_v5_summary.txt` | text summary of all numbers above |
| `data/derived_hawkes_v5_declustering.csv` | per-event background probability under v5 |
| `data/run_v5_log.txt` | full console log of the v5 run |
| `figures/12_v5_mu_of_t.png` | μ(t) overlaid on G4+ event ticks and SSN |
| `figures/13_v5_residuals.png` | QQ-plot, v5 vs v4 residuals against Exp(1) |
| `figures/14_v5_decadal_counts.png` | decadal G4+ count distributions for all four models |

## Data and reproducibility

- Kp/ap: GFZ Potsdam, [`Kp_ap_since_1932.txt`](https://kp.gfz.de/app/files/Kp_ap_since_1932.txt)
- SSN: SILSO Royal Observatory of Belgium, [`SN_m_tot_V2.0.txt`](https://www.sidc.be/SILSO/DATA/SN_m_tot_V2.0.txt)
- Random seed: 20260523 (every stochastic step)
- License: MIT for code, CC0 for derived data and figures

## What's next (offered, not yet done)

- **v6 — power-law kernel** $\alpha (c + (t-t_i))^{-(p+1)}$ to test whether
  exponential is the right kernel shape
- **v6′ — marked Hawkes** with G4 vs G5 marks (19 G5 days available 1932–2025)
- **v7 — pre-1932 extension** via the aa-index (1868→) and the auroral record
  to triple the observation window and break ties between SSN regressors
- **A cross-check using Dst** instead of Kp to confirm the result is not
  specific to the Kp definition of "G4+"

— Diatom Sky R&D, 2026-05-23
