#!/usr/bin/env python3
"""
v8: Out-of-sample test of the v7 marked Hawkes model.

Strategy
--------
Fit the model on 1868-01-01 → 2015-12-31 only ("training set"), freeze all
parameters, then evaluate predictive performance on 2016-01-01 → present
("test set"). Solar Cycle 25 ramp + maximum (2020-2026) is entirely held out,
including the May 2024 Gannon storm. This is the gold-standard credibility
check for the framework.

Diagnostics
-----------
  1.  Held-out log-likelihood vs Poisson null (rate from training)
  2.  Time-rescaling theorem applied to TEST τ's only; KS test, QQ plot
  3.  Cumulative event count: predicted E[N(t)] vs observed N(t) over test window
  4.  Conditional intensity λ*(t) trace through test period with realized
      events overlaid -- did the model "see" the May 2024 cluster?
  5.  Forecast a 30-day-ahead probability of ≥1 G4+ event at each test day,
      then compute Brier score and reliability diagram against realized 30-day
      outcomes.
  6.  Mark predictions: observed G5 fraction vs model-implied G5 fraction.

All comparisons are zero-tuning: the model parameters never see 2016-2025 data.
"""

import os, time, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize, stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FIG  = os.path.join(ROOT, "figures")

rng = np.random.default_rng(20260523)

# ----------------------------------------------------------------------
# 1. Load merged 1868-2025 events (from v7) + Kp daily + SSN
# ----------------------------------------------------------------------
events = pd.read_csv(os.path.join(DATA, "derived_events_extended_1868_2025.csv"),
                     parse_dates=["date"])
events = events.sort_values("date").reset_index(drop=True)
kp = pd.read_csv(os.path.join(DATA, "derived_daily_with_phase.csv"), parse_dates=["date"])

T0 = pd.Timestamp("1868-01-01")
SPLIT = pd.Timestamp("2016-01-01")
T_END = kp.date.max()  # last day with Kp data
T_END_days = (T_END - T0).days
SPLIT_days = (SPLIT - T0).days
print(f"[load] events: N={len(events)}, range {events.date.min().date()} → {events.date.max().date()}")
print(f"[load] full window: {T_END_days} days ({T_END_days/365.25:.1f} yr)")
print(f"[split] train: 1868-01-01 → 2015-12-31  ({SPLIT_days} days, {SPLIT_days/365.25:.1f} yr)")
print(f"        test : 2016-01-01 → {T_END.date()}  ({T_END_days-SPLIT_days} days, "
      f"{(T_END_days-SPLIT_days)/365.25:.2f} yr)")

# Jitter ties (same scheme as v7)
t_int = (events.date - T0).dt.days.values.astype(float)
jitter = rng.uniform(0.0, 1.0, size=len(t_int))
order = np.argsort(t_int + jitter)
t_all = (t_int + jitter)[order]
m_all = events.mark.values[order].astype(float)
src_all = events.source.values[order]

# Split
train_mask = t_all < SPLIT_days
test_mask  = ~train_mask
t_tr = t_all[train_mask]; m_tr = m_all[train_mask]
t_te = t_all[test_mask];  m_te = m_all[test_mask]
N_tr, N_te = len(t_tr), len(t_te)
T_tr = float(SPLIT_days)
T_te = float(T_END_days - SPLIT_days)
print(f"[split] train events: {N_tr}  (G4 {int((m_tr<9).sum())}, G5 {int((m_tr>=9).sum())})")
print(f"        test  events: {N_te}  (G4 {int((m_te<9).sum())}, G5 {int((m_te>=9).sum())})")

# ----------------------------------------------------------------------
# 2. SSN modulation (same as v7) -- needs full window for evaluation
# ----------------------------------------------------------------------
ssn = pd.read_csv(os.path.join(DATA, "SN_m_tot_V2.0.txt"),
                  sep=r"\s+", header=None,
                  names=["year","month","yfrac","sn","sd","n","prov"],
                  engine="python")
ssn = ssn[["yfrac","sn"]].copy()
ssn.loc[ssn.sn < 0, "sn"] = np.nan
ssn["sn_smooth"] = ssn.sn.rolling(window=13, center=True, min_periods=7).mean()
ssn = ssn.dropna().reset_index(drop=True)
T0_year = T0.year + (T0.month - 1)/12.0
ssn_t = (ssn.yfrac.values - T0_year) * 365.25
ssn_s = ssn.sn_smooth.values

def S(td):
    return np.interp(td, ssn_t, ssn_s, left=ssn_s[0], right=ssn_s[-1])

# CRITICAL: S_bar must be the TRAINING-period mean to avoid leakage
mask_tr_w = (ssn_t >= 0) & (ssn_t <= T_tr)
S_bar_train = ssn_s[mask_tr_w].mean()
print(f"[ssn] Training-period mean smoothed SSN (1868-2015) = {S_bar_train:.2f}")

# ----------------------------------------------------------------------
# 3. Marked Hawkes log-likelihood (same form as v7)
# ----------------------------------------------------------------------
m0 = 8.0
def loglike(params, t, m, T_obs, t_start, S_bar):
    """Log-likelihood on window [t_start, t_start+T_obs]. All events must lie in window."""
    mu0, gamma, alpha, beta, kappa = params
    if mu0 <= 0 or alpha < 0 or beta <= 0: return -1e18
    S_events = S(t)
    mu_events = mu0 * (S_events / S_bar) ** gamma
    g = np.exp(kappa * (m - m0))
    R = 0.0; s_log = 0.0
    for i in range(len(t)):
        if i > 0:
            R = np.exp(-beta * (t[i] - t[i-1])) * (R + g[i-1])
        rate = mu_events[i] + alpha * R
        if rate <= 0: return -1e18
        s_log += np.log(rate)
    # Integral of mu over window
    grid = np.arange(t_start, t_start + T_obs + 1.0, 1.0)
    Sg = S(grid)
    mu_grid = mu0 * (Sg / S_bar) ** gamma
    s_int_mu = np.trapezoid(mu_grid, grid)
    T_end = t_start + T_obs
    s_comp = (alpha/beta) * np.sum(g * (1.0 - np.exp(-beta*(T_end - t))))
    return s_log - s_int_mu - s_comp

def neg_ll(p, *a): return -loglike(p, *a)

# ----------------------------------------------------------------------
# 4. MLE on TRAINING set only
# ----------------------------------------------------------------------
print("\n[train] multi-start MLE on 1868-2015 ONLY…")
starts = [
    (0.00549, 0.845, 0.111, 0.652, 1.00),
    (0.005,   0.7,   0.10,  0.50,  0.5),
    (0.006,   1.0,   0.15,  0.70,  1.2),
    (0.004,   1.2,   0.08,  0.30,  0.8),
    (0.007,   0.5,   0.18,  0.80,  0.3),
    (0.005,   0.9,   0.12,  0.60,  1.5),
]
best = None
t0 = time.time()
for k, x0 in enumerate(starts):
    res = optimize.minimize(neg_ll, x0,
                            args=(t_tr, m_tr, T_tr, 0.0, S_bar_train),
                            method="Nelder-Mead",
                            options={"xatol":1e-7,"fatol":1e-7,"maxiter":60000})
    if res.fun < 1e17:
        mu0r,gr,ar,br,kr = res.x
        print(f"  trial {k+1}: μ0={mu0r:.5f} γ={gr:.3f} α={ar:.4f} β={br:.4f} κ={kr:+.3f}  -LL={res.fun:.2f}")
        if best is None or res.fun < best.fun:
            best = res
print(f"  fit time: {time.time()-t0:.1f}s")

mu0_h, gamma_h, alpha_h, beta_h, kappa_h = best.x
ll_train = -best.fun
print(f"\n[train MLE]")
print(f"  μ0 = {mu0_h:.5f}/d = {mu0_h*365.25:.3f}/yr at S̄_train")
print(f"  γ  = {gamma_h:.4f}")
print(f"  1/β= {1/beta_h:.2f} d")
print(f"  exp(κ) = {np.exp(kappa_h):.3f}×")
print(f"  η(G4) = {alpha_h/beta_h:.3f},  η(G5) = {alpha_h*np.exp(kappa_h)/beta_h:.3f}")
print(f"  log-L (train) = {ll_train:.2f}")

# Compare with v7 (full-data) parameters
v7_full = (0.00443, 0.995, 0.1133, 0.6424, 0.899)
print(f"\n[compare] v7 full-data fit (for reference -- NOT used below):")
print(f"  μ0={v7_full[0]:.5f} γ={v7_full[1]:.3f} α={v7_full[2]:.4f} β={v7_full[3]:.4f} κ={v7_full[4]:+.3f}")

# ----------------------------------------------------------------------
# 5. HELD-OUT LOG-LIKELIHOOD on the test period
#    -- key trick: history before SPLIT contributes to λ*(t) via the kernel
# ----------------------------------------------------------------------
print("\n[test] held-out log-likelihood on 2016-2025")

def loglike_test(params, t_pre, m_pre, t_test, m_test, t_start, T_obs, S_bar):
    """LL on test window, with t_pre/m_pre history before t_start carried into λ*."""
    mu0, gamma, alpha, beta, kappa = params
    g_pre = np.exp(kappa * (m_pre - m0))
    g_te  = np.exp(kappa * (m_test - m0))
    # Initialise R at t_start using full pre-test history
    if len(t_pre):
        R = float(np.sum(g_pre * np.exp(-beta * (t_start - t_pre))))
    else:
        R = 0.0
    s_log = 0.0
    last_t = t_start
    for i, ti in enumerate(t_test):
        R = np.exp(-beta * (ti - last_t)) * R  # decay to ti
        S_i = S(np.array([ti]))[0]
        mu_i = mu0 * (S_i / S_bar) ** gamma
        rate = mu_i + alpha * R
        if rate <= 0: return -1e18
        s_log += np.log(rate)
        R = R + g_te[i]   # accumulate after observation
        last_t = ti
    # Compensator integral over [t_start, t_start+T_obs]
    T_end = t_start + T_obs
    grid = np.arange(t_start, T_end + 1.0, 1.0)
    Sg = S(grid)
    mu_grid = mu0 * (Sg / S_bar) ** gamma
    s_int_mu = np.trapezoid(mu_grid, grid)
    # Excitation from pre-test history over [t_start, T_end]
    if len(t_pre):
        s_comp_pre = (alpha/beta) * np.sum(g_pre *
                                           (np.exp(-beta*(t_start - t_pre)) -
                                            np.exp(-beta*(T_end   - t_pre))))
    else:
        s_comp_pre = 0.0
    # Excitation from test events
    s_comp_te = (alpha/beta) * np.sum(g_te * (1.0 - np.exp(-beta*(T_end - t_test))))
    return s_log - s_int_mu - s_comp_pre - s_comp_te

ll_te_hawkes = loglike_test(best.x, t_tr, m_tr, t_te, m_te, T_tr, T_te, S_bar_train)
print(f"  held-out log-L (frozen v8 Hawkes): {ll_te_hawkes:.3f}")

# Null Poisson: constant rate = N_tr/T_tr from training set
lam_null = N_tr / T_tr
ll_te_null = N_te * np.log(lam_null) - lam_null * T_te
print(f"  null Poisson rate (train): λ = {lam_null*365.25:.3f}/yr")
print(f"  held-out log-L (Poisson):  {ll_te_null:.3f}")
print(f"  ΔlogL (Hawkes − Poisson): {ll_te_hawkes - ll_te_null:+.3f}  "
      f"({(ll_te_hawkes - ll_te_null)/N_te:+.3f} per event)")

# Null Poisson with SSN modulation only (no excitation)
def loglike_test_ssn(params, t_test, t_start, T_obs, S_bar):
    mu0, gamma = params
    if mu0 <= 0: return -1e18
    S_i = S(t_test)
    mu_i = mu0 * (S_i / S_bar) ** gamma
    s_log = float(np.sum(np.log(mu_i)))
    grid = np.arange(t_start, t_start + T_obs + 1.0, 1.0)
    Sg = S(grid)
    mu_grid = mu0 * (Sg / S_bar) ** gamma
    s_int_mu = np.trapezoid(mu_grid, grid)
    return s_log - s_int_mu
ll_te_ssn = loglike_test_ssn((mu0_h, gamma_h), t_te, T_tr, T_te, S_bar_train)
print(f"  SSN-Poisson (frozen μ0,γ): logL = {ll_te_ssn:.3f}  "
      f"ΔlogL Hawkes − SSN-Poisson: {ll_te_hawkes - ll_te_ssn:+.3f}")

# ----------------------------------------------------------------------
# 6. Time-rescaling on the TEST window
# ----------------------------------------------------------------------
print("\n[rescale] time-rescaling theorem on test period")
grid_te = np.arange(T_tr, T_tr + T_te + 1.0, 1.0)
Sg_te = S(grid_te)
mu_grid_te = mu0_h * (Sg_te / S_bar_train) ** gamma_h
cum_mu_te = np.concatenate([[0.0],
                            np.cumsum(0.5*(mu_grid_te[:-1]+mu_grid_te[1:])*np.diff(grid_te))])
def int_mu(tq): return np.interp(tq, grid_te, cum_mu_te)

# Pre-test history contribution to Λ(ti)
g_pre = np.exp(kappa_h * (m_tr - m0))
def excite_pre(tq):
    # ∫_{T_tr}^{tq} α Σ_pre g_i exp(-β(s - t_i)) ds
    A = (alpha_h/beta_h) * np.sum(g_pre *
                                  (np.exp(-beta_h*(T_tr - t_tr)) -
                                   np.exp(-beta_h*(tq   - t_tr))))
    return A

g_te_arr = np.exp(kappa_h * (m_te - m0))
Lam_te = np.zeros(len(t_te))
G = 0.0; A = 0.0; last_t = T_tr
# A here will track Σ_pre g_i exp(-β(t - t_i)) summed only over pre-test events at time `last_t`
A_pre = float(np.sum(g_pre * np.exp(-beta_h * (T_tr - t_tr))))
for i, ti in enumerate(t_te):
    # Compensator contribution from baseline
    base = int_mu(ti)
    # Pre-test excitation cumulative from T_tr to ti
    pre_exc = excite_pre(ti)
    # Test-event excitation cumulative
    te_exc = (alpha_h/beta_h) * np.sum(g_te_arr[:i] *
                                       (1.0 - np.exp(-beta_h*(ti - t_te[:i]))))
    Lam_te[i] = base + pre_exc + te_exc

# τ_i = Λ(t_i) - Λ(t_{i-1}); under H_0 these should be i.i.d. Exp(1)
# First τ uses Λ(t_1) - Λ(T_tr); rest are differences.
tau_te = np.diff(np.concatenate([[0.0], Lam_te]))
ks_te = stats.kstest(tau_te, "expon", args=(0, 1.0))
if len(tau_te) >= 3:
    lag1_te = stats.pearsonr(tau_te[:-1], tau_te[1:])
    lag1_te_stat = lag1_te.statistic; lag1_te_p = lag1_te.pvalue
else:
    lag1_te_stat = np.nan; lag1_te_p = np.nan
print(f"  N_test = {N_te}, τ mean={tau_te.mean():.3f}, var={tau_te.var():.3f}")
print(f"  KS p (Exp(1)): {ks_te.pvalue:.3e}   lag-1 r: {lag1_te_stat:+.3f}  p={lag1_te_p:.2e}")

# ----------------------------------------------------------------------
# 7. Cumulative count: predicted E[N(t)] vs observed N(t) on test window
# ----------------------------------------------------------------------
print("\n[cumcount] expected vs observed cumulative count on test window")
# Use the compensator Λ*(t) absorbing history -- it equals E[N((T_tr,t])] under H_0.
def Lambda_test_at(tq):
    base = int_mu(tq)
    pre_exc = excite_pre(tq)
    # contribution of test-events that already occurred
    mask = t_te < tq
    te_exc = (alpha_h/beta_h) * np.sum(g_te_arr[mask] *
                                       (1.0 - np.exp(-beta_h*(tq - t_te[mask]))))
    return base + pre_exc + te_exc

tq_grid = np.arange(T_tr, T_tr + T_te + 1.0, 7.0)  # weekly
Lam_grid = np.array([Lambda_test_at(x) for x in tq_grid])
# Observed N(t): step at each event
obs_n = np.searchsorted(t_te, tq_grid, side="right")

# 95% prediction band on N(t) under inhomogeneous Poisson with intensity λ*(t):
# N(t) ~ Poisson(Λ*(t)). Use Λ* (compensator)
low = stats.poisson.ppf(0.025, Lam_grid)
high = stats.poisson.ppf(0.975, Lam_grid)
print(f"  At end of test window: E[N]={Lam_grid[-1]:.1f}, observed N={obs_n[-1]}")
print(f"  Poisson 95% band: [{low[-1]:.0f}, {high[-1]:.0f}]")

# Total compensator over test window vs total events
total_Lam = Lambda_test_at(T_tr + T_te)
p_one_sided = stats.poisson.sf(N_te - 1, total_Lam) if N_te >= total_Lam else stats.poisson.cdf(N_te, total_Lam)
print(f"  Total Λ over test = {total_Lam:.2f}, N_te = {N_te}")
print(f"  Poisson test p (two-sided, deviation from mean): "
      f"{2*min(stats.poisson.cdf(N_te,total_Lam), stats.poisson.sf(N_te-1,total_Lam)):.3f}")

# ----------------------------------------------------------------------
# 8. 30-day-ahead rolling forecast: probability of ≥1 G4+ event
# ----------------------------------------------------------------------
print("\n[rolling] 30-day-ahead forecasts (daily)")
WIN = 30.0
forecast_days = np.arange(T_tr, T_tr + T_te - WIN + 1, 1.0)
def expected_count(t_start, t_end):
    # E[N((t_start, t_end])] using compensator
    return Lambda_test_at(t_end) - Lambda_test_at(t_start)
p_at_least_1 = np.zeros(len(forecast_days))
obs_at_least_1 = np.zeros(len(forecast_days), dtype=int)
for k, td in enumerate(forecast_days):
    lam_win = expected_count(td, td + WIN)
    p_at_least_1[k] = 1.0 - np.exp(-lam_win)
    obs_at_least_1[k] = int(((t_te >= td) & (t_te < td + WIN)).any())

# Brier score
brier = float(np.mean((p_at_least_1 - obs_at_least_1)**2))
# Climatological Brier (probability = fraction of windows with ≥1)
base_rate = float(obs_at_least_1.mean())
brier_clim = base_rate * (1 - base_rate)
brier_skill = 1.0 - brier / brier_clim
# Marginal calibration: predicted mean vs observed mean
print(f"  N forecast days: {len(forecast_days)}")
print(f"  mean predicted P(≥1 in 30d): {p_at_least_1.mean():.3f}")
print(f"  observed rate of ≥1 in 30d:  {obs_at_least_1.mean():.3f}")
print(f"  Brier score: {brier:.4f}")
print(f"  Brier (climatology): {brier_clim:.4f}")
print(f"  Brier skill score: {brier_skill:+.3f}")

# Reliability diagram bins
bins = np.linspace(0, 1, 11)
mids = 0.5*(bins[:-1]+bins[1:])
rel_obs = np.zeros(10); rel_pred = np.zeros(10); counts = np.zeros(10)
for i in range(10):
    sel = (p_at_least_1 >= bins[i]) & (p_at_least_1 < bins[i+1])
    if sel.any():
        rel_obs[i] = obs_at_least_1[sel].mean()
        rel_pred[i] = p_at_least_1[sel].mean()
        counts[i] = sel.sum()
print("  reliability bins (pred mid → obs freq, n):")
for i in range(10):
    if counts[i] > 0:
        print(f"    [{bins[i]:.1f},{bins[i+1]:.1f}]: pred={rel_pred[i]:.2f} obs={rel_obs[i]:.2f} n={int(counts[i])}")

# ----------------------------------------------------------------------
# 9. Mark prediction: expected G5 fraction
# ----------------------------------------------------------------------
print("\n[marks] observed G5 fraction in test vs model implication")
g5_obs = float((m_te >= 9).sum())
g4_obs = float((m_te <  9).sum())
print(f"  observed: {int(g5_obs)} G5, {int(g4_obs)} G4-class  →  G5 frac = {g5_obs/N_te:.3f}")
# Training G5 fraction (climatology)
g5_tr = float((m_tr >= 9).sum())
g4_tr = float((m_tr <  9).sum())
print(f"  train climatology: {int(g5_tr)}/{int(N_tr)} = {g5_tr/N_tr:.3f}")
# Test of difference (binomial)
p_clim = g5_tr/N_tr
pval_g5 = 2*min(stats.binom.cdf(int(g5_obs), N_te, p_clim),
                stats.binom.sf(int(g5_obs)-1, N_te, p_clim))
print(f"  binomial test (test G5 frac vs train clim): p = {pval_g5:.3f}")

# ----------------------------------------------------------------------
# 10. Save summary + plots
# ----------------------------------------------------------------------
summary = {
    "split_date": str(SPLIT.date()),
    "train": {"N": int(N_tr), "T_days": float(T_tr), "T_years": T_tr/365.25,
              "G4": int((m_tr<9).sum()), "G5": int((m_tr>=9).sum())},
    "test":  {"N": int(N_te), "T_days": float(T_te), "T_years": T_te/365.25,
              "G4": int((m_te<9).sum()), "G5": int((m_te>=9).sum())},
    "train_mle": {"mu0": float(mu0_h), "gamma": float(gamma_h),
                  "alpha": float(alpha_h), "beta": float(beta_h),
                  "kappa": float(kappa_h),
                  "eta_g4": float(alpha_h/beta_h),
                  "eta_g5": float(alpha_h*np.exp(kappa_h)/beta_h),
                  "exp_kappa": float(np.exp(kappa_h)),
                  "inv_beta": float(1/beta_h),
                  "mu0_per_yr": float(mu0_h*365.25)},
    "v7_full": {"mu0": v7_full[0], "gamma": v7_full[1], "alpha": v7_full[2],
                "beta": v7_full[3], "kappa": v7_full[4]},
    "holdout_loglik": {"hawkes": float(ll_te_hawkes),
                       "poisson_null": float(ll_te_null),
                       "ssn_poisson": float(ll_te_ssn),
                       "delta_vs_poisson": float(ll_te_hawkes - ll_te_null),
                       "delta_vs_ssn_poisson": float(ll_te_hawkes - ll_te_ssn),
                       "delta_per_event_vs_poisson": float((ll_te_hawkes - ll_te_null)/N_te)},
    "time_rescaling": {"N": int(N_te), "tau_mean": float(tau_te.mean()),
                       "tau_var": float(tau_te.var()),
                       "ks_p": float(ks_te.pvalue),
                       "lag1_r": float(lag1_te_stat),
                       "lag1_p": float(lag1_te_p)},
    "cumcount": {"expected_total": float(total_Lam),
                 "observed": int(N_te),
                 "pois_band_low": float(low[-1]),
                 "pois_band_high": float(high[-1])},
    "rolling30d": {"n_days": int(len(forecast_days)),
                   "mean_pred": float(p_at_least_1.mean()),
                   "obs_rate": float(obs_at_least_1.mean()),
                   "brier": float(brier),
                   "brier_clim": float(brier_clim),
                   "brier_skill": float(brier_skill)},
    "marks": {"g5_obs": int(g5_obs), "n_test": int(N_te),
              "g5_frac_obs": float(g5_obs/N_te),
              "g5_frac_train": float(p_clim),
              "binom_p": float(pval_g5)},
}
with open(os.path.join(DATA, "v8_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print(f"\n[save] data/v8_summary.json written")

# --- Figure 19: cumulative E[N] vs observed N -------------------------
fig, ax = plt.subplots(figsize=(11, 5.5))
yrs = T0 + pd.to_timedelta(tq_grid, unit="D")
ax.fill_between(yrs, low, high, color="0.85", label="Poisson 95% prediction band")
ax.plot(yrs, Lam_grid, color="C0", lw=2, label="Predicted $\\Lambda^*(t)$  (frozen, fit 1868-2015)")
# Step plot of observed cumulative count
obs_steps_x = list(T0 + pd.to_timedelta(t_te, unit="D")) + [T0 + pd.to_timedelta(T_tr+T_te, unit="D")]
obs_steps_y = list(np.arange(1, N_te+1)) + [N_te]
ax.step([T0 + pd.to_timedelta(T_tr, unit="D")] + obs_steps_x,
        [0] + obs_steps_y, where="post", color="C3", lw=2,
        label=f"Observed $N(t)$  (held out, N={N_te})")
for tt, mm in zip(t_te, m_te):
    color = "darkred" if mm >= 9 else "C3"
    ax.axvline(T0 + pd.to_timedelta(tt, unit="D"), color=color, alpha=0.18, lw=0.8)
ax.set_xlabel("Date")
ax.set_ylabel("Cumulative G4+ events")
ax.set_title("v8 — out-of-sample forecast: predicted vs observed cumulative count, 2016-2025")
ax.legend(loc="upper left")
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "19_v8_cumulative_count.png"), dpi=140)
plt.close()
print("[plot] figures/19_v8_cumulative_count.png")

# --- Figure 20: QQ plot of test τ's vs Exp(1) -------------------------
fig, ax = plt.subplots(figsize=(6.5, 6.5))
n = len(tau_te)
theo = -np.log(1 - (np.arange(1, n+1) - 0.5)/n)
ax.plot(theo, np.sort(tau_te), "o", color="C0", ms=7)
mx = max(theo.max(), np.sort(tau_te)[-1])
ax.plot([0, mx], [0, mx], "k--", alpha=0.5)
ax.set_xlabel("Exp(1) theoretical quantiles")
ax.set_ylabel("Held-out τ quantiles")
ax.set_title(f"v8 — time-rescaling QQ plot (test 2016-2025, N={n})\n"
             f"KS p = {ks_te.pvalue:.3f}, lag-1 r = {lag1_te_stat:+.3f}")
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "20_v8_qq_holdout.png"), dpi=140)
plt.close()
print("[plot] figures/20_v8_qq_holdout.png")

# --- Figure 21: 30-day rolling forecast + observed events -------------
fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                         gridspec_kw={"height_ratios":[3,1]})
ax = axes[0]
dates_fc = T0 + pd.to_timedelta(forecast_days, unit="D")
ax.plot(dates_fc, p_at_least_1, color="C0", lw=1.2,
        label="Model P(≥1 G4+ in next 30 d)")
ax.axhline(base_rate, color="0.5", ls="--",
           label=f"Climatology (test base rate = {base_rate:.2f})")
ax.set_ylabel("P(≥1 G4+ in 30 d)")
ax.set_title(f"v8 — rolling 30-day forecast vs realized events (Brier skill = {brier_skill:+.3f})")
ax.legend(loc="upper left"); ax.grid(alpha=0.3); ax.set_ylim(0, 1.05)
ax2 = axes[1]
for tt, mm in zip(t_te, m_te):
    d = T0 + pd.to_timedelta(tt, unit="D")
    c = "darkred" if mm >= 9 else "C3"
    h = 1.0 if mm >= 9 else 0.7
    ax2.vlines(d, 0, h, colors=c, lw=2)
ax2.set_ylim(0, 1.1); ax2.set_yticks([])
ax2.set_xlabel("Date")
ax2.set_ylabel("Observed events", rotation=0, ha="right", va="center")
plt.tight_layout()
plt.savefig(os.path.join(FIG, "21_v8_rolling30d.png"), dpi=140)
plt.close()
print("[plot] figures/21_v8_rolling30d.png")

# --- Figure 22: reliability diagram -----------------------------------
fig, ax = plt.subplots(figsize=(6.5, 6.5))
mask_rel = counts > 0
ax.plot([0,1],[0,1], "k--", alpha=0.4, label="perfect calibration")
ax.plot(rel_pred[mask_rel], rel_obs[mask_rel], "o-", color="C0", lw=2, ms=9,
        label="model")
for i in range(10):
    if counts[i] > 0:
        ax.annotate(f"n={int(counts[i])}",
                    (rel_pred[i], rel_obs[i]),
                    xytext=(6,-12), textcoords="offset points", fontsize=8)
ax.set_xlabel("Predicted P(≥1 in 30d)")
ax.set_ylabel("Observed frequency")
ax.set_title(f"v8 — reliability diagram (Brier {brier:.3f}, skill {brier_skill:+.3f})")
ax.legend(); ax.grid(alpha=0.3)
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "22_v8_reliability.png"), dpi=140)
plt.close()
print("[plot] figures/22_v8_reliability.png")

print("\n[done] v8 out-of-sample evaluation complete.")
print("       summary JSON: data/v8_summary.json")
print("       plots: figures/19-22")
