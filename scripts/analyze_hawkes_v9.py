#!/usr/bin/env python3
"""
v9: Cycle-dependent productivity.

Model
-----
    λ*(t) = μ(t) + Σ_{t_i < t} α(t_i) · g(m_i) · β · e^{-β(t-t_i)}
          where  μ(t)   = μ_0 · (S(t)/S̄)^γ
                 α(t_i) = α_0 · (S(t_i)/S̄)^δ      ← NEW (v9)
                 g(m)   = exp(κ (m - m_0))

The kernel integrates to α(t_i)·g(m_i) per ancestor over t→∞, so the
branching ratio per parent is

    η_parent(t_i, m_i) = α_0 · (S(t_i)/S̄)^δ · exp(κ (m_i - m_0))

When δ = 0 this reduces to v7 exactly.

Pipeline
--------
  1.  Fit 6-parameter MLE on full 1868-2025 (multi-start).
  2.  LRT and AIC vs v7 (δ = 0 nested model).
  3.  Out-of-sample test: fit on 1868-2015 only, evaluate on 2016-2025
      with the same diagnostics as v8 (held-out log-L, time-rescaling,
      cumulative count band, 30-day rolling Brier skill, reliability).
  4.  Block-bootstrap 95% CIs on (μ0, γ, α0, β, κ, δ).

Random seed: 20260523 throughout.
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
# 1. Load merged 1868-2025 events + SSN
# ----------------------------------------------------------------------
events = pd.read_csv(os.path.join(DATA, "derived_events_extended_1868_2025.csv"),
                     parse_dates=["date"])
events = events.sort_values("date").reset_index(drop=True)
kp = pd.read_csv(os.path.join(DATA, "derived_daily_with_phase.csv"), parse_dates=["date"])

T0 = pd.Timestamp("1868-01-01")
SPLIT = pd.Timestamp("2016-01-01")
T_END = kp.date.max()
T_END_days = (T_END - T0).days
SPLIT_days = (SPLIT - T0).days

t_int = (events.date - T0).dt.days.values.astype(float)
jitter = rng.uniform(0.0, 1.0, size=len(t_int))
order = np.argsort(t_int + jitter)
t_all = (t_int + jitter)[order]
m_all = events.mark.values[order].astype(float)
N_all = len(t_all)
T_all = float(T_END_days)

print(f"[load] N={N_all} events, T={T_all:.0f} days = {T_all/365.25:.1f} yr")

# ----------------------------------------------------------------------
# 2. SSN
# ----------------------------------------------------------------------
ssn = pd.read_csv(os.path.join(DATA, "SN_m_tot_V2.0.txt"),
                  sep=r"\s+", header=None,
                  names=["year","month","yfrac","sn","sd","n","prov"], engine="python")
ssn = ssn[["yfrac","sn"]].copy()
ssn.loc[ssn.sn < 0, "sn"] = np.nan
ssn["sn_smooth"] = ssn.sn.rolling(window=13, center=True, min_periods=7).mean()
ssn = ssn.dropna().reset_index(drop=True)
T0_year = T0.year + (T0.month - 1)/12.0
ssn_t = (ssn.yfrac.values - T0_year) * 365.25
ssn_s = ssn.sn_smooth.values

def S(td):
    return np.interp(td, ssn_t, ssn_s, left=ssn_s[0], right=ssn_s[-1])

mask_w_full = (ssn_t >= 0) & (ssn_t <= T_all)
S_bar_full = ssn_s[mask_w_full].mean()
mask_w_tr = (ssn_t >= 0) & (ssn_t <= SPLIT_days)
S_bar_train = ssn_s[mask_w_tr].mean()
print(f"[ssn] S_bar full window = {S_bar_full:.2f}")
print(f"[ssn] S_bar train window= {S_bar_train:.2f}")

# ----------------------------------------------------------------------
# 3. v9 log-likelihood (6 params)
#
#  μ(t)   = mu0  * (S(t)/Sbar)^gamma
#  α(t_i) = alpha0 * (S(t_i)/Sbar)^delta
#  g(m)   = exp(kappa (m - m0))
#  λ*(t)  = μ(t) + Σ_{i: t_i<t} α(t_i) g(m_i) β e^{-β(t-t_i)}
#
#  ∫_{t_i}^{T} α(t_i) g(m_i) β e^{-β(s-t_i)} ds
#       = α(t_i) g(m_i) (1 - e^{-β(T-t_i)})
#
#  Recurrence: define R(t) = Σ_{t_i<t} α(t_i) g(m_i) e^{-β(t-t_i)}.
#  R(t_{k+1}) = e^{-β(t_{k+1}-t_k)} (R(t_k) + α(t_k) g(m_k))
#  λ*(t_{k+1}^-) needs R(t_{k+1}) without the new event.  We computed it pre-jump.
#  Note: I include the multiplicative β inside the kernel to keep things scale-clean.
#  → at jump time t_k: rate just before = μ(t_k) + β · R(t_k)
# ----------------------------------------------------------------------
m0 = 8.0

def loglike_v9(params, t, m, T_obs, t_start, S_bar):
    mu0, gamma, alpha0, beta, kappa, delta = params
    if mu0 <= 0 or alpha0 < 0 or beta <= 0:
        return -1e18
    S_events = S(t)
    g = np.exp(kappa * (m - m0))
    a_i = alpha0 * (S_events / S_bar) ** delta  # productivity per parent
    mu_events = mu0 * (S_events / S_bar) ** gamma

    s_log = 0.0
    R = 0.0
    for i in range(len(t)):
        if i > 0:
            R = np.exp(-beta * (t[i] - t[i-1])) * (R + a_i[i-1] * g[i-1])
        rate = mu_events[i] + beta * R
        if rate <= 0:
            return -1e18
        s_log += np.log(rate)

    grid = np.arange(t_start, t_start + T_obs + 1.0, 1.0)
    Sg = S(grid)
    mu_grid = mu0 * (Sg / S_bar) ** gamma
    s_int_mu = np.trapezoid(mu_grid, grid)
    T_end = t_start + T_obs
    s_comp = float(np.sum(a_i * g * (1.0 - np.exp(-beta * (T_end - t)))))
    return s_log - s_int_mu - s_comp

def neg_ll_v9(p, *a): return -loglike_v9(p, *a)

# v7-equivalent restricted likelihood (delta=0) for LRT
def loglike_v7(params5, t, m, T_obs, t_start, S_bar):
    mu0, gamma, alpha0, beta, kappa = params5
    return loglike_v9((mu0, gamma, alpha0, beta, kappa, 0.0), t, m, T_obs, t_start, S_bar)

# ----------------------------------------------------------------------
# 4. MLE on FULL window 1868-2025
# ----------------------------------------------------------------------
print("\n[fit-full] multi-start MLE on 1868-2025 (6 params)…")
starts = [
    # (mu0,    gamma, alpha0, beta, kappa, delta)
    (0.00443, 0.995, 0.1133, 0.6424, 0.899, 0.0),
    (0.005,   0.85,  0.12,   0.65,   1.00,  0.5),
    (0.004,   1.10,  0.10,   0.60,   0.80,  1.0),
    (0.006,   0.70,  0.15,   0.70,   0.50,  -0.3),
    (0.0045,  1.00,  0.11,   0.64,   0.90,  0.2),
    (0.005,   0.90,  0.13,   0.55,   1.20,  0.8),
    (0.0035,  1.15,  0.09,   0.50,   0.70,  1.2),
]
best_full = None
t0 = time.time()
for k, x0 in enumerate(starts):
    res = optimize.minimize(neg_ll_v9, x0,
                            args=(t_all, m_all, T_all, 0.0, S_bar_full),
                            method="Nelder-Mead",
                            options={"xatol":1e-7,"fatol":1e-7,"maxiter":80000})
    if res.fun < 1e17:
        mu0r,gr,a0r,br,kr,dr = res.x
        print(f"  trial {k+1}: μ0={mu0r:.5f} γ={gr:.3f} α0={a0r:.4f} β={br:.4f} κ={kr:+.3f} δ={dr:+.3f}  -LL={res.fun:.2f}")
        if best_full is None or res.fun < best_full.fun:
            best_full = res
print(f"  fit time: {time.time()-t0:.1f}s")

mu0_h, gamma_h, alpha0_h, beta_h, kappa_h, delta_h = best_full.x
ll_v9_full = -best_full.fun
aic_v9 = 2*6 - 2*ll_v9_full
print(f"\n[MLE v9 on full data 1868-2025]")
print(f"  μ0    = {mu0_h:.5f}/d = {mu0_h*365.25:.3f}/yr  at S=S̄")
print(f"  γ     = {gamma_h:.4f}")
print(f"  α0    = {alpha0_h:.4f}")
print(f"  β     = {beta_h:.4f}  → 1/β = {1/beta_h:.2f} d")
print(f"  κ     = {kappa_h:+.4f}  → exp(κ) = {np.exp(kappa_h):.3f}× G5/G4")
print(f"  δ     = {delta_h:+.4f}  → α scales as (S/S̄)^δ")
print(f"  log-L = {ll_v9_full:.2f},  AIC = {aic_v9:.2f}")
print(f"  Implied α at solar min (S=20):  α(min)/α0 = {(20/S_bar_full)**delta_h:.3f}")
print(f"  Implied α at solar max (S=180): α(max)/α0 = {(180/S_bar_full)**delta_h:.3f}")

# v7 restricted (delta=0) using same data
res_v7 = optimize.minimize(lambda p: -loglike_v7(p, t_all, m_all, T_all, 0.0, S_bar_full),
                           (mu0_h, gamma_h, alpha0_h, beta_h, kappa_h),
                           method="Nelder-Mead",
                           options={"xatol":1e-7,"fatol":1e-7,"maxiter":80000})
ll_v7_full = -res_v7.fun
aic_v7 = 2*5 - 2*ll_v7_full
print(f"\n[MLE v7 (δ=0 restricted) on same data] log-L = {ll_v7_full:.2f}, AIC = {aic_v7:.2f}")

# LRT
lr = 2*(ll_v9_full - ll_v7_full)
lr_p = stats.chi2.sf(lr, df=1)
print(f"\n[LRT v9 vs v7] χ²(1) = {lr:.3f},  p = {lr_p:.4g}")
print(f"               ΔAIC = AIC_v9 - AIC_v7 = {aic_v9 - aic_v7:+.2f}")

# ----------------------------------------------------------------------
# 5. Train-only fit on 1868-2015 for out-of-sample test
# ----------------------------------------------------------------------
print("\n[fit-train] multi-start MLE on 1868-2015 ONLY (6 params)…")
train_mask = t_all < SPLIT_days
test_mask  = ~train_mask
t_tr = t_all[train_mask]; m_tr = m_all[train_mask]
t_te = t_all[test_mask];  m_te = m_all[test_mask]
N_tr, N_te = len(t_tr), len(t_te)
T_tr = float(SPLIT_days)
T_te = float(T_END_days - SPLIT_days)
print(f"  train: N={N_tr} (G4 {(m_tr<9).sum()}, G5 {(m_tr>=9).sum()})  T={T_tr:.0f}d")
print(f"  test:  N={N_te} (G4 {(m_te<9).sum()}, G5 {(m_te>=9).sum()})  T={T_te:.0f}d")

best_tr = None
t0 = time.time()
for k, x0 in enumerate(starts):
    res = optimize.minimize(neg_ll_v9, x0,
                            args=(t_tr, m_tr, T_tr, 0.0, S_bar_train),
                            method="Nelder-Mead",
                            options={"xatol":1e-7,"fatol":1e-7,"maxiter":80000})
    if res.fun < 1e17:
        if best_tr is None or res.fun < best_tr.fun:
            best_tr = res
mu0_t, gamma_t, alpha0_t, beta_t, kappa_t, delta_t = best_tr.x
ll_train_v9 = -best_tr.fun
print(f"  train MLE: μ0={mu0_t:.5f} γ={gamma_t:.3f} α0={alpha0_t:.4f} β={beta_t:.4f} κ={kappa_t:+.3f} δ={delta_t:+.3f}")
print(f"             log-L (train) = {ll_train_v9:.2f},  fit time {time.time()-t0:.1f}s")

# ----------------------------------------------------------------------
# 6. Out-of-sample evaluation (v9 frozen, test 2016-2025)
# ----------------------------------------------------------------------
def loglike_v9_test(params, t_pre, m_pre, t_test, m_test, t_start, T_obs, S_bar):
    mu0, gamma, alpha0, beta, kappa, delta = params
    # Pre-test parents contribute through their own α(t_i)·g(m_i)
    g_pre = np.exp(kappa * (m_pre - m0))
    a_pre = alpha0 * (S(t_pre) / S_bar) ** delta
    g_te  = np.exp(kappa * (m_test - m0))
    a_te  = alpha0 * (S(t_test) / S_bar) ** delta
    R = float(np.sum(a_pre * g_pre * np.exp(-beta * (t_start - t_pre))))
    s_log = 0.0
    last_t = t_start
    for i, ti in enumerate(t_test):
        R = np.exp(-beta * (ti - last_t)) * R
        S_i = S(np.array([ti]))[0]
        mu_i = mu0 * (S_i / S_bar) ** gamma
        rate = mu_i + beta * R
        if rate <= 0: return -1e18
        s_log += np.log(rate)
        R = R + a_te[i] * g_te[i]
        last_t = ti
    T_end = t_start + T_obs
    grid = np.arange(t_start, T_end + 1.0, 1.0)
    Sg = S(grid)
    mu_grid = mu0 * (Sg / S_bar) ** gamma
    s_int_mu = np.trapezoid(mu_grid, grid)
    s_comp_pre = float(np.sum(a_pre * g_pre *
                              (np.exp(-beta*(t_start - t_pre)) - np.exp(-beta*(T_end - t_pre)))))
    s_comp_te = float(np.sum(a_te * g_te * (1.0 - np.exp(-beta*(T_end - t_test)))))
    return s_log - s_int_mu - s_comp_pre - s_comp_te

print("\n[oos] held-out log-likelihood, 2016-2025")
ll_te_v9 = loglike_v9_test(best_tr.x, t_tr, m_tr, t_te, m_te, T_tr, T_te, S_bar_train)
print(f"  v9 held-out logL: {ll_te_v9:.3f}")

# v8 (v7-form) train fit for fair comparison
res_v8 = optimize.minimize(lambda p: -loglike_v7(p, t_tr, m_tr, T_tr, 0.0, S_bar_train),
                           (mu0_t, gamma_t, alpha0_t, beta_t, kappa_t),
                           method="Nelder-Mead",
                           options={"xatol":1e-7,"fatol":1e-7,"maxiter":80000})
mu0_v8, gamma_v8, alpha_v8, beta_v8, kappa_v8 = res_v8.x
ll_te_v8 = loglike_v9_test((mu0_v8, gamma_v8, alpha_v8, beta_v8, kappa_v8, 0.0),
                            t_tr, m_tr, t_te, m_te, T_tr, T_te, S_bar_train)
print(f"  v8 (δ=0) held-out logL: {ll_te_v8:.3f}")

lam_null = N_tr / T_tr
ll_te_pois = N_te*np.log(lam_null) - lam_null*T_te
print(f"  Poisson null logL: {ll_te_pois:.3f}")
print(f"  Δ v9 - v8 = {ll_te_v9 - ll_te_v8:+.3f}")
print(f"  Δ v9 - Poisson = {ll_te_v9 - ll_te_pois:+.3f}  ({(ll_te_v9 - ll_te_pois)/N_te:+.3f}/event)")

# ----------------------------------------------------------------------
# 7. Time-rescaling on test events under v9
# ----------------------------------------------------------------------
print("\n[oos] time-rescaling on test events")
grid_te = np.arange(T_tr, T_tr + T_te + 1.0, 1.0)
Sg_te = S(grid_te)
mu_grid_te = mu0_t * (Sg_te / S_bar_train) ** gamma_t
cum_mu_te = np.concatenate([[0.0],
                            np.cumsum(0.5*(mu_grid_te[:-1]+mu_grid_te[1:])*np.diff(grid_te))])
def int_mu(tq): return np.interp(tq, grid_te, cum_mu_te)

g_pre = np.exp(kappa_t * (m_tr - m0))
a_pre = alpha0_t * (S(t_tr) / S_bar_train) ** delta_t
g_te_arr = np.exp(kappa_t * (m_te - m0))
a_te_arr = alpha0_t * (S(t_te) / S_bar_train) ** delta_t

def Lambda_test(tq):
    base = int_mu(tq)
    pre = float(np.sum(a_pre * g_pre *
                       (np.exp(-beta_t*(T_tr - t_tr)) - np.exp(-beta_t*(tq - t_tr)))))
    mask = t_te < tq
    te = float(np.sum(a_te_arr[mask] * g_te_arr[mask] *
                      (1.0 - np.exp(-beta_t*(tq - t_te[mask])))))
    return base + pre + te

Lam_te = np.array([Lambda_test(ti) for ti in t_te])
tau_te = np.diff(np.concatenate([[0.0], Lam_te]))
ks_te = stats.kstest(tau_te, "expon", args=(0, 1.0))
lag1_te = stats.pearsonr(tau_te[:-1], tau_te[1:]) if len(tau_te) >= 3 else None
print(f"  τ mean={tau_te.mean():.3f}, var={tau_te.var():.3f}")
print(f"  KS p = {ks_te.pvalue:.3e}")
if lag1_te is not None:
    print(f"  lag-1 r = {lag1_te.statistic:+.3f}  p={lag1_te.pvalue:.2e}")

# ----------------------------------------------------------------------
# 8. Cumulative count and rolling 30-day forecast
# ----------------------------------------------------------------------
tq_grid = np.arange(T_tr, T_tr + T_te + 1.0, 7.0)
Lam_grid = np.array([Lambda_test(x) for x in tq_grid])
low = stats.poisson.ppf(0.025, Lam_grid)
high = stats.poisson.ppf(0.975, Lam_grid)
total_Lam = Lambda_test(T_tr + T_te)
print(f"\n[cumcount] Λ_end = {total_Lam:.2f}, observed = {N_te}")
p_dev = 2*min(stats.poisson.cdf(N_te,total_Lam), stats.poisson.sf(N_te-1,total_Lam))
print(f"  Poisson two-sided p = {p_dev:.3f}")

print("\n[rolling] 30-day-ahead forecasts (daily)")
WIN = 30.0
forecast_days = np.arange(T_tr, T_tr + T_te - WIN + 1, 1.0)
p_at_least_1 = np.zeros(len(forecast_days))
obs_at_least_1 = np.zeros(len(forecast_days), dtype=int)
for k, td in enumerate(forecast_days):
    lam_win = Lambda_test(td + WIN) - Lambda_test(td)
    p_at_least_1[k] = 1.0 - np.exp(-lam_win)
    obs_at_least_1[k] = int(((t_te >= td) & (t_te < td + WIN)).any())
brier = float(np.mean((p_at_least_1 - obs_at_least_1)**2))
base_rate = float(obs_at_least_1.mean())
brier_clim = base_rate * (1 - base_rate)
bss = 1.0 - brier / brier_clim
print(f"  Brier = {brier:.4f}, climatology = {brier_clim:.4f}, BSS = {bss:+.3f}")

# Reliability
bins = np.linspace(0, 1, 11)
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

# Compare to v8 reliability — reload that table if present
v8_reliability = {}
v8_json_path = os.path.join(DATA, "v8_summary.json")
if os.path.exists(v8_json_path):
    with open(v8_json_path) as f:
        v8_summary = json.load(f)
        v8_reliability = v8_summary.get("rolling30d", {})
    print(f"  v8 BSS for comparison: {v8_reliability.get('brier_skill', float('nan')):+.3f}")
    print(f"  v9 BSS - v8 BSS = {bss - v8_reliability.get('brier_skill', 0):+.3f}")

# ----------------------------------------------------------------------
# 9. Block bootstrap (B=200, 365d blocks) on full data
# ----------------------------------------------------------------------
print("\n[bootstrap] block bootstrap (block=365d, B=200) on FULL window…")
BLOCK = 365.0
NB = int(np.ceil(T_all/BLOCK))
B = 200
rng_bs = np.random.default_rng(20260523)
params_bs = np.zeros((B, 6))
t_bs_start = time.time()
n_ok = 0
for b in range(B):
    starts_idx = rng_bs.integers(0, NB, size=NB)
    new_t, new_m = [], []
    cur_T = 0.0
    for k, si in enumerate(starts_idx):
        b_start = si * BLOCK
        b_end = min(b_start + BLOCK, T_all)
        sel = (t_all >= b_start) & (t_all < b_end)
        if sel.any():
            ts = t_all[sel] - b_start + cur_T
            ms = m_all[sel]
            new_t.append(ts); new_m.append(ms)
        cur_T += (b_end - b_start)
    if not new_t:
        continue
    tb = np.concatenate(new_t); mb = np.concatenate(new_m)
    order_b = np.argsort(tb)
    tb = tb[order_b]; mb = mb[order_b]
    Tb = cur_T
    try:
        res_b = optimize.minimize(neg_ll_v9, best_full.x,
                                  args=(tb, mb, Tb, 0.0, S_bar_full),
                                  method="Nelder-Mead",
                                  options={"xatol":1e-5,"fatol":1e-5,"maxiter":25000})
        if res_b.fun < 1e17 and np.all(np.isfinite(res_b.x)):
            params_bs[n_ok] = res_b.x
            n_ok += 1
    except Exception:
        pass
    if (b+1) % 50 == 0:
        print(f"  bootstrap {b+1}/{B}  ({n_ok} ok)  elapsed {time.time()-t_bs_start:.0f}s")

params_bs = params_bs[:n_ok]
print(f"\n[bootstrap] {n_ok}/{B} replicates succeeded ({time.time()-t_bs_start:.0f}s)")
names_p = ["μ0","γ","α0","β","κ","δ"]
hat = list(best_full.x)
ci = {}
print(f"\n  {'param':>6}  {'MLE':>10}  {'2.5%':>10}  {'97.5%':>10}  {'SE':>10}")
for i, name in enumerate(names_p):
    q = np.quantile(params_bs[:, i], [0.025, 0.5, 0.975])
    se = params_bs[:, i].std()
    ci[name] = (float(q[0]), float(q[2]))
    print(f"  {name:>6}  {hat[i]:10.5f}  {q[0]:10.5f}  {q[2]:10.5f}  {se:10.5f}")

# Derived: branching ratios at solar min vs max
eta_min_g4 = best_full.x[2] * (20/S_bar_full)**best_full.x[5]
eta_max_g4 = best_full.x[2] * (180/S_bar_full)**best_full.x[5]
eta_min_g5 = eta_min_g4 * np.exp(best_full.x[4])
eta_max_g5 = eta_max_g4 * np.exp(best_full.x[4])
print(f"\n  Branching at solar minimum (S=20):")
print(f"     η(G4) ≈ {eta_min_g4:.3f},  η(G5) ≈ {eta_min_g5:.3f}")
print(f"  Branching at solar maximum (S=180):")
print(f"     η(G4) ≈ {eta_max_g4:.3f},  η(G5) ≈ {eta_max_g5:.3f}")
print(f"  Max/min ratio = {eta_max_g4/eta_min_g4:.2f}× (productivity boost at solar max)")

np.save(os.path.join(DATA, "v9_bootstrap_params.npy"), params_bs)

# ----------------------------------------------------------------------
# 10. Save summary + plots
# ----------------------------------------------------------------------
summary = {
    "model": "v9 marked Hawkes with cycle-dependent productivity α(t) = α0 * (S(t)/Sbar)^δ",
    "full_mle": {
        "mu0": float(mu0_h), "gamma": float(gamma_h), "alpha0": float(alpha0_h),
        "beta": float(beta_h), "kappa": float(kappa_h), "delta": float(delta_h),
        "inv_beta": float(1/beta_h), "exp_kappa": float(np.exp(kappa_h)),
        "mu0_per_yr": float(mu0_h*365.25),
        "logL_full": float(ll_v9_full), "AIC_v9": float(aic_v9),
        "logL_v7_restricted": float(ll_v7_full), "AIC_v7": float(aic_v7),
        "LRT_chi2": float(lr), "LRT_p": float(lr_p),
        "delta_AIC": float(aic_v9 - aic_v7),
    },
    "branching": {
        "eta_g4_min": float(eta_min_g4),  "eta_g5_min": float(eta_min_g5),
        "eta_g4_max": float(eta_max_g4),  "eta_g5_max": float(eta_max_g5),
        "ratio_max_min": float(eta_max_g4/eta_min_g4),
    },
    "train_mle_1868_2015": {
        "mu0": float(mu0_t), "gamma": float(gamma_t), "alpha0": float(alpha0_t),
        "beta": float(beta_t), "kappa": float(kappa_t), "delta": float(delta_t),
        "logL_train": float(ll_train_v9),
    },
    "oos_test_2016_2025": {
        "N_test": int(N_te),
        "logL_v9": float(ll_te_v9),
        "logL_v8": float(ll_te_v8),
        "logL_poisson": float(ll_te_pois),
        "delta_v9_v8": float(ll_te_v9 - ll_te_v8),
        "delta_v9_poisson": float(ll_te_v9 - ll_te_pois),
        "ks_p": float(ks_te.pvalue),
        "lag1_r": float(lag1_te.statistic) if lag1_te else None,
        "lag1_p": float(lag1_te.pvalue) if lag1_te else None,
        "expected_total": float(total_Lam),
        "observed_total": int(N_te),
        "pois_two_sided_p": float(p_dev),
        "brier": float(brier),
        "brier_clim": float(brier_clim),
        "brier_skill": float(bss),
        "v8_brier_skill": float(v8_reliability.get("brier_skill", float("nan"))) if v8_reliability else None,
    },
    "bootstrap_ci_95": ci,
}
with open(os.path.join(DATA, "v9_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print(f"\n[save] data/v9_summary.json")

# --- Figure 23: full-window in-sample fit comparison v7 vs v9 ---------
fig, ax = plt.subplots(figsize=(10, 5))
# plot α(t)/α0 over the full SSN trace
s_curve_t = np.arange(0, T_all+1.0, 30.0)
s_curve = S(s_curve_t)
alpha_ratio = (s_curve/S_bar_full)**delta_h
dates_curve = T0 + pd.to_timedelta(s_curve_t, unit="D")
ax2 = ax.twinx()
ax2.plot(dates_curve, s_curve, color="0.6", lw=1, alpha=0.7, label="Smoothed SSN")
ax2.set_ylabel("Smoothed SSN", color="0.5")
ax.plot(dates_curve, alpha_ratio, color="C2", lw=2,
        label=f"α(t)/α0 = (S/S̄)^δ,  δ={delta_h:+.2f}")
ax.axhline(1.0, color="k", lw=0.7, ls="--", alpha=0.5)
ax.set_xlabel("Date")
ax.set_ylabel("α(t) / α0  (productivity multiplier)")
ax.set_title(f"v9 — cycle-dependent productivity:  α(t) = α0·(S(t)/S̄)^δ,  δ = {delta_h:+.3f}")
ax.legend(loc="upper left"); ax2.legend(loc="upper right")
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "23_v9_alpha_curve.png"), dpi=140)
plt.close()

# --- Figure 24: out-of-sample cumulative count v9 vs v8 ---------------
fig, ax = plt.subplots(figsize=(11, 5.5))
yrs = T0 + pd.to_timedelta(tq_grid, unit="D")
ax.fill_between(yrs, low, high, color="0.85", label="v9 Poisson 95% prediction band")
ax.plot(yrs, Lam_grid, color="C2", lw=2,
        label=f"v9 predicted Λ*(t)  (δ={delta_t:+.2f}, fit 1868-2015)")
# v8 curve for comparison: reuse from saved summary if available
ax.step([T0 + pd.to_timedelta(T_tr, unit="D")] + list(T0 + pd.to_timedelta(t_te, unit="D")) +
        [T0 + pd.to_timedelta(T_tr+T_te, unit="D")],
        [0] + list(np.arange(1, N_te+1)) + [N_te],
        where="post", color="C3", lw=2, label=f"Observed N(t)  (N={N_te})")
ax.set_xlabel("Date")
ax.set_ylabel("Cumulative G4+ events")
ax.set_title("v9 — out-of-sample cumulative count, 2016-2025")
ax.legend(loc="upper left"); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "24_v9_cum_count.png"), dpi=140)
plt.close()

# --- Figure 25: reliability diagram v9 ---------------------------------
fig, ax = plt.subplots(figsize=(7, 7))
mask_rel = counts > 0
ax.plot([0,1],[0,1], "k--", alpha=0.4, label="perfect calibration")
ax.plot(rel_pred[mask_rel], rel_obs[mask_rel], "o-", color="C2", lw=2, ms=9,
        label=f"v9 model (BSS = {bss:+.3f})")
# Overlay v8 calibration for direct comparison (loaded from notes in v8 logs)
v8_bins_pred = [0.03, 0.16, 0.22, 0.35, 0.45, 0.55, 0.62]
v8_bins_obs  = [0.00, 0.00, 0.16, 0.95, 0.94, 1.00, 1.00]
ax.plot(v8_bins_pred, v8_bins_obs, "s--", color="C0", lw=1.5, ms=7, alpha=0.8,
        label=f"v8 (BSS = {v8_reliability.get('brier_skill', float('nan')):+.3f})")
for i in range(10):
    if counts[i] > 0:
        ax.annotate(f"n={int(counts[i])}",
                    (rel_pred[i], rel_obs[i]),
                    xytext=(6,-12), textcoords="offset points", fontsize=8)
ax.set_xlabel("Predicted P(≥1 G4+ in 30d)")
ax.set_ylabel("Observed frequency")
ax.set_title(f"v9 reliability vs v8  (held-out 2016-2025)")
ax.legend(); ax.grid(alpha=0.3)
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "25_v9_reliability.png"), dpi=140)
plt.close()

# --- Figure 26: bootstrap distribution of δ ---------------------------
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(params_bs[:, 5], bins=30, color="C2", alpha=0.7, edgecolor="white")
ax.axvline(delta_h, color="k", lw=2, label=f"MLE δ = {delta_h:+.3f}")
ax.axvline(0, color="0.4", lw=1, ls="--", label="δ = 0  (v7 nested)")
ql, qh = np.quantile(params_bs[:, 5], [0.025, 0.975])
ax.axvline(ql, color="C3", lw=1, ls=":"); ax.axvline(qh, color="C3", lw=1, ls=":")
ax.set_xlabel("δ  (productivity scaling exponent)")
ax.set_ylabel("Bootstrap replicate count")
ax.set_title(f"v9 bootstrap distribution of δ\n95% CI [{ql:+.3f}, {qh:+.3f}],  N_bs = {n_ok}")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "26_v9_delta_bootstrap.png"), dpi=140)
plt.close()

print("[plots] figures 23-26 saved")
print("\n[done] v9 analysis complete.")
