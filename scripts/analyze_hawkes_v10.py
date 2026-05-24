#!/usr/bin/env python3
"""
v10: Rolling-origin out-of-sample test.

v8 reported a Brier skill score of +0.426 for held-out 30-day forecasts on
2016-2025 after fitting v7 on 1868-2015. v10 asks: was that a single lucky
split (the SC25 maximum just happened to fall in the held-out window in a
way the model could anticipate), or is the forecast skill robust?

For each split-year s ∈ {1980, 1985, ..., 2015}:
  - Train v7 marked Hawkes on 1868-01-01 → (s)-12-31
  - Freeze parameters
  - Compute held-out logL, time-rescaling KS p, and rolling 30-day Brier
    skill score on (s+1)-01-01 → 2025-11-30
  - Repeat the v9 calendar-window post-mortem on the largest test windows

We then plot the BSS distribution across splits and compute its mean/IQR.

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
T_END = kp.date.max()
T_END_days = (T_END - T0).days

t_int = (events.date - T0).dt.days.values.astype(float)
jitter = rng.uniform(0.0, 1.0, size=len(t_int))
order = np.argsort(t_int + jitter)
t_all = (t_int + jitter)[order]
m_all = events.mark.values[order].astype(float)
N_all = len(t_all)

print(f"[load] N={N_all} events, T={T_END_days:.0f}d = {T_END_days/365.25:.1f}yr")

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

# ----------------------------------------------------------------------
# 2. v7 log-likelihood (5 params)
# ----------------------------------------------------------------------
m0 = 8.0

def loglike(params, t, m, T_obs, t_start, S_bar):
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
    grid = np.arange(t_start, t_start + T_obs + 1.0, 1.0)
    Sg = S(grid)
    mu_grid = mu0 * (Sg / S_bar) ** gamma
    s_int_mu = np.trapezoid(mu_grid, grid)
    T_end = t_start + T_obs
    s_comp = (alpha/beta) * np.sum(g * (1.0 - np.exp(-beta*(T_end - t))))
    return s_log - s_int_mu - s_comp

def loglike_test(params, t_pre, m_pre, t_test, m_test, t_start, T_obs, S_bar):
    mu0, gamma, alpha, beta, kappa = params
    g_pre = np.exp(kappa * (m_pre - m0))
    g_te  = np.exp(kappa * (m_test - m0))
    R = float(np.sum(g_pre * np.exp(-beta * (t_start - t_pre)))) if len(t_pre) else 0.0
    s_log = 0.0; last_t = t_start
    for i, ti in enumerate(t_test):
        R = np.exp(-beta * (ti - last_t)) * R
        S_i = S(np.array([ti]))[0]
        mu_i = mu0 * (S_i / S_bar) ** gamma
        rate = mu_i + alpha * R
        if rate <= 0: return -1e18
        s_log += np.log(rate)
        R = R + g_te[i]; last_t = ti
    T_end = t_start + T_obs
    grid = np.arange(t_start, T_end + 1.0, 1.0)
    Sg = S(grid)
    mu_grid = mu0 * (Sg / S_bar) ** gamma
    s_int_mu = np.trapezoid(mu_grid, grid)
    s_comp_pre = (alpha/beta) * np.sum(g_pre *
                                       (np.exp(-beta*(t_start - t_pre)) -
                                        np.exp(-beta*(T_end   - t_pre)))) if len(t_pre) else 0.0
    s_comp_te = (alpha/beta) * np.sum(g_te * (1.0 - np.exp(-beta*(T_end - t_test)))) if len(t_test) else 0.0
    return s_log - s_int_mu - s_comp_pre - s_comp_te

# ----------------------------------------------------------------------
# 3. Run rolling-origin OOS
# ----------------------------------------------------------------------
SPLIT_YEARS = [1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015]
results = []

start_params = [
    (0.005, 0.85, 0.12, 0.65, 1.00),
    (0.004, 1.10, 0.10, 0.60, 0.80),
    (0.006, 0.70, 0.15, 0.70, 0.50),
    (0.0045, 1.00, 0.11, 0.64, 0.90),
    (0.005, 0.90, 0.13, 0.55, 1.20),
]

print("\n[rolling] computing OOS for each split year…")
t_total = time.time()
for yr in SPLIT_YEARS:
    split_date = pd.Timestamp(f"{yr+1}-01-01")
    split_days = (split_date - T0).days
    train_mask = t_all < split_days
    test_mask  = ~train_mask
    t_tr = t_all[train_mask]; m_tr = m_all[train_mask]
    t_te = t_all[test_mask];  m_te = m_all[test_mask]
    N_tr, N_te = len(t_tr), len(t_te)
    T_tr = float(split_days)
    T_te = float(T_END_days - split_days)

    mask_w_tr = (ssn_t >= 0) & (ssn_t <= T_tr)
    S_bar_train = ssn_s[mask_w_tr].mean()

    # Fit on train (multi-start)
    t_fit = time.time()
    best = None
    for x0 in start_params:
        res = optimize.minimize(lambda p: -loglike(p, t_tr, m_tr, T_tr, 0.0, S_bar_train),
                                x0, method="Nelder-Mead",
                                options={"xatol":1e-7,"fatol":1e-7,"maxiter":60000})
        if res.fun < 1e17 and (best is None or res.fun < best.fun):
            best = res
    if best is None:
        print(f"  split {yr}: FIT FAILED")
        continue
    mu0_h, gamma_h, alpha_h, beta_h, kappa_h = best.x
    ll_train = -best.fun

    # Held-out logL
    ll_te = loglike_test(best.x, t_tr, m_tr, t_te, m_te, T_tr, T_te, S_bar_train)
    lam_null = N_tr / T_tr
    ll_te_pois = N_te*np.log(lam_null) - lam_null*T_te if N_te > 0 else 0.0

    # Time-rescaling on test
    grid_te = np.arange(T_tr, T_tr + T_te + 1.0, 1.0)
    Sg_te = S(grid_te)
    mu_grid = mu0_h * (Sg_te / S_bar_train) ** gamma_h
    cum_mu = np.concatenate([[0.0], np.cumsum(0.5*(mu_grid[:-1]+mu_grid[1:])*np.diff(grid_te))])
    def int_mu(tq): return np.interp(tq, grid_te, cum_mu)
    g_pre = np.exp(kappa_h * (m_tr - m0))
    g_te_arr = np.exp(kappa_h * (m_te - m0))
    def Lambda_test(tq):
        base = int_mu(tq)
        pre = (alpha_h/beta_h) * np.sum(g_pre *
                                        (np.exp(-beta_h*(T_tr - t_tr)) - np.exp(-beta_h*(tq - t_tr))))
        mask = t_te < tq
        te = (alpha_h/beta_h) * np.sum(g_te_arr[mask] *
                                       (1.0 - np.exp(-beta_h*(tq - t_te[mask]))))
        return base + pre + te
    Lam_te = np.array([Lambda_test(ti) for ti in t_te])
    tau_te = np.diff(np.concatenate([[0.0], Lam_te])) if N_te else np.array([])
    ks_p = stats.kstest(tau_te, "expon", args=(0,1.0)).pvalue if len(tau_te) >= 3 else np.nan
    lag1_r = stats.pearsonr(tau_te[:-1], tau_te[1:]).statistic if len(tau_te) >= 3 else np.nan

    # 30-day rolling Brier
    WIN = 30.0
    fc = np.arange(T_tr, T_tr + T_te - WIN + 1, 1.0)
    pred = np.zeros(len(fc)); obs = np.zeros(len(fc), dtype=int)
    for k, td in enumerate(fc):
        lam_win = Lambda_test(td + WIN) - Lambda_test(td)
        pred[k] = 1.0 - np.exp(-lam_win)
        obs[k]  = int(((t_te >= td) & (t_te < td + WIN)).any())
    brier = float(np.mean((pred - obs)**2)) if len(fc) else np.nan
    base_rate = float(obs.mean()) if len(fc) else np.nan
    brier_clim = base_rate*(1-base_rate) if base_rate not in (0,1) else np.nan
    bss = (1 - brier/brier_clim) if (brier_clim and brier_clim > 0) else np.nan

    # Expected vs observed count
    total_Lam = float(Lambda_test(T_tr + T_te)) if N_te else 0.0
    pois_p = 2*min(stats.poisson.cdf(N_te, total_Lam),
                   stats.poisson.sf(N_te-1, total_Lam)) if total_Lam > 0 else np.nan

    fit_time = time.time() - t_fit
    print(f"  split {yr}:  N_tr={N_tr:3d} N_te={N_te:3d} | "
          f"μ0={mu0_h*365.25:.2f}/yr γ={gamma_h:.2f} 1/β={1/beta_h:.2f}d "
          f"exp(κ)={np.exp(kappa_h):.2f}× | "
          f"Λ_pred={total_Lam:.1f} vs obs={N_te} | KS p={ks_p:.2f} | "
          f"BSS={bss:+.3f} | {fit_time:.0f}s")
    results.append({
        "split_year": yr, "split_date": str(split_date.date()),
        "N_train": int(N_tr), "N_test": int(N_te),
        "T_train_yr": T_tr/365.25, "T_test_yr": T_te/365.25,
        "mle": {"mu0": float(mu0_h), "gamma": float(gamma_h),
                "alpha": float(alpha_h), "beta": float(beta_h),
                "kappa": float(kappa_h),
                "inv_beta": float(1/beta_h),
                "exp_kappa": float(np.exp(kappa_h)),
                "mu0_per_yr": float(mu0_h*365.25)},
        "ll_train": float(ll_train),
        "ll_test_hawkes": float(ll_te),
        "ll_test_poisson": float(ll_te_pois),
        "delta_logl_per_event": float((ll_te - ll_te_pois)/N_te) if N_te else np.nan,
        "ks_p_test": float(ks_p),
        "lag1_r_test": float(lag1_r),
        "lambda_predicted": float(total_Lam),
        "pois_p_count": float(pois_p) if pois_p is not None else np.nan,
        "brier": float(brier), "brier_clim": float(brier_clim) if brier_clim else None,
        "bss": float(bss) if not np.isnan(bss) else None,
        "base_rate_test": float(base_rate),
        "mean_predicted": float(pred.mean()),
    })

print(f"\n[done] {len(results)}/{len(SPLIT_YEARS)} splits done in {time.time()-t_total:.0f}s")

# ----------------------------------------------------------------------
# 4. Aggregate
# ----------------------------------------------------------------------
df = pd.DataFrame(results)
bss_vals = df.bss.dropna().values
print("\n[summary]")
print(f"  n splits with valid BSS = {len(bss_vals)}")
if len(bss_vals):
    print(f"  BSS mean   = {bss_vals.mean():+.3f}")
    print(f"  BSS median = {np.median(bss_vals):+.3f}")
    print(f"  BSS IQR    = [{np.quantile(bss_vals,0.25):+.3f}, {np.quantile(bss_vals,0.75):+.3f}]")
    print(f"  BSS range  = [{bss_vals.min():+.3f}, {bss_vals.max():+.3f}]")
    print(f"  fraction BSS > 0 = {(bss_vals>0).mean():.2f}")
    print(f"  fraction BSS > 0.2 = {(bss_vals>0.2).mean():.2f}")

# Δ logL / event vs Poisson
dlpe = df.delta_logl_per_event.dropna().values
if len(dlpe):
    print(f"  ΔlogL/event vs Poisson, median = {np.median(dlpe):+.3f}, IQR "
          f"[{np.quantile(dlpe,0.25):+.3f}, {np.quantile(dlpe,0.75):+.3f}]")

# Parameter stability across splits
print("\n[parameter stability across splits]")
for k in ["mu0_per_yr", "gamma", "inv_beta", "exp_kappa"]:
    series = df.mle.apply(lambda d: d[k]).values
    print(f"  {k:>12}: median={np.median(series):.3f}  range [{series.min():.3f}, {series.max():.3f}]")

# Save
with open(os.path.join(DATA, "v10_rolling_summary.json"), "w") as f:
    json.dump({"splits": results,
               "aggregate": {
                   "n_splits": int(len(bss_vals)),
                   "bss_mean": float(np.nanmean(bss_vals)) if len(bss_vals) else None,
                   "bss_median": float(np.median(bss_vals)) if len(bss_vals) else None,
                   "bss_q25": float(np.quantile(bss_vals,0.25)) if len(bss_vals) else None,
                   "bss_q75": float(np.quantile(bss_vals,0.75)) if len(bss_vals) else None,
                   "bss_min": float(bss_vals.min()) if len(bss_vals) else None,
                   "bss_max": float(bss_vals.max()) if len(bss_vals) else None,
                   "frac_positive": float((bss_vals>0).mean()) if len(bss_vals) else None,
                   "frac_above_p2": float((bss_vals>0.2).mean()) if len(bss_vals) else None,
                   "dlpe_median": float(np.median(dlpe)) if len(dlpe) else None,
                }
              }, f, indent=2)
print("\n[save] data/v10_rolling_summary.json")

# ----------------------------------------------------------------------
# 5. Plots
# ----------------------------------------------------------------------
# Figure 28: BSS vs split year
fig, axes = plt.subplots(2, 2, figsize=(12, 9))
ax = axes[0,0]
ax.plot(df.split_year, df.bss, "o-", color="C0", lw=2, ms=10)
ax.axhline(0, color="0.5", ls="--", alpha=0.5)
ax.axhline(0.426, color="C3", ls=":", alpha=0.7, label="v8 reference (+0.426)")
ax.fill_between(df.split_year, 0, df.bss.values, alpha=0.15, color="C0")
ax.set_xlabel("Train/test split year"); ax.set_ylabel("Brier skill score (30-day)")
ax.set_title("v10 — BSS across rolling-origin splits")
ax.legend(); ax.grid(alpha=0.3)

ax = axes[0,1]
ax.plot(df.split_year, df.delta_logl_per_event, "s-", color="C2", lw=2, ms=10)
ax.axhline(0, color="0.5", ls="--", alpha=0.5)
ax.set_xlabel("Train/test split year"); ax.set_ylabel("ΔlogL / event vs Poisson")
ax.set_title("Held-out logL gain over constant-rate Poisson")
ax.grid(alpha=0.3)

ax = axes[1,0]
for k, color in [("mu0_per_yr", "C0"), ("gamma", "C1")]:
    series = df.mle.apply(lambda d: d[k]).values
    ax.plot(df.split_year, series, "o-", color=color, lw=2, ms=8, label=k)
ax.set_xlabel("Train end year"); ax.set_ylabel("parameter value")
ax.set_title("Parameter stability: μ₀ (events/yr) and γ")
ax.legend(); ax.grid(alpha=0.3)

ax = axes[1,1]
for k, color in [("inv_beta", "C3"), ("exp_kappa", "C2")]:
    series = df.mle.apply(lambda d: d[k]).values
    ax.plot(df.split_year, series, "o-", color=color, lw=2, ms=8, label=k)
ax.set_xlabel("Train end year"); ax.set_ylabel("parameter value")
ax.set_title("Parameter stability: 1/β (d) and exp(κ)")
ax.legend(); ax.grid(alpha=0.3)

plt.suptitle("v10 — rolling-origin out-of-sample diagnostics (v7 marked Hawkes)",
             fontsize=12, y=1.0)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "28_v10_rolling_oos.png"), dpi=140)
plt.close()

# Figure 29: predicted vs observed count for each split
fig, ax = plt.subplots(figsize=(9, 6))
ax.plot([0, max(df.lambda_predicted.max(), df.N_test.max())*1.1],
        [0, max(df.lambda_predicted.max(), df.N_test.max())*1.1],
        "k--", alpha=0.5, label="perfect")
for _, row in df.iterrows():
    ax.scatter(row.lambda_predicted, row.N_test, s=120, color="C0", zorder=3)
    ax.annotate(f"{row.split_year}", (row.lambda_predicted, row.N_test),
                xytext=(8, 8), textcoords="offset points", fontsize=10)
ax.set_xlabel("Predicted total Λ over test window")
ax.set_ylabel("Observed N over test window")
ax.set_title("v10 — calibration across rolling splits")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "29_v10_calibration.png"), dpi=140)
plt.close()

print("[plots] figures 28-29 saved")
print("\n[done] v10 complete")
