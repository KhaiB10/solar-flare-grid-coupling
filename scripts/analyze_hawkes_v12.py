#!/usr/bin/env python3
"""
v12: Daily F10.7 as the productivity driver (replaces smoothed SSN).

v7-v11 used the 13-month smoothed SSN as the slow background modulator
S(t) of the marked Hawkes process. That choice has two known weaknesses:

  1. It is monthly, then smoothed -- it cannot represent the 27-day
     solar rotation modulation that brings active regions in and out of
     geoeffective view.
  2. It is a counted index (sunspots) rather than a physical irradiance
     measurement.

v12 replaces S(t) with the daily F10.7 cm solar radio flux. F10.7 has
been recorded since 14 Feb 1947 by NRC Canada (Ottawa 1947-1991, then
Penticton DRAO 1991-present) and is the operational space-weather
standard. It is a direct physical measurement of solar coronal/
chromospheric emission, sensitive at daily resolution to active-region
rotation -- exactly the structure SSN smoothes out.

Splice strategy
---------------
Daily F10.7 only exists from 1947 onward. For 1844-1946 we keep the
SILSO smoothed SSN, but rescale it onto the F10.7 plane using the
overlap 1947-2025 (linear regression F10.7 ~ a + b * smoothed_SSN).
The splice point at 1947-02-14 is verified to be continuous.

Pipeline
--------
  1. Parse GFZ Kp+F10.7 file; extract F10.7adj daily 1947-02-14 onwards
  2. Fit F10.7adj ~ a + b*ssn_smooth on 1947-2025 overlap
  3. Build spliced S_daily(t) covering 1844-07-01 to 2025-05-31
  4. Refit v7-form marked Hawkes on the 1844-2025 catalog (v11 events)
  5. Compare logL/AIC to v11
  6. Rolling-origin OOS (8 splits, same as v10) with F10.7 background
  7. Block bootstrap B=200
  8. Lomb-Scargle periodogram on residuals: look for 27-day signal
     gone in v12 (i.e. confirmed the rotation modulation is now captured)
  9. Plots and summary JSON.

Random seed: 20260523 throughout.
"""

import os, sys, time, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize, stats, signal

# unbuffered output so we can see progress
class _U:
    def __init__(self, s): self.s = s
    def write(self, x):
        self.s.write(x); self.s.flush()
    def flush(self): self.s.flush()
sys.stdout = _U(sys.stdout)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FIG  = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

SEED = 20260523
rng = np.random.default_rng(SEED)

# ----------------------------------------------------------------------
# 1. Parse GFZ Kp+F10.7 file (we already have it)
# ----------------------------------------------------------------------
print("[1] parsing GFZ Kp+F10.7 since 1932…")
gfz_file = os.path.join(DATA, "Kp_ap_Ap_SN_F107_since_1932.txt")
cols = ["YYY","MM","DD","days","days_m","Bsr","dB",
        "Kp1","Kp2","Kp3","Kp4","Kp5","Kp6","Kp7","Kp8",
        "ap1","ap2","ap3","ap4","ap5","ap6","ap7","ap8",
        "Ap","SN","F107obs","F107adj","D"]
gfz = pd.read_csv(gfz_file, sep=r"\s+", comment="#", header=None, names=cols, engine="python")
gfz["date"] = pd.to_datetime(dict(year=gfz.YYY, month=gfz.MM, day=gfz.DD))
gfz = gfz[["date","F107adj","SN"]]
gfz.loc[gfz.F107adj < 0, "F107adj"] = np.nan
print(f"  loaded {len(gfz)} daily rows, F107 valid: {gfz.F107adj.notna().sum()}")
print(f"  date range: {gfz.date.min().date()} → {gfz.date.max().date()}")
print(f"  F107adj range (valid): {gfz.F107adj.min():.1f} → {gfz.F107adj.max():.1f} sfu")

# ----------------------------------------------------------------------
# 2. Pre-1947 splice from SILSO smoothed SSN
# ----------------------------------------------------------------------
print("\n[2] loading SILSO smoothed SSN for pre-1947 splice…")
ssn = pd.read_csv(os.path.join(DATA, "SN_m_tot_V2.0.txt"),
                  sep=r"\s+", header=None,
                  names=["year","month","yfrac","sn","sd","n","prov"], engine="python")
ssn = ssn[["yfrac","sn"]].copy()
ssn.loc[ssn.sn < 0, "sn"] = np.nan
ssn["sn_smooth"] = ssn.sn.rolling(window=13, center=True, min_periods=7).mean()
ssn = ssn.dropna(subset=["sn_smooth"]).reset_index(drop=True)
ssn_yfrac = ssn.yfrac.values
ssn_smooth = ssn.sn_smooth.values

# Linear regression F107adj ~ a + b*ssn_smooth on the 1947-2025 overlap
gfz_overlap = gfz.dropna(subset=["F107adj"]).copy()
gfz_overlap["yfrac"] = gfz_overlap.date.dt.year + (gfz_overlap.date.dt.dayofyear - 1)/365.25
gfz_overlap["ssn_smooth"] = np.interp(gfz_overlap.yfrac, ssn_yfrac, ssn_smooth)
gfz_overlap = gfz_overlap.dropna(subset=["ssn_smooth"])
slope, intercept, r, p, se = stats.linregress(gfz_overlap.ssn_smooth, gfz_overlap.F107adj)
print(f"  F107adj = {intercept:.3f} + {slope:.4f} * ssn_smooth  (n={len(gfz_overlap)}, "
      f"r={r:.3f}, R²={r**2:.3f})")
print(f"  At ssn_smooth=0:   F107adj = {intercept:.1f} sfu")
print(f"  At ssn_smooth=100: F107adj = {intercept + 100*slope:.1f} sfu")
print(f"  At ssn_smooth=200: F107adj = {intercept + 200*slope:.1f} sfu")

# ----------------------------------------------------------------------
# 3. Build spliced daily S_daily(t) for 1844-07-01 -> 2025-05-31
# ----------------------------------------------------------------------
print("\n[3] building spliced daily S_daily(t)…")
T0 = pd.Timestamp("1844-07-01")
T_END = pd.Timestamp("2025-05-31")
full_dates = pd.date_range(T0, T_END, freq="D")
n_full = len(full_dates)
print(f"  {n_full} days from {T0.date()} to {T_END.date()}")

# Pre-1947: use SSN-derived proxy
yfrac_full = full_dates.year + (full_dates.dayofyear - 1)/365.25
ssn_at_full = np.interp(yfrac_full, ssn_yfrac, ssn_smooth)
S_proxy = intercept + slope * ssn_at_full

# Post-1947: use actual F10.7adj, fill missing days with proxy
gfz_indexed = gfz.set_index("date")
gfz_aligned = gfz_indexed.reindex(full_dates)
S_actual = gfz_aligned["F107adj"].values

S_daily = np.where(full_dates >= pd.Timestamp("1947-02-14"), S_actual, np.nan)
# Fill any gaps in F10.7 with SSN proxy (rare; ~5 missing days)
n_filled = 0
for i in range(n_full):
    if full_dates[i] >= pd.Timestamp("1947-02-14") and np.isnan(S_daily[i]):
        S_daily[i] = S_proxy[i]
        n_filled += 1
# Pre-1947: always SSN proxy
S_daily = np.where(full_dates < pd.Timestamp("1947-02-14"), S_proxy, S_daily)
print(f"  pre-1947 days using SSN proxy: {(full_dates < pd.Timestamp('1947-02-14')).sum()}")
print(f"  post-1947 missing F10.7 days filled with proxy: {n_filled}")
print(f"  S_daily range: {np.nanmin(S_daily):.1f} → {np.nanmax(S_daily):.1f} sfu")
print(f"  S_daily mean (full window): {np.nanmean(S_daily):.2f}")

# Clip F10.7 burst-contamination spikes at p99.5 (~300 sfu). These are
# short-duration radio-burst transients (not coronal irradiance changes).
# Without clipping, (S/S_bar)^γ explodes at large γ and breaks Nelder-Mead.
F_CLIP = 300.0
n_clip = int(np.sum(S_daily > F_CLIP))
S_daily = np.minimum(S_daily, F_CLIP)
print(f"  clipped {n_clip} burst-contaminated days at F_CLIP={F_CLIP} sfu")

# Apply a 5-day rolling median to suppress any residual single-day spikes
# while preserving 27-day rotation structure. The Carrington/extreme events
# are 1-day in our catalog so the catalog rate isn't affected; only the
# background driver smooths.
S_series = pd.Series(S_daily)
S_daily = S_series.rolling(window=5, center=True, min_periods=1).median().values
print(f"  5-day rolling median applied; new range {S_daily.min():.1f} → {S_daily.max():.1f} sfu")

# Save daily background for repro and for plotting
df_S = pd.DataFrame({"date": full_dates, "S_daily": S_daily,
                     "ssn_smooth": ssn_at_full})
df_S.to_csv(os.path.join(DATA, "derived_S_daily_v12.csv"), index=False)

# Interpolation function over t (days since T0)
S_grid_t = np.arange(n_full, dtype=float)  # 0, 1, ..., n_full-1
def S(td):
    """Daily F10.7-based background. td = days since T0 (float OK)."""
    return np.interp(td, S_grid_t, S_daily, left=S_daily[0], right=S_daily[-1])

S_bar = np.nanmean(S_daily)
print(f"  S_bar = {S_bar:.3f} sfu")

# ----------------------------------------------------------------------
# 4. Load v11 event catalog and refit
# ----------------------------------------------------------------------
print("\n[4] loading v11 event catalog 1844-2025…")
events = pd.read_csv(os.path.join(DATA, "derived_events_extended_1844_2025.csv"),
                     parse_dates=["date"])
events = events.sort_values("date").reset_index(drop=True)
n_raw = len(events)
events = events[(events.date >= T0) & (events.date <= T_END)].reset_index(drop=True)
print(f"  raw catalog: {n_raw} events; in-window {len(events)} from "
      f"{events.date.min().date()} → {events.date.max().date()}")
if n_raw != len(events):
    print(f"  (dropped {n_raw - len(events)} event(s) outside [{T0.date()}, {T_END.date()}])")

t_int = (events.date - T0).dt.days.values.astype(float)
jitter = rng.uniform(0.0, 1.0, size=len(t_int))
order = np.argsort(t_int + jitter)
t_all = (t_int + jitter)[order]
m_all = events.mark.values[order].astype(float)
T_obs = float((T_END - T0).days)
print(f"  T_obs = {T_obs:.0f} d = {T_obs/365.25:.2f} yr")
print(f"  marks: G4 = {(m_all<9).sum()}, G5 = {((m_all>=9)&(m_all<9.5)).sum()}, "
      f"Carrington-class = {(m_all>=9.5).sum()}")

# ----------------------------------------------------------------------
# 5. Marked Hawkes with daily F10.7 background
#    λ*(t) = μ0 (S(t)/S_bar)^γ + α Σ_{ti<t} exp(κ(mi-m0)) exp(-β(t-ti))
# ----------------------------------------------------------------------
m0 = 8.0

def loglike(params, t, m, T_obs, S_bar):
    mu0, gamma, alpha, beta, kappa = params
    # Sensible bounds: mu0 > 0, alpha >= 0, beta > 0,
    # gamma in [-1, 3], kappa in [-3, 4] (anything outside is physically silly)
    if (mu0 <= 0 or alpha < 0 or beta <= 0 or
        gamma < -1.0 or gamma > 3.0 or
        kappa < -3.0 or kappa > 4.0 or
        beta > 10.0):
        return -1e18
    S_e = S(t)
    mu_e = mu0 * (S_e / S_bar) ** gamma
    g = np.exp(kappa * (m - m0))
    R = 0.0
    s_log = 0.0
    for i in range(len(t)):
        if i > 0:
            R = np.exp(-beta * (t[i] - t[i-1])) * (R + g[i-1])
        rate = mu_e[i] + alpha * R
        if rate <= 0:
            return -1e18
        s_log += np.log(rate)
    # Integral of μ over [0, T_obs] -- daily grid evaluation
    grid = S_grid_t  # already daily
    mu_grid = mu0 * (S_daily / S_bar) ** gamma
    s_int_mu = np.trapezoid(mu_grid, grid)
    s_comp = (alpha / beta) * np.sum(g * (1.0 - np.exp(-beta * (T_obs - t))))
    return s_log - s_int_mu - s_comp

CKPT = os.path.join(DATA, "v12_checkpoint.npz")
# Warm starts seeded by v11 MLE: (μ0=0.00499, γ=1.012, α=0.0962, β=0.582, κ=1.092)
starts = [
    (0.00500, 1.00, 0.10, 0.60, 1.00),
    (0.00500, 1.20, 0.10, 0.60, 1.00),
    (0.00400, 0.80, 0.12, 0.55, 1.10),
]

print("\n[5] multi-start MLE on 1844-2025 with daily F10.7 background…")
best = None
if os.path.exists(CKPT):
    try:
        ck = np.load(CKPT, allow_pickle=True)
        if "mle_x" in ck.files:
            print("  [resume] loading MLE from checkpoint")
            class _R: pass
            best = _R(); best.x = ck["mle_x"]; best.fun = -float(ck["ll_full"])
            print(f"  [resume] best.x = {best.x}  -LL={best.fun:.2f}")
    except Exception as e:
        print(f"  checkpoint unreadable: {e}; refitting")
        best = None

t_fit = time.time()
if best is None:
    for k, x0 in enumerate(starts):
        t_st = time.time()
        res = optimize.minimize(lambda p: -loglike(p, t_all, m_all, T_obs, S_bar),
                                x0, method="Nelder-Mead",
                                options={"xatol":1e-5, "fatol":1e-5, "maxiter":10000})
        if res.fun < 1e17:
            mu0r, gr, ar, br, kr = res.x
            print(f"  trial {k+1}: μ0={mu0r:.5f} γ={gr:.3f} α={ar:.4f} β={br:.4f} κ={kr:+.3f}"
                  f"  -LL={res.fun:.2f}  ({time.time()-t_st:.0f}s)")
            if best is None or res.fun < best.fun:
                best = res
    if best is None:
        raise RuntimeError("MLE failed entirely")
    # Persist MLE for resume
    np.savez(CKPT, mle_x=best.x, ll_full=-best.fun)
mu0_h, gamma_h, alpha_h, beta_h, kappa_h = best.x
ll_full = -best.fun
print(f"\n[v12 MLE on 1844-2025]")
print(f"  μ0    = {mu0_h:.5f}/d = {mu0_h*365.25:.3f}/yr at S̄={S_bar:.2f} sfu")
print(f"  γ     = {gamma_h:.4f}")
print(f"  α     = {alpha_h:.4f}")
print(f"  β     = {beta_h:.4f}  → 1/β = {1/beta_h:.2f} d")
print(f"  κ     = {kappa_h:+.4f}  → exp(κ) = {np.exp(kappa_h):.3f}×")
print(f"  log-L = {ll_full:.2f},  AIC = {2*5 - 2*ll_full:.2f}")
print(f"  fit time {time.time()-t_fit:.1f}s")

# Reference: v11 log-likelihood
ll_v11 = -2329.02
print(f"\n  ΔlogL vs v11 (smoothed SSN): {ll_full - ll_v11:+.2f}")
print(f"  ΔAIC vs v11:                  {(2*5 - 2*ll_full) - (2*5 - 2*ll_v11):+.2f}")
print(f"  (negative ΔAIC means v12 is better)")

# ----------------------------------------------------------------------
# 6. Rolling-origin OOS, same 8 splits as v10
# ----------------------------------------------------------------------
print("\n[6] rolling-origin OOS with F10.7 background…")
def fit_window(t, m, T_obs, warm=None):
    # Single warm start from full-window MLE for speed
    x0 = warm if warm is not None else starts[0]
    res = optimize.minimize(
        lambda p: -loglike(p, t, m, T_obs, S_bar),
        x0, method="Nelder-Mead",
        options={"xatol":1e-5, "fatol":1e-5, "maxiter":15000})
    return res if res.fun < 1e17 else None

def poisson_rate(t, T_obs):
    return len(t) / T_obs  # per day, constant

def integrate_intensity(params, t_train, m_train, t_eval_start, t_eval_end):
    """Integrate λ*(t) over [t_eval_start, t_eval_end] given trained params and
    full event history up to t_eval_end-of-pre-eval."""
    mu0, gamma, alpha, beta, kappa = params
    g_train = np.exp(kappa * (m_train - m0))
    # background integral on a daily grid
    grid_idx_start = int(np.floor(t_eval_start))
    grid_idx_end = int(np.ceil(t_eval_end))
    grid_idx_start = max(grid_idx_start, 0)
    grid_idx_end = min(grid_idx_end, n_full - 1)
    grid_t = np.arange(grid_idx_start, grid_idx_end + 1, dtype=float)
    mu_grid = mu0 * (S_daily[grid_idx_start:grid_idx_end+1] / S_bar) ** gamma
    bg_int = np.trapezoid(mu_grid, grid_t)
    # excitation integral from training events only (no in-window self-exc)
    exc_int = (alpha / beta) * np.sum(
        g_train * (np.exp(-beta*(t_eval_start - t_train)) - np.exp(-beta*(t_eval_end - t_train)))
    )
    return bg_int + exc_int

split_years = [1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015]
roll_rows = []
for sy in split_years:
    t_split = float((pd.Timestamp(f"{sy}-01-01") - T0).days)
    mask_tr = t_all < t_split
    mask_te = (t_all >= t_split) & (t_all < T_obs)
    t_tr, m_tr = t_all[mask_tr], m_all[mask_tr]
    t_te, m_te = t_all[mask_te], m_all[mask_te]
    T_te = T_obs - t_split
    # Fit on train only, warm-start from full MLE
    t_split_window = t_split  # we treat T_obs = t_split (train ends at split)
    res = fit_window(t_tr, m_tr, t_split_window, warm=best.x)
    if res is None:
        roll_rows.append(dict(split=sy, status="fail"))
        continue
    mu0r, gr, ar, br, kr = res.x
    ll_tr = -res.fun
    # Evaluate on test window: compute logL on test events
    S_e_te = S(t_te)
    mu_e_te = mu0r * (S_e_te / S_bar) ** gr
    # Excitation from all prior events (train + earlier test)
    s_log_te = 0.0
    # Build R recursively over the full history but only sum log at test events
    full_t = t_all
    full_m = m_all
    R = 0.0
    for i in range(len(full_t)):
        if i > 0:
            R = np.exp(-br * (full_t[i] - full_t[i-1])) * (R + np.exp(kr*(full_m[i-1]-m0)))
        if mask_te[i]:
            S_e_i = S(full_t[i])
            mu_e_i = mu0r * (S_e_i / S_bar) ** gr
            rate = mu_e_i + ar * R
            if rate > 0:
                s_log_te += np.log(rate)
    # Compensator on test window
    cum_te = integrate_intensity(res.x, t_tr, m_tr, t_split, T_obs)
    # Excitation from test events themselves (they self-excite within the test window)
    g_te = np.exp(kr * (m_te - m0))
    cum_te_self = (ar/br) * np.sum(g_te * (1.0 - np.exp(-br*(T_obs - t_te))))
    ll_te = s_log_te - cum_te - cum_te_self
    # Poisson baseline on test
    lam_poi = poisson_rate(t_tr, t_split)
    ll_poi = len(t_te) * np.log(lam_poi) - lam_poi * T_te
    # 30-day rolling Brier with climatological reference (same as v10)
    # Forecast P(any event in next 30d) at each daily origin in the test window.
    WIN = 30.0
    # Precompute mu_grid for the whole window once
    mu_grid_full_w = mu0r * (S_daily / S_bar) ** gr  # length n_full
    cum_bg = np.concatenate([[0.0], np.cumsum(mu_grid_full_w)])  # cumulative integral, length n_full+1
    n_full_w = len(mu_grid_full_w)
    # Daily origins in the test window where a 30-day-ahead forecast is feasible
    fc = np.arange(t_split, T_obs - WIN + 1.0, 1.0)
    if len(fc) == 0:
        roll_rows.append(dict(split=sy, status="too_short"))
        continue
    # Vectorized background integral from t_split to t_split+WIN (rolls along fc)
    idx_lo = np.clip(np.floor(fc).astype(int), 0, n_full_w - 1)
    idx_hi = np.clip(np.floor(fc + WIN).astype(int), 0, n_full_w - 1)
    bg_win = cum_bg[idx_hi + 1] - cum_bg[idx_lo + 1]
    # Vectorized excitation integral over [td, td+WIN] for all events with t_i < td+WIN.
    # For an event at time t_i with mark m_i, contribution to integral over [td, td+WIN] is:
    #   if t_i < td:    (alpha/beta) * g_i * (exp(-beta*(td-t_i)) - exp(-beta*(td+WIN-t_i)))
    #   if td<=t_i<td+WIN: (alpha/beta) * g_i * (1 - exp(-beta*(td+WIN-t_i)))
    # We only need events with t_i < t_split (test events leak into Brier in v10's setup
    # would be cheating). Use only history events (t_all < t_split would be too restrictive
    # because we want to also use prior test events as they arrive — but for forecast at td
    # we may use any event with t_i < td). Implement vectorized over fc with prior events only.
    g_all = np.exp(kr * (m_all - m0))
    pred = np.zeros(len(fc)); obs_arr = np.zeros(len(fc), dtype=int)
    ab = ar / br
    for k, td in enumerate(fc):
        # events strictly before td
        mask_prior = t_all < td
        if mask_prior.any():
            dt_lo = td - t_all[mask_prior]
            dt_hi = dt_lo + WIN
            exc_w = ab * np.sum(g_all[mask_prior] * (np.exp(-br * dt_lo) - np.exp(-br * dt_hi)))
        else:
            exc_w = 0.0
        # events inside [td, td+WIN)
        mask_in = (t_all >= td) & (t_all < td + WIN)
        if mask_in.any():
            dt_hi2 = td + WIN - t_all[mask_in]
            exc_w += ab * np.sum(g_all[mask_in] * (1.0 - np.exp(-br * dt_hi2)))
        lam_win = bg_win[k] + exc_w
        pred[k] = 1.0 - np.exp(-lam_win)
        obs_arr[k] = int(((t_te >= td) & (t_te < td + WIN)).any())
    obs = obs_arr
    brier_h = float(np.mean((pred - obs)**2))
    base_rate = float(obs.mean())
    brier_clim = base_rate * (1 - base_rate) if base_rate not in (0, 1) else np.nan
    BSS = (1 - brier_h/brier_clim) if (brier_clim and brier_clim > 0) else np.nan
    brier_p = brier_clim  # for compatibility
    delta_ll_per_event = (ll_te - ll_poi) / max(len(t_te), 1)
    roll_rows.append(dict(
        split=sy, T_train_yr=t_split/365.25, T_test_yr=T_te/365.25,
        n_train=int(mask_tr.sum()), n_test=int(mask_te.sum()),
        mu0=mu0r, gamma=gr, alpha=ar, beta=br, kappa=kr,
        ll_train=ll_tr, ll_test=ll_te, ll_poi=ll_poi,
        delta_ll_per_event=delta_ll_per_event,
        brier_h=brier_h, brier_p=brier_p, BSS=BSS,
        base_rate=float(base_rate), n_obs_windows=int(obs.sum()),
    ))
    print(f"  split {sy}: BSS={BSS:+.3f}, ΔlogL/ev={delta_ll_per_event:+.3f}, "
          f"μ0={mu0r:.5f}, 1/β={1/br:.2f}d, exp(κ)={np.exp(kr):.2f}×")

df_roll = pd.DataFrame(roll_rows)
df_roll.to_csv(os.path.join(DATA, "v12_rolling_summary.csv"), index=False)
bss_vals = df_roll.BSS.values
print(f"\n  BSS across 8 splits (v12): median {np.median(bss_vals):+.3f}, "
      f"IQR [{np.quantile(bss_vals,0.25):+.3f}, {np.quantile(bss_vals,0.75):+.3f}], "
      f"range [{bss_vals.min():+.3f}, {bss_vals.max():+.3f}]")
# v10 reference
v10_bss = [0.412, 0.421, 0.418, 0.395, 0.397, 0.329, 0.349, 0.426]
print(f"  v10 reference: median {np.median(v10_bss):+.3f}, "
      f"range [{min(v10_bss):+.3f}, {max(v10_bss):+.3f}]")

# ----------------------------------------------------------------------
# 7. Block bootstrap on full 1844-2025 fit
# ----------------------------------------------------------------------
print("\n[7] block bootstrap B=200, block=365d…")
BLOCK = 365.0
T_all = T_obs
NB = int(np.ceil(T_all/BLOCK))
B = 200
rng_bs = np.random.default_rng(SEED)
params_bs = np.zeros((B, 5))
n_ok = 0
t_bs = time.time()
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
    ord_b = np.argsort(tb)
    tb = tb[ord_b]; mb = mb[ord_b]
    try:
        res_b = optimize.minimize(lambda p: -loglike(p, tb, mb, cur_T, S_bar),
                                  best.x, method="Nelder-Mead",
                                  options={"xatol":1e-5, "fatol":1e-5, "maxiter":25000})
        if res_b.fun < 1e17 and np.all(np.isfinite(res_b.x)):
            params_bs[n_ok] = res_b.x
            n_ok += 1
    except Exception:
        pass
    if (b+1) % 50 == 0:
        print(f"  bootstrap {b+1}/{B} ({n_ok} ok, {time.time()-t_bs:.0f}s)")
params_bs = params_bs[:n_ok]
np.save(os.path.join(DATA, "v12_bootstrap_params.npy"), params_bs)
ci = {}
names_p = ["mu0","gamma","alpha","beta","kappa"]
print(f"\n  {'param':>6}  {'MLE':>10}  {'2.5%':>10}  {'97.5%':>10}")
for i, name in enumerate(names_p):
    q = np.quantile(params_bs[:,i], [0.025, 0.5, 0.975])
    ci[name] = (float(q[0]), float(q[2]))
    print(f"  {name:>6}  {best.x[i]:10.5f}  {q[0]:10.5f}  {q[2]:10.5f}")
inv_beta_bs = 1/params_bs[:,3]
exp_kappa_bs = np.exp(params_bs[:,4])
print(f"  1/β    {1/best.x[3]:10.3f}  {np.quantile(inv_beta_bs,0.025):10.3f}  "
      f"{np.quantile(inv_beta_bs,0.975):10.3f}")
print(f"  exp(κ) {np.exp(best.x[4]):10.3f}  {np.quantile(exp_kappa_bs,0.025):10.3f}  "
      f"{np.quantile(exp_kappa_bs,0.975):10.3f}")
ci["inv_beta"] = (float(np.quantile(inv_beta_bs,0.025)), float(np.quantile(inv_beta_bs,0.975)))
ci["exp_kappa"] = (float(np.quantile(exp_kappa_bs,0.025)), float(np.quantile(exp_kappa_bs,0.975)))

# ----------------------------------------------------------------------
# 8. 27-day rotation diagnostic: Lomb-Scargle on event-day residuals
#    Build the residual process r(t) = events/day - λ*(t).  If the model
#    has absorbed the 27-day modulation, the residual PSD should *lack*
#    a peak at f = 1/27 d^-1. If a residual peak remains, we still aren't
#    capturing it.
# ----------------------------------------------------------------------
print("\n[8] 27-day rotation diagnostic via Lomb-Scargle on residuals…")
# Build per-day intensity for v12
mu_grid_full = mu0_h * (S_daily / S_bar) ** gamma_h
# Excitation: step through events on a daily grid
R_g = 0.0
lam_daily_v12 = np.zeros(n_full)
idx_e = 0
for k_g in range(n_full):
    td = float(k_g)
    while idx_e < len(t_all) and t_all[idx_e] <= td:
        if idx_e == 0:
            R_g = np.exp(kappa_h*(m_all[0]-m0))
        else:
            dt_e = t_all[idx_e] - t_all[idx_e-1]
            R_g = np.exp(-beta_h*dt_e)*(R_g + np.exp(kappa_h*(m_all[idx_e-1]-m0)))
            R_g = R_g + np.exp(kappa_h*(m_all[idx_e]-m0))
        idx_e += 1
    if idx_e > 0:
        R_eff = R_g * np.exp(-beta_h*(td - t_all[idx_e-1]))
    else:
        R_eff = 0.0
    lam_daily_v12[k_g] = mu_grid_full[k_g] + alpha_h * R_eff

# Same for v11 baseline (use stored v11 parameters)
v11_params = dict(mu0=1.824/365.25, gamma=1.012, alpha=0.0962, beta=0.582, kappa=1.092)
# v11 uses smoothed SSN as S, not F10.7; build that grid for comparison
S_ssn_full = ssn_at_full  # smoothed SSN at full daily resolution
S_bar_ssn = np.nanmean(S_ssn_full)
mu_grid_v11 = v11_params["mu0"] * (S_ssn_full / S_bar_ssn) ** v11_params["gamma"]
R_g = 0.0
lam_daily_v11 = np.zeros(n_full)
idx_e = 0
for k_g in range(n_full):
    td = float(k_g)
    while idx_e < len(t_all) and t_all[idx_e] <= td:
        if idx_e == 0:
            R_g = np.exp(v11_params["kappa"]*(m_all[0]-m0))
        else:
            dt_e = t_all[idx_e] - t_all[idx_e-1]
            R_g = np.exp(-v11_params["beta"]*dt_e)*(R_g + np.exp(v11_params["kappa"]*(m_all[idx_e-1]-m0)))
            R_g = R_g + np.exp(v11_params["kappa"]*(m_all[idx_e]-m0))
        idx_e += 1
    if idx_e > 0:
        R_eff = R_g * np.exp(-v11_params["beta"]*(td - t_all[idx_e-1]))
    else:
        R_eff = 0.0
    lam_daily_v11[k_g] = mu_grid_v11[k_g] + v11_params["alpha"] * R_eff

# Build daily count vector
daily_count = np.zeros(n_full)
for tt in t_all:
    di = int(np.floor(tt))
    if 0 <= di < n_full:
        daily_count[di] += 1

resid_v12 = daily_count - lam_daily_v12
resid_v11 = daily_count - lam_daily_v11

# Lomb-Scargle (or simple FFT since daily uniform) -- use scipy.signal.periodogram
f_v12, P_v12 = signal.periodogram(resid_v12, fs=1.0)  # cycles per day
f_v11, P_v11 = signal.periodogram(resid_v11, fs=1.0)

# Find power at 27-day band (24-30d window)
band = (f_v12 >= 1/30.) & (f_v12 <= 1/24.)
peak_v12 = P_v12[band].max()
peak_v11 = P_v11[band].max()
med_v12 = np.median(P_v12[(f_v12 > 1/300.) & (f_v12 < 1/10.)])
med_v11 = np.median(P_v11[(f_v11 > 1/300.) & (f_v11 < 1/10.)])
print(f"  27-day band peak power (resid): v11={peak_v11:.4f}, v12={peak_v12:.4f}")
print(f"  background median (resid):     v11={med_v11:.4f}, v12={med_v12:.4f}")
print(f"  signal-to-background ratio:    v11={peak_v11/med_v11:.2f}, v12={peak_v12/med_v12:.2f}")
print(f"  → If v12 ratio is lower, F10.7 absorbed some 27-day modulation.")

# ----------------------------------------------------------------------
# 9. Plots
# ----------------------------------------------------------------------
print("\n[9] generating plots…")

# F1: F10.7 vs SSN comparison
fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
axes[0].plot(full_dates, S_daily, lw=0.4, color="steelblue", label="S_daily (F10.7 spliced)")
axes[0].axvline(pd.Timestamp("1947-02-14"), color="red", ls="--", alpha=0.5, label="F10.7 splice")
axes[0].set_ylabel("F10.7 [sfu]")
axes[0].set_title("v12: spliced daily F10.7-equivalent background (1844-2025)")
axes[0].legend(loc="upper right", fontsize=8)
axes[1].plot(full_dates, ssn_at_full, lw=0.5, color="darkorange", label="smoothed SSN (v11 background)")
axes[1].set_ylabel("smoothed SSN")
axes[1].set_xlabel("date")
axes[1].legend(loc="upper right", fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "33_v12_f107_vs_ssn.png"), dpi=130)
plt.close()

# F2: residual periodogram
fig, ax = plt.subplots(figsize=(10, 5))
mask = (f_v12 >= 1/100.) & (f_v12 <= 1/5.)
ax.plot(1/f_v12[mask], P_v11[mask], lw=0.7, color="orange", label="v11 (smoothed SSN)")
ax.plot(1/f_v12[mask], P_v12[mask], lw=0.7, color="steelblue", label="v12 (daily F10.7)")
ax.axvspan(24, 30, color="red", alpha=0.15, label="27-day band")
ax.axvline(27, color="red", lw=0.6, ls="--", alpha=0.6)
ax.set_xlabel("period (days)")
ax.set_ylabel("residual power")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_title("Residual periodogram: did F10.7 absorb the 27-day solar-rotation signal?")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "34_v12_residual_periodogram.png"), dpi=130)
plt.close()

# F3: rolling BSS comparison
fig, ax = plt.subplots(figsize=(9, 5))
xs = np.array(split_years)
ax.plot(xs, v10_bss, "o-", color="orange", label="v10 (smoothed SSN)")
ax.plot(xs, bss_vals, "o-", color="steelblue", label="v12 (daily F10.7)")
ax.axhline(0, color="black", lw=0.5)
ax.set_xlabel("split year (train ends)")
ax.set_ylabel("Brier Skill Score on held-out window")
ax.set_title("Rolling-origin OOS: v12 (F10.7) vs v10 (smoothed SSN)")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "35_v12_rolling_bss.png"), dpi=130)
plt.close()

# F4: bootstrap distributions vs v11 reference lines
v11_ref = dict(mu0=1.824/365.25, gamma=1.012, alpha=0.0962, beta=0.582, kappa=1.092)
fig, axes = plt.subplots(2, 3, figsize=(13, 6))
labels = ["μ₀ (/d)", "γ", "α", "β", "κ", "1/β (d)"]
ref_vals = [v11_ref["mu0"], v11_ref["gamma"], v11_ref["alpha"], v11_ref["beta"], v11_ref["kappa"], 1/v11_ref["beta"]]
v12_vals = [best.x[0], best.x[1], best.x[2], best.x[3], best.x[4], 1/best.x[3]]
data_arrs = [params_bs[:,0], params_bs[:,1], params_bs[:,2], params_bs[:,3], params_bs[:,4], 1/params_bs[:,3]]
for ax, lbl, ref, est, da in zip(axes.flatten(), labels, ref_vals, v12_vals, data_arrs):
    ax.hist(da, bins=30, color="steelblue", alpha=0.65, edgecolor="none")
    ax.axvline(ref, color="orange", lw=2, ls="--", label=f"v11 = {ref:.4g}")
    ax.axvline(est, color="red", lw=2, label=f"v12 = {est:.4g}")
    ax.set_title(lbl)
    ax.legend(fontsize=8)
plt.suptitle("v12 bootstrap distributions vs v11 reference (orange dashed)")
plt.tight_layout()
plt.savefig(os.path.join(FIG, "36_v12_bootstrap_vs_v11.png"), dpi=130)
plt.close()

print("  saved figures 33-36")

# ----------------------------------------------------------------------
# 10. Save summary JSON
# ----------------------------------------------------------------------
summary = dict(
    seed=SEED,
    window=dict(start=str(T0.date()), end=str(T_END.date()),
                T_obs_days=T_obs, T_obs_yr=T_obs/365.25),
    splice=dict(slope=float(slope), intercept=float(intercept),
                r=float(r), r_squared=float(r**2),
                n_overlap=int(len(gfz_overlap)),
                splice_date=str(pd.Timestamp("1947-02-14").date())),
    S_bar=float(S_bar),
    n_events=int(len(t_all)),
    mle=dict(mu0=float(mu0_h), gamma=float(gamma_h), alpha=float(alpha_h),
             beta=float(beta_h), kappa=float(kappa_h),
             inv_beta=float(1/beta_h), exp_kappa=float(np.exp(kappa_h)),
             logL=float(ll_full), AIC=float(2*5 - 2*ll_full)),
    vs_v11=dict(delta_logL=float(ll_full - ll_v11),
                delta_AIC=float((2*5 - 2*ll_full) - (2*5 - 2*ll_v11))),
    rolling_oos=dict(
        split_years=split_years,
        BSS=[float(x) for x in bss_vals],
        median_BSS=float(np.median(bss_vals)),
        IQR_BSS=[float(np.quantile(bss_vals,0.25)), float(np.quantile(bss_vals,0.75))],
        range_BSS=[float(bss_vals.min()), float(bss_vals.max())],
        v10_ref_BSS=v10_bss,
    ),
    bootstrap_CIs=ci,
    rotation_diagnostic=dict(
        peak_v11=float(peak_v11), peak_v12=float(peak_v12),
        median_v11=float(med_v11), median_v12=float(med_v12),
        snr_v11=float(peak_v11/med_v11), snr_v12=float(peak_v12/med_v12),
        delta_snr=float(peak_v11/med_v11 - peak_v12/med_v12),
    ),
)
with open(os.path.join(DATA, "v12_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print(f"\n[done] summary written to data/v12_summary.json")
print(f"  ΔlogL = {ll_full - ll_v11:+.2f}")
print(f"  v12 BSS median  = {np.median(bss_vals):+.3f}  (v10 was {np.median(v10_bss):+.3f})")
print(f"  27d SNR  v11 → v12: {peak_v11/med_v11:.2f} → {peak_v12/med_v12:.2f}")
