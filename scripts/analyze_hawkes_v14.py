#!/usr/bin/env python3
"""
v14: Omori-Utsu power-law excitation kernel (replaces v12's exponential).

Background (S(t) = daily F10.7) is kept identical to v12. The only structural
change is the excitation triggering function:

    v12 (exponential):   phi(tau) = alpha * exp(-beta * tau)
    v14 (Omori-Utsu):    phi(tau) = alpha * (tau + c) ** (-p)         for tau > 0

So v14 adds one parameter (the power-law exponent p; c replaces 1/beta).
Marked-Hawkes form is:

    lambda*(t) = mu0 * (S(t)/S_bar)**gamma
                + alpha * sum_{ti < t} exp(kappa*(mi - m0)) * (t - ti + c)**(-p)

Parameters: (mu0, gamma, alpha, c, p, kappa)   -> 6 params (vs v12's 5)

Motivation
----------
v12 found that the 27-day Carrington-rotation residual signal was NOT
absorbed by swapping smoothed-SSN for daily F10.7 (SNR went 10.34 -> 10.39,
unchanged). That points to the rotation modulation living in the EXCITATION
KERNEL DECAY, not the background. An exponential decay with 1/beta = 1.7 d
piles up to a Fourier line near 27 d when many CMEs are spaced ~27 d apart.

The Omori-Utsu law is the earthquake-seismology standard for aftershock
decay and has a fatter tail than exponential -- residual rate at tau = 27 d
will be much smaller relative to tau = 1 d than for exponential. This should
either:
  (a) absorb the 27-day residual signal -> SNR drops materially
  (b) leave it untouched -> evidence the 27-day signal is in the background
      after all, and we need finer EUV proxies (future v16)
Either outcome is scientifically useful.

Random seed: 20260523 throughout.
"""

import os, sys, time, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize, signal

# unbuffered output for background runs
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
# 1. Load shared v12 background S(t) and event catalog
# ----------------------------------------------------------------------
print("[1] loading v12 spliced daily background S(t)…")
S_df = pd.read_csv(os.path.join(DATA, "derived_S_daily_v12.csv"), parse_dates=["date"])
S_df = S_df.sort_values("date").reset_index(drop=True)
T0 = S_df.date.iloc[0]
T_END = S_df.date.iloc[-1]
S_daily = S_df["S_daily"].values.astype(float)
n_full = len(S_daily)
S_grid_t = np.arange(n_full, dtype=float)  # days since T0
S_bar = float(np.mean(S_daily))
print(f"  S_daily: {n_full} days, {T0.date()} → {T_END.date()}, S̄={S_bar:.2f} sfu")

def S(t):
    """Linear interpolation of the daily background; clipped at endpoints."""
    return np.interp(t, S_grid_t, S_daily)

print("\n[2] loading event catalog 1844-2025…")
ev = pd.read_csv(os.path.join(DATA, "derived_events_extended_1844_2025.csv"),
                 parse_dates=["date"])
ev = ev.sort_values("date").reset_index(drop=True)
# Restrict to [T0, T_END] same as v12
mask_w = (ev["date"] >= T0) & (ev["date"] <= T_END)
n_dropped = int((~mask_w).sum())
ev = ev[mask_w].reset_index(drop=True)
t_all = ((ev["date"] - T0).dt.days).values.astype(float)
m_all = ev["mark"].values.astype(float)
T_obs = float((T_END - T0).days)
N = len(t_all)
print(f"  {N} events in window {T0.date()} → {T_END.date()} (dropped {n_dropped})")
print(f"  T_obs = {T_obs:.0f} d = {T_obs/365.25:.2f} yr")

# ----------------------------------------------------------------------
# 3. Omori-Utsu marked Hawkes
#    lambda*(t) = mu0 (S(t)/Sbar)^gamma + alpha sum_{ti<t} g_i (t-ti+c)^(-p)
#    Integral over [0,T] of phi(t-ti):
#      if p != 1: ((T-ti+c)^(1-p) - c^(1-p)) / (1-p)
#      if p == 1: log((T-ti+c)/c)
# ----------------------------------------------------------------------
m0 = 8.0

def omori_int(tau_end, c, p):
    """Integral of (s+c)^(-p) from s=0 to s=tau_end. tau_end is array-like or scalar."""
    tau_end = np.asarray(tau_end)
    if abs(p - 1.0) < 1e-9:
        return np.log((tau_end + c) / c)
    return ((tau_end + c) ** (1.0 - p) - c ** (1.0 - p)) / (1.0 - p)

def loglike(params, t, m, T_obs_w, S_bar_w):
    mu0, gamma, alpha, c, p, kappa = params
    # Bounds: mu0>0, alpha>=0, c>0, p>1 (so tail integral is finite as T->inf;
    # but for finite T any p>0 works), gamma in [-1,3], kappa in [-3,4]
    if (mu0 <= 0 or alpha < 0 or c <= 1e-4 or c > 30.0 or
        p <= 0.05 or p > 4.0 or
        gamma < -1.0 or gamma > 4.0 or
        kappa < -3.0 or kappa > 4.0):
        return -1e18
    S_e = S(t)
    mu_e = mu0 * (S_e / S_bar_w) ** gamma
    g = np.exp(kappa * (m - m0))
    # Log-sum at events: for each event i, lambda(ti) = mu_e[i] + alpha * sum_{j<i} g_j (ti-tj+c)^(-p)
    s_log = 0.0
    for i in range(len(t)):
        if i == 0:
            R = 0.0
        else:
            tau = t[i] - t[:i]
            R = float(np.sum(g[:i] * (tau + c) ** (-p)))
        rate = mu_e[i] + alpha * R
        if rate <= 0:
            return -1e18
        s_log += np.log(rate)
    # Background integral over [0, T_obs] on daily grid
    mu_grid = mu0 * (S_daily / S_bar_w) ** gamma
    s_int_mu = float(np.trapezoid(mu_grid, S_grid_t))
    # Excitation compensator: alpha * sum_i g_i * Phi(T - t_i; c, p)
    tau_remain = T_obs_w - t
    s_comp = alpha * float(np.sum(g * omori_int(tau_remain, c, p)))
    return s_log - s_int_mu - s_comp

# Warm starts. We start near v12 in terms of background; alpha/c/p chosen so
# the kernel's "characteristic" tail (where 50% of triggering integral lives)
# is ~1-2 days, matching v12's 1/beta=1.72 d.
# For Omori with p~1.5, c~0.5 d, the half-integral time is ~1-2 d.
starts = [
    # mu0,   gamma,  alpha, c,    p,    kappa
    (0.00450, 1.80, 0.080, 0.50, 1.40, 1.00),
    (0.00450, 1.80, 0.100, 0.30, 1.30, 1.00),
    (0.00450, 2.00, 0.120, 0.80, 1.50, 1.00),
    (0.00450, 1.50, 0.060, 0.20, 1.20, 1.10),
]

CKPT = os.path.join(DATA, "v14_checkpoint.npz")
print("\n[3] multi-start MLE on 1844-2025 with Omori-Utsu kernel…")

best = None
if os.path.exists(CKPT):
    try:
        ck = np.load(CKPT, allow_pickle=True)
        if "mle_x" in ck.files:
            class _R: pass
            best = _R()
            best.x = ck["mle_x"]
            best.fun = -float(ck["ll_full"])
            print(f"  [resume] best.x = {best.x}  -LL={best.fun:.2f}")
    except Exception as e:
        print(f"  ckpt unreadable: {e}; refitting")
        best = None

t_fit = time.time()
if best is None:
    for k, x0 in enumerate(starts):
        t_st = time.time()
        res = optimize.minimize(lambda p: -loglike(p, t_all, m_all, T_obs, S_bar),
                                x0, method="Nelder-Mead",
                                options={"xatol":1e-5, "fatol":1e-5, "maxiter":20000})
        if res.fun < 1e17:
            mu0r, gr, ar, cr, pr, kr = res.x
            print(f"  trial {k+1}: μ0={mu0r:.5f} γ={gr:.3f} α={ar:.4f} c={cr:.3f} p={pr:.3f} κ={kr:+.3f}"
                  f"  -LL={res.fun:.2f}  ({time.time()-t_st:.0f}s)")
            if best is None or res.fun < best.fun:
                best = res
    if best is None:
        raise RuntimeError("MLE failed entirely")
    np.savez(CKPT, mle_x=best.x, ll_full=-best.fun)

mu0_h, gamma_h, alpha_h, c_h, p_h, kappa_h = best.x
ll_full = -best.fun
AIC = 2*6 - 2*ll_full
print(f"\n[v14 MLE on 1844-2025, Omori kernel]")
print(f"  μ0    = {mu0_h:.5f}/d = {mu0_h*365.25:.3f}/yr at S̄={S_bar:.2f} sfu")
print(f"  γ     = {gamma_h:.4f}")
print(f"  α     = {alpha_h:.4f}")
print(f"  c     = {c_h:.4f} d")
print(f"  p     = {p_h:.4f}")
print(f"  κ     = {kappa_h:+.4f}  → exp(κ) = {np.exp(kappa_h):.3f}×")
print(f"  log-L = {ll_full:.2f},  AIC = {AIC:.2f}")
print(f"  fit time {time.time()-t_fit:.1f}s")

# Reference: v12 log-likelihood and AIC
ll_v12 = -2318.61
AIC_v12 = 2*5 - 2*ll_v12
print(f"\n  ΔlogL vs v12 (exp kernel):    {ll_full - ll_v12:+.2f}")
print(f"  ΔAIC vs v12 (penalizes p):     {AIC - AIC_v12:+.2f}")
print(f"  (negative ΔAIC means v14 is better despite extra parameter)")

# Characteristic decay time: time s* where (s*+c)^(-p) = 0.5 * c^(-p)
# i.e., s* = c * (2^(1/p) - 1)
half_decay = c_h * (2.0 ** (1.0 / p_h) - 1.0)
# Median triggering time (where Phi(s; c, p) reaches half of Phi(inf; c, p))
# For p > 1, Phi(inf) = c^(1-p)/(p-1). Half = c^(1-p)/(2(p-1)) which is reached
# when (s+c)^(1-p) = c^(1-p)/2 -> s+c = c * 2^(1/(p-1))
if p_h > 1.0:
    s_med = c_h * (2.0 ** (1.0 / (p_h - 1.0)) - 1.0)
else:
    s_med = np.inf
print(f"  kernel half-amplitude time: {half_decay:.2f} d")
print(f"  kernel half-integral time:  {s_med:.2f} d  (v12 exp half-life: 1.72 d)")

# Branching ratio (for p>1): n = alpha * <g> * c^(1-p)/(p-1)
g_all = np.exp(kappa_h * (m_all - m0))
g_bar = float(np.mean(g_all))
if p_h > 1.0:
    branching_v14 = alpha_h * g_bar * c_h ** (1.0 - p_h) / (p_h - 1.0)
else:
    branching_v14 = np.nan
print(f"  branching ratio n = α<g> ∫₀^∞ (s+c)^(-p) ds = {branching_v14:.3f}")

# ----------------------------------------------------------------------
# 4. Rolling-origin OOS, same 8 splits as v10 and v12
# ----------------------------------------------------------------------
print("\n[4] rolling-origin OOS with Omori kernel…")

def fit_window(t_tr, m_tr, T_w, warm=None, multi=False):
    """Refit Omori-Hawkes on a train window. Returns best result or None."""
    if warm is None:
        starts_w = starts
    elif not multi:
        starts_w = [warm]
    else:
        starts_w = [warm] + starts[:2]
    bestr = None
    for x0 in starts_w:
        try:
            r = optimize.minimize(lambda p: -loglike(p, t_tr, m_tr, T_w, S_bar),
                                  x0, method="Nelder-Mead",
                                  options={"xatol":1e-5, "fatol":1e-5, "maxiter":20000})
            if r.fun < 1e17 and (bestr is None or r.fun < bestr.fun):
                bestr = r
        except Exception:
            continue
    return bestr

split_years = [1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015]
roll_rows = []
for sy in split_years:
    t_split = float((pd.Timestamp(f"{sy}-01-01") - T0).days)
    mask_tr = t_all < t_split
    mask_te = (t_all >= t_split) & (t_all < T_obs)
    t_tr, m_tr = t_all[mask_tr], m_all[mask_tr]
    t_te, m_te = t_all[mask_te], m_all[mask_te]
    T_te = T_obs - t_split

    res = fit_window(t_tr, m_tr, t_split, warm=best.x)
    if res is None:
        roll_rows.append(dict(split=sy, status="fail"))
        continue
    mu0r, gr, ar, cr, pr, kr = res.x
    ll_tr = -res.fun

    # Test-window log-likelihood: events in test, with kernel from all prior events
    s_log_te = 0.0
    g_full = np.exp(kr * (m_all - m0))
    for idx in np.where(mask_te)[0]:
        ti = t_all[idx]
        prior = t_all < ti
        if prior.any():
            tau = ti - t_all[prior]
            R = float(np.sum(g_full[prior] * (tau + cr) ** (-pr)))
        else:
            R = 0.0
        S_e_i = S(np.array([ti]))[0]
        mu_e_i = mu0r * (S_e_i / S_bar) ** gr
        rate = mu_e_i + ar * R
        if rate > 0:
            s_log_te += np.log(rate)
    # Compensator over test window:
    # background on test
    mu_grid = mu0r * (S_daily / S_bar) ** gr
    # indices for [t_split, T_obs] on daily grid
    i_lo = int(np.floor(t_split)); i_hi = n_full
    bg_te = float(np.trapezoid(mu_grid[i_lo:i_hi], S_grid_t[i_lo:i_hi]))
    # excitation compensator from ALL events in [0, T_obs], integrated over test only:
    # integral_{t_split}^{T} (s - ti + c)^(-p) ds  where s>ti
    # For ti < t_split: range is [t_split, T] -> Phi(T-ti) - Phi(t_split-ti)
    # For ti in [t_split, T]: range is [ti, T] -> Phi(T-ti) - Phi(0)=0
    Phi_T   = omori_int(T_obs - t_all, cr, pr)
    Phi_lo  = np.where(t_all < t_split,
                       omori_int(t_split - t_all, cr, pr),
                       0.0)
    # Mask: only events with ti < T_obs contribute (all of them do)
    exc_te = ar * float(np.sum(g_full * (Phi_T - Phi_lo)))
    ll_te = s_log_te - bg_te - exc_te

    # Poisson baseline on test
    lam_poi = len(t_tr) / t_split
    ll_poi = len(t_te) * np.log(lam_poi) - lam_poi * T_te

    # 30-day rolling Brier with climatological reference (same as v12/v10)
    WIN = 30.0
    cum_bg = np.concatenate([[0.0], np.cumsum(mu_grid)])
    n_full_w = len(mu_grid)
    fc = np.arange(t_split, T_obs - WIN + 1.0, 1.0)
    if len(fc) == 0:
        roll_rows.append(dict(split=sy, status="too_short"))
        continue
    idx_lo_a = np.clip(np.floor(fc).astype(int), 0, n_full_w - 1)
    idx_hi_a = np.clip(np.floor(fc + WIN).astype(int), 0, n_full_w - 1)
    bg_win = cum_bg[idx_hi_a + 1] - cum_bg[idx_lo_a + 1]
    pred = np.zeros(len(fc)); obs_arr = np.zeros(len(fc), dtype=int)
    for k_i, td in enumerate(fc):
        # Excitation integral over [td, td+WIN] for all events with ti < td+WIN
        mask_prior = t_all < td
        if mask_prior.any():
            tau_lo = td - t_all[mask_prior]
            tau_hi = tau_lo + WIN
            exc_w = ar * float(np.sum(g_full[mask_prior] *
                                       (omori_int(tau_hi, cr, pr) -
                                        omori_int(tau_lo, cr, pr))))
        else:
            exc_w = 0.0
        mask_in = (t_all >= td) & (t_all < td + WIN)
        if mask_in.any():
            tau_hi2 = td + WIN - t_all[mask_in]
            exc_w += ar * float(np.sum(g_full[mask_in] *
                                        omori_int(tau_hi2, cr, pr)))
        lam_win = bg_win[k_i] + exc_w
        pred[k_i] = 1.0 - np.exp(-max(lam_win, 0.0))
        obs_arr[k_i] = int(((t_te >= td) & (t_te < td + WIN)).any())
    obs = obs_arr
    brier_h = float(np.mean((pred - obs) ** 2))
    base_rate = float(obs.mean())
    brier_clim = base_rate * (1 - base_rate) if base_rate not in (0, 1) else np.nan
    BSS = (1 - brier_h / brier_clim) if (brier_clim and brier_clim > 0) else np.nan
    delta_ll_per_event = (ll_te - ll_poi) / max(len(t_te), 1)

    roll_rows.append(dict(
        split=sy, T_train_yr=t_split/365.25, T_test_yr=T_te/365.25,
        n_train=int(mask_tr.sum()), n_test=int(mask_te.sum()),
        mu0=mu0r, gamma=gr, alpha=ar, c=cr, p=pr, kappa=kr,
        ll_train=ll_tr, ll_test=ll_te, ll_poi=ll_poi,
        delta_ll_per_event=delta_ll_per_event,
        brier_h=brier_h, brier_clim=brier_clim, BSS=BSS,
        base_rate=base_rate, n_obs_windows=int(obs.sum()),
    ))
    print(f"  split {sy}: BSS={BSS:+.3f}, ΔlogL/ev={delta_ll_per_event:+.3f}, "
          f"μ0={mu0r:.5f}, c={cr:.2f}d, p={pr:.2f}")

df_roll = pd.DataFrame(roll_rows)
df_roll.to_csv(os.path.join(DATA, "v14_rolling_summary.csv"), index=False)
bss_vals = df_roll.BSS.values
print(f"\n  BSS across 8 splits (v14): median {np.median(bss_vals):+.3f}, "
      f"IQR [{np.quantile(bss_vals,0.25):+.3f}, {np.quantile(bss_vals,0.75):+.3f}], "
      f"range [{bss_vals.min():+.3f}, {bss_vals.max():+.3f}]")
# Reference: v12 and v10
v12_bss = [0.436, 0.438, 0.430, 0.422, 0.432, 0.351, 0.311, 0.406]
v10_bss = [0.412, 0.421, 0.418, 0.395, 0.397, 0.329, 0.349, 0.426]
print(f"  v12 ref: median {np.median(v12_bss):+.3f}, range [{min(v12_bss):+.3f}, {max(v12_bss):+.3f}]")
print(f"  v10 ref: median {np.median(v10_bss):+.3f}, range [{min(v10_bss):+.3f}, {max(v10_bss):+.3f}]")

# ----------------------------------------------------------------------
# 5. Block bootstrap on full fit (B=200, 365-d blocks)
# ----------------------------------------------------------------------
print("\n[5] block bootstrap B=200, block=365d…")
BLOCK = 365.0
NB = int(np.ceil(T_obs / BLOCK))
B = 200
rng_bs = np.random.default_rng(SEED)
params_bs = np.zeros((B, 6))
n_ok = 0
t_bs = time.time()
for b in range(B):
    starts_idx = rng_bs.integers(0, NB, size=NB)
    new_t, new_m = [], []
    cur_T = 0.0
    for kk, si in enumerate(starts_idx):
        b_start = si * BLOCK
        b_end = min(b_start + BLOCK, T_obs)
        sel = (t_all >= b_start) & (t_all < b_end)
        new_t.extend(t_all[sel] - b_start + cur_T)
        new_m.extend(m_all[sel])
        cur_T += (b_end - b_start)
    new_t = np.array(sorted(new_t))
    if len(new_t) < 30:
        params_bs[b] = np.nan
        continue
    # re-sort marks alongside times
    order = np.argsort(new_t)
    new_t = np.array(new_t)[order] if False else new_t  # already sorted
    new_m_arr = np.array(new_m)
    # We sorted only times; align marks: zip then sort. Do it cleanly here.
    pair = sorted(zip(new_t, new_m_arr))
    nt = np.array([p[0] for p in pair])
    nm = np.array([p[1] for p in pair])
    rr = fit_window(nt, nm, float(cur_T), warm=best.x)
    if rr is None or rr.fun > 1e17:
        params_bs[b] = np.nan
        continue
    params_bs[b] = rr.x
    n_ok += 1
    if (b+1) % 50 == 0:
        print(f"  bootstrap {b+1}/{B} ({n_ok} ok, {time.time()-t_bs:.0f}s)")

np.save(os.path.join(DATA, "v14_bootstrap_params.npy"), params_bs)
print(f"  bootstrap complete: {n_ok}/{B} ok")

# Compute CIs (ignoring NaN)
mle_vec = np.array([mu0_h, gamma_h, alpha_h, c_h, p_h, kappa_h])
ci = np.nanpercentile(params_bs, [2.5, 97.5], axis=0)
names = ["mu0", "gamma", "alpha", "c", "p", "kappa"]
print(f"\n   {'param':>8} {'MLE':>12} {'2.5%':>12} {'97.5%':>12}")
for i, nm in enumerate(names):
    print(f"   {nm:>8} {mle_vec[i]:>12.5f} {ci[0,i]:>12.5f} {ci[1,i]:>12.5f}")

# ----------------------------------------------------------------------
# 6. 27-day residual periodogram (the actual test of v12's hypothesis)
# ----------------------------------------------------------------------
print("\n[6] 27-day rotation diagnostic via Lomb-Scargle on residuals…")
# Build per-day intensity for v14
days = np.arange(int(T_obs))
S_d = S_daily[:len(days)]
mu_d = mu0_h * (S_d / S_bar) ** gamma_h
g_all = np.exp(kappa_h * (m_all - m0))
# Excitation at each day d: alpha * sum_{ti<d} g_i * (d-ti+c)^(-p)
exc_d_v14 = np.zeros(len(days))
for i, ti in enumerate(t_all):
    sel = days > ti
    exc_d_v14[sel] += alpha_h * g_all[i] * (days[sel] - ti + c_h) ** (-p_h)
lam_v14 = mu_d + exc_d_v14

# Build v12 intensity for comparison (exp kernel parameters from v12_summary.json)
with open(os.path.join(DATA, "v12_summary.json")) as f:
    v12 = json.load(f)
mu0_v12 = v12["mle"]["mu0"]; gamma_v12 = v12["mle"]["gamma"]
alpha_v12 = v12["mle"]["alpha"]; beta_v12 = v12["mle"]["beta"]; kappa_v12 = v12["mle"]["kappa"]
g_v12 = np.exp(kappa_v12 * (m_all - m0))
mu_d_v12 = mu0_v12 * (S_d / S_bar) ** gamma_v12
exc_d_v12 = np.zeros(len(days))
for i, ti in enumerate(t_all):
    sel = days > ti
    exc_d_v12[sel] += alpha_v12 * g_v12[i] * np.exp(-beta_v12 * (days[sel] - ti))
lam_v12 = mu_d_v12 + exc_d_v12

# Daily count vector
counts = np.zeros(len(days))
for ti in t_all:
    di = int(np.floor(ti))
    if 0 <= di < len(days):
        counts[di] += 1

resid_v14 = counts - lam_v14
resid_v12 = counts - lam_v12

# Periodogram (uniform daily sampling -> use Welch for robust averaging)
freqs_v12, P_v12 = signal.welch(resid_v12, fs=1.0, nperseg=4096)
freqs_v14, P_v14 = signal.welch(resid_v14, fs=1.0, nperseg=4096)
# Find power at 27-day band (24-30 d)
band_lo, band_hi = 1.0/30.0, 1.0/24.0
band_mask_v12 = (freqs_v12 >= band_lo) & (freqs_v12 <= band_hi)
band_mask_v14 = (freqs_v14 >= band_lo) & (freqs_v14 <= band_hi)
peak_v12 = float(np.max(P_v12[band_mask_v12]))
peak_v14 = float(np.max(P_v14[band_mask_v14]))
med_v12 = float(np.median(P_v12))
med_v14 = float(np.median(P_v14))
snr_v12 = peak_v12 / med_v12
snr_v14 = peak_v14 / med_v14
print(f"  27-day band peak power (resid): v12={peak_v12:.4f}, v14={peak_v14:.4f}")
print(f"  background median (resid):     v12={med_v12:.4f}, v14={med_v14:.4f}")
print(f"  signal-to-background ratio:    v12={snr_v12:.2f}, v14={snr_v14:.2f}")
print(f"  → If v14 SNR < v12 SNR, the Omori kernel absorbed the 27-day signal.")

# ----------------------------------------------------------------------
# 7. Plots
# ----------------------------------------------------------------------
print("\n[7] generating plots…")

# F1: Kernel shape comparison (v12 exp vs v14 Omori)
fig, ax = plt.subplots(figsize=(9,5))
tau = np.linspace(0.1, 30, 600)
# normalize both to the same integral over [0, 30] for shape comparison
exp_k = np.exp(-beta_v12 * tau)
omo_k = (tau + c_h) ** (-p_h)
exp_k /= np.trapezoid(exp_k, tau)
omo_k /= np.trapezoid(omo_k, tau)
ax.semilogy(tau, exp_k, 'C0-', lw=2, label=f"v12 exponential, 1/β={1/beta_v12:.2f} d")
ax.semilogy(tau, omo_k, 'C3-', lw=2, label=f"v14 Omori, c={c_h:.2f} d, p={p_h:.2f}")
ax.axvline(27, color="gray", ls=":", alpha=0.6, label="27-day rotation")
ax.set_xlabel("τ since trigger (days)")
ax.set_ylabel("normalised kernel density")
ax.set_title("v14 Omori-Utsu vs v12 exponential excitation kernel\n(same integral over [0, 30 d])")
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "40_v14_kernel_shape.png"), dpi=150)
plt.close()

# F2: Residual periodogram comparison
fig, ax = plt.subplots(figsize=(9,5))
periods_v12 = 1.0/np.where(freqs_v12>0, freqs_v12, np.nan)
periods_v14 = 1.0/np.where(freqs_v14>0, freqs_v14, np.nan)
ax.loglog(periods_v12, P_v12, 'C0-', alpha=0.7, lw=1.2, label="v12 (exp kernel)")
ax.loglog(periods_v14, P_v14, 'C3-', alpha=0.7, lw=1.2, label="v14 (Omori kernel)")
ax.axvspan(24, 30, color="orange", alpha=0.20, label="27-d band")
ax.set_xlim(2, 500)
ax.set_xlabel("period (days)")
ax.set_ylabel("PSD of residual rate")
ax.set_title(f"Residual periodogram: SNR v12={snr_v12:.2f} → v14={snr_v14:.2f} at 27 d")
ax.legend()
ax.grid(alpha=0.3, which="both")
plt.tight_layout()
plt.savefig(os.path.join(FIG, "41_v14_residual_periodogram.png"), dpi=150)
plt.close()

# F3: Rolling BSS comparison (v10, v12, v14)
fig, ax = plt.subplots(figsize=(9,5))
splits_arr = np.array(split_years)
ax.plot(splits_arr, v10_bss, 'C2--o', label="v10 (SSN, exp)")
ax.plot(splits_arr, v12_bss, 'C0--o', label="v12 (F10.7, exp)")
ax.plot(splits_arr, bss_vals, 'C3-o', lw=2, label="v14 (F10.7, Omori)")
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("train/test split year")
ax.set_ylabel("Brier Skill Score (30-day rolling)")
ax.set_title("Out-of-sample BSS: v14 Omori kernel vs v12 exp")
ax.grid(alpha=0.3); ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(FIG, "42_v14_rolling_bss.png"), dpi=150)
plt.close()

# F4: Bootstrap distributions vs MLE
fig, axes = plt.subplots(2, 3, figsize=(13, 7))
v12_refs = {"mu0": v12["mle"]["mu0"], "gamma": v12["mle"]["gamma"],
            "alpha": v12["mle"]["alpha"], "c": None, "p": None,
            "kappa": v12["mle"]["kappa"]}
for ax, name, mle in zip(axes.flat, names, mle_vec):
    vals = params_bs[:, names.index(name)]
    vals = vals[np.isfinite(vals)]
    ax.hist(vals, bins=30, color="C3", alpha=0.65)
    ax.axvline(mle, color="k", lw=2, label=f"v14 MLE = {mle:.4f}")
    ref = v12_refs[name]
    if ref is not None:
        ax.axvline(ref, color="C0", lw=2, ls=":", label=f"v12 ref = {ref:.4f}")
    ax.set_title(name)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "43_v14_bootstrap.png"), dpi=150)
plt.close()

print("  saved figures 40-43")

# ----------------------------------------------------------------------
# 8. Summary JSON
# ----------------------------------------------------------------------
summary = {
    "seed": SEED,
    "window": {"start": str(T0.date()), "end": str(T_END.date()),
               "T_obs_days": T_obs, "T_obs_yr": T_obs/365.25},
    "S_bar": S_bar,
    "n_events": int(N),
    "kernel": "Omori-Utsu (t-ti+c)^(-p)",
    "mle": {
        "mu0": float(mu0_h), "gamma": float(gamma_h),
        "alpha": float(alpha_h), "c": float(c_h), "p": float(p_h),
        "kappa": float(kappa_h),
        "half_amplitude_d": float(half_decay),
        "half_integral_d": float(s_med) if np.isfinite(s_med) else None,
        "branching_ratio": float(branching_v14) if np.isfinite(branching_v14) else None,
        "logL": float(ll_full), "AIC": float(AIC),
    },
    "vs_v12": {
        "delta_logL": float(ll_full - ll_v12),
        "delta_AIC": float(AIC - AIC_v12),
    },
    "rolling_oos": {
        "split_years": split_years,
        "BSS_v14": [float(x) for x in bss_vals],
        "BSS_v12_ref": v12_bss,
        "BSS_v10_ref": v10_bss,
        "median_BSS_v14": float(np.median(bss_vals)),
        "median_BSS_v12": float(np.median(v12_bss)),
        "median_BSS_v10": float(np.median(v10_bss)),
    },
    "bootstrap_CIs": {nm: [float(ci[0,i]), float(ci[1,i])]
                      for i, nm in enumerate(names)},
    "rotation_diagnostic": {
        "peak_v12": peak_v12, "peak_v14": peak_v14,
        "median_v12": med_v12, "median_v14": med_v14,
        "snr_v12": snr_v12, "snr_v14": snr_v14,
        "delta_snr": snr_v14 - snr_v12,
    },
}
with open(os.path.join(DATA, "v14_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)

print("\n[done] summary written to data/v14_summary.json")
print(f"  ΔlogL vs v12 = {ll_full - ll_v12:+.2f}")
print(f"  ΔAIC vs v12  = {AIC - AIC_v12:+.2f}")
print(f"  v14 BSS median = {np.median(bss_vals):+.3f}  (v12 was {np.median(v12_bss):+.3f})")
print(f"  27d SNR  v12 → v14: {snr_v12:.2f} → {snr_v14:.2f}")
