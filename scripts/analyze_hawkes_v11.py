#!/usr/bin/env python3
"""
v11: Pre-1844 extension via Helsinki magnetic observatory (Nevanlinna).

We add the Helsinki H-component K/Ak record (1 July 1844 - 31 May 1897) to
the front of v7's event series. This is the Lockwood/Nevanlinna composite
data used to construct the homogeneous IDV(1d) index. It brings the
Carrington Event of 1-2 September 1859 -- the largest geomagnetic storm in
the recorded history of Earth -- into the marked Hawkes fit.

Pipeline
--------
  1. Parse Helsinki Ak file (fixed-width). Mask missing (9999 or **** K).
  2. Calibrate Helsinki Ak threshold against aa-index G4+ proxy
     (aa_max >= 212) in the 1868-1880 overlap window.
  3. Detect Carrington and pre-1868 events; assign Kp-equivalent marks
     using a monotone Ak->Kp map (Ak=400 -> Kp=9.5 special: extreme).
  4. Build extended event series 1844-2025 (~181 yr) and refit v7 marked
     Hawkes on the full window. Evaluate:
       - Parameter stability vs v7 (1868-2025)
       - Tail-conditional logL: how surprised is the model by Carrington?
       - Recurrence probability for an Ak>=400 event under the fit
  5. Block-bootstrap CIs.

Random seed: 20260523 throughout.
"""

import os, time, json, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize, stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FIG  = os.path.join(ROOT, "figures")

rng = np.random.default_rng(20260523)

# ----------------------------------------------------------------------
# 1. Parse the Helsinki H-component K/Ak file
#    Format (fixed-width):
#       cols 1-4: year, 5-6: month, 7-8: day  (compact YYYYMMDD)
#       col 9: space
#       cols 10-17: K-values (8 digits, '*' = missing)
#       cols 18-19: spaces
#       cols 20-22: daily Ak (right-justified, missing = 9999)
# ----------------------------------------------------------------------
print("[helsinki] parsing 1844-1897 K/Ak record…")
rows = []
with open(os.path.join(DATA, "helsinki_H_K_1844-1897.txt")) as f:
    for line in f:
        line = line.rstrip("\n")
        if len(line) < 22 or not line[:8].isdigit():
            continue
        date_str = line[:8]
        try:
            d = pd.Timestamp(year=int(date_str[:4]),
                             month=int(date_str[4:6]),
                             day=int(date_str[6:8]))
        except ValueError:
            continue
        k_str = line[9:17]   # 8 chars
        ak_str = line[18:].strip()  # take everything after col 17 (incl. trailing 9999)
        # Count missing K values
        n_missing = sum(1 for c in k_str if c == "*")
        # max K (treat '*' as missing for non-extremes)
        k_digits = [int(c) for c in k_str if c.isdigit()]
        kmax = max(k_digits) if k_digits else np.nan
        try:
            ak = int(ak_str)
        except ValueError:
            ak = np.nan
        # Helsinki encoding: 9999 = whole-day missing, 999 = unmeasurable.
        # On extreme storm days some K's are written as '*' when the
        # magnetometer pegged out -- those K=* should be treated as K>=9.
        # We only drop the Ak when the file flags it explicitly (9999/999)
        # OR when ALL 8 K's are missing AND the daily Ak is also missing.
        if ak == 9999 or ak == 999:
            ak = np.nan
        if n_missing == 8 and np.isnan(ak):
            ak = np.nan  # no info
        rows.append((d, kmax, ak, n_missing))
df_h = pd.DataFrame(rows, columns=["date", "k_max_helsinki", "ak_helsinki", "n_missing_k"])
df_h = df_h.sort_values("date").reset_index(drop=True)
n_total = len(df_h)
n_valid = df_h.ak_helsinki.notna().sum()
print(f"  rows: {n_total}, valid Ak: {n_valid}  ({n_valid/n_total*100:.1f}%)")
print(f"  date range: {df_h.date.min().date()} → {df_h.date.max().date()}")
print(f"  Ak range (valid only): {df_h.ak_helsinki.min():.0f} → {df_h.ak_helsinki.max():.0f}")
print(f"  K_max range: {int(df_h.k_max_helsinki.min())} → {int(df_h.k_max_helsinki.max())}")

# Spotlight Carrington
print("\n  Sept 1859 (Carrington Event):")
print(df_h[(df_h.date >= "1859-08-28") & (df_h.date <= "1859-09-10")].to_string(index=False))

df_h.to_csv(os.path.join(DATA, "derived_helsinki_daily.csv"), index=False)

# ----------------------------------------------------------------------
# 2. Calibrate Helsinki Ak vs aa-index in 1868-1880 overlap
# ----------------------------------------------------------------------
aa = pd.read_csv(os.path.join(DATA, "derived_aa_daily.csv"), parse_dates=["date"])
overlap = pd.merge(df_h, aa[["date","aa_max"]], on="date", how="inner").dropna()
overlap = overlap[(overlap.date >= "1868-01-01") & (overlap.date <= "1880-06-30")]
print(f"\n[overlap] Helsinki × aa 1868-1880: {len(overlap)} days")
print(f"  Helsinki Ak vs aa_max correlation: r = {overlap.ak_helsinki.corr(overlap.aa_max):.3f}")

# Pick Helsinki threshold matching G4+ event rate from aa
# v7 calibration: aa_max >= 212 = G4+
aa_g4 = (overlap.aa_max >= 212).sum()
print(f"  aa-defined G4+ events in overlap: {aa_g4}")
candidate_thr = sorted(overlap.ak_helsinki.dropna().unique(), reverse=True)
# Find the Helsinki Ak threshold that gives ~same event count
target = aa_g4
best_thr = None; best_diff = 1e9
for thr in range(40, 250):
    n = (overlap.ak_helsinki >= thr).sum()
    diff = abs(n - target)
    if diff < best_diff:
        best_diff = diff; best_thr = thr
print(f"  Helsinki Ak >= {best_thr} gives {(overlap.ak_helsinki >= best_thr).sum()} days "
      f"(target {target})")

# Day-level agreement
helsinki_flag = overlap.ak_helsinki >= best_thr
aa_flag = overlap.aa_max >= 212
both = (helsinki_flag & aa_flag).sum()
agree = ((helsinki_flag & aa_flag) | (~helsinki_flag & ~aa_flag)).mean()
print(f"  day-level both-flag agreement: {agree*100:.0f}% ({both} co-occurring G4+)")

# Helsinki Ak >= 400 = "extreme" (Carrington-class) — special handling for G5+ tier
# Calibrate Kp-equivalent mark using Ak: rough monotone map
# Ak ≈ 50-100  → Kp 7-7.67  (G3)
# Ak ≈ 100-180 → Kp 7.67-8.33 (G4)
# Ak ≈ 180-300 → Kp 8.33-9.0 (G4/G5 boundary)
# Ak ≈ 300+    → Kp 9.0+    (G5)
# Ak ≈ 400+    → Kp 9.5     (extreme — assign as Carrington-class)
def ak_to_kp(ak):
    if pd.isna(ak): return np.nan
    if ak >= 400: return 9.5     # Carrington / extreme
    if ak >= 250: return 9.0     # G5
    if ak >= 180: return 8.667   # high G4
    if ak >= 130: return 8.333   # G4
    if ak >= best_thr: return 8.0  # G4 low
    return np.nan

# ----------------------------------------------------------------------
# 3. Build pre-1868 event series from Helsinki
# ----------------------------------------------------------------------
pre_1868 = df_h[(df_h.date < "1868-01-01") & (df_h.ak_helsinki >= best_thr)].copy()
pre_1868["mark"] = pre_1868.ak_helsinki.apply(ak_to_kp)
pre_1868 = pre_1868.dropna(subset=["mark"])
print(f"\n[pre-1868] {len(pre_1868)} events found at Ak >= {best_thr}")
print(f"  mark distribution:")
print(pre_1868.mark.value_counts().sort_index().to_string())
print(f"  largest events (top 10):")
print(pre_1868.nlargest(10, "ak_helsinki")[["date","ak_helsinki","mark"]].to_string(index=False))

# Carrington must be in
assert (pre_1868.date == "1859-09-03").any(), "Carrington missing!"
carrington_row = pre_1868[pre_1868.date == "1859-09-03"].iloc[0]
print(f"\n  ★ Carrington Event captured: {carrington_row.date.date()}, "
      f"Ak={carrington_row.ak_helsinki:.0f}, mark={carrington_row.mark}")

pre_1868_events = pre_1868[["date", "mark"]].copy()
pre_1868_events["source"] = "helsinki"

# Load existing v7 event series (1868-2025)
existing = pd.read_csv(os.path.join(DATA, "derived_events_extended_1868_2025.csv"),
                       parse_dates=["date"])
print(f"\n[v7 events] {len(existing)} events 1868-2025")

# Merge
all_events = pd.concat([pre_1868_events, existing], ignore_index=True)
all_events = all_events.sort_values("date").reset_index(drop=True)
print(f"[v11 events] total {len(all_events)} events, "
      f"{all_events.date.min().date()} → {all_events.date.max().date()}")
all_events.to_csv(os.path.join(DATA, "derived_events_extended_1844_2025.csv"), index=False)

# ----------------------------------------------------------------------
# 4. SSN — we have SILSO back to 1818, plenty enough for 1844
# ----------------------------------------------------------------------
T0 = pd.Timestamp("1844-07-01")
T_END = existing.date.max()  # use Kp daily max as before
T_END_days = (T_END - T0).days
print(f"\n[window] T0 = {T0.date()}, T_end = {T_END.date()}, "
      f"T = {T_END_days/365.25:.1f} yr")

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

mask_w = (ssn_t >= 0) & (ssn_t <= T_END_days)
S_bar = ssn_s[mask_w].mean()
print(f"  S_bar (1844-2025) = {S_bar:.2f}")

# ----------------------------------------------------------------------
# 5. Refit v7 marked Hawkes on full 1844-2025 series
# ----------------------------------------------------------------------
m0 = 8.0
def loglike(params, t, m, T_obs, t_start, S_bar):
    mu0, gamma, alpha, beta, kappa = params
    if mu0 <= 0 or alpha < 0 or beta <= 0: return -1e18
    S_e = S(t)
    mu_e = mu0 * (S_e / S_bar) ** gamma
    g = np.exp(kappa * (m - m0))
    R = 0.0; s_log = 0.0
    for i in range(len(t)):
        if i > 0:
            R = np.exp(-beta * (t[i] - t[i-1])) * (R + g[i-1])
        rate = mu_e[i] + alpha * R
        if rate <= 0: return -1e18
        s_log += np.log(rate)
    grid = np.arange(t_start, t_start + T_obs + 1.0, 1.0)
    Sg = S(grid)
    mu_grid = mu0 * (Sg / S_bar) ** gamma
    s_int_mu = np.trapezoid(mu_grid, grid)
    T_end = t_start + T_obs
    s_comp = (alpha/beta) * np.sum(g * (1.0 - np.exp(-beta*(T_end - t))))
    return s_log - s_int_mu - s_comp

t_int = (all_events.date - T0).dt.days.values.astype(float)
jitter = rng.uniform(0.0, 1.0, size=len(t_int))
order = np.argsort(t_int + jitter)
t_all = (t_int + jitter)[order]
m_all = all_events.mark.values[order].astype(float)
print(f"\n[fit] N = {len(t_all)} events over {T_END_days:.0f}d = {T_END_days/365.25:.1f}yr")
print(f"  marks: G4 = {(m_all<9).sum()}, G5 = {((m_all>=9)&(m_all<9.5)).sum()}, "
      f"Carrington-class (m=9.5): {(m_all>=9.5).sum()}")

starts = [
    (0.00443, 0.995, 0.1133, 0.6424, 0.899),
    (0.005, 0.85, 0.12, 0.65, 1.00),
    (0.004, 1.10, 0.10, 0.60, 0.80),
    (0.006, 0.70, 0.15, 0.70, 0.50),
    (0.0045, 1.00, 0.11, 0.64, 0.90),
    (0.005, 0.90, 0.13, 0.55, 1.20),
]
best = None
t_fit = time.time()
print("\n[mle] multi-start MLE on 1844-2025…")
for k, x0 in enumerate(starts):
    res = optimize.minimize(lambda p: -loglike(p, t_all, m_all, float(T_END_days), 0.0, S_bar),
                            x0, method="Nelder-Mead",
                            options={"xatol":1e-7,"fatol":1e-7,"maxiter":80000})
    if res.fun < 1e17:
        mu0r,gr,ar,br,kr = res.x
        print(f"  trial {k+1}: μ0={mu0r:.5f} γ={gr:.3f} α={ar:.4f} β={br:.4f} κ={kr:+.3f}  -LL={res.fun:.2f}")
        if best is None or res.fun < best.fun:
            best = res
mu0_h, gamma_h, alpha_h, beta_h, kappa_h = best.x
ll_full = -best.fun
print(f"\n[MLE v11 on 1844-2025]")
print(f"  μ0    = {mu0_h:.5f}/d = {mu0_h*365.25:.3f}/yr at S̄")
print(f"  γ     = {gamma_h:.4f}")
print(f"  α     = {alpha_h:.4f}")
print(f"  β     = {beta_h:.4f}  → 1/β = {1/beta_h:.2f} d")
print(f"  κ     = {kappa_h:+.4f}  → exp(κ) = {np.exp(kappa_h):.3f}×")
print(f"  log-L = {ll_full:.2f},  AIC = {2*5 - 2*ll_full:.2f}")
print(f"  η(G4) = {alpha_h/beta_h:.3f}, η(G5) = {alpha_h*np.exp(kappa_h)/beta_h:.3f}, "
      f"η(Carrington-class m=9.5) = {alpha_h*np.exp(kappa_h*1.5)/beta_h:.3f}")
print(f"  fit time {time.time()-t_fit:.1f}s")

# ----------------------------------------------------------------------
# 6. Tail-conditional analysis: how surprised is the model by Carrington?
# ----------------------------------------------------------------------
print("\n[tail] Probability that an extreme G5 occurs in any given year")
# Per-year expected number of background G5 events at the mean SSN level
mu_g5_yr = mu0_h * 365.25  # at S̄, the background rate of G4+
# Per-event probability of being marked as G5+ requires the mark distribution
# Compute empirical mark distribution above G4 threshold
n_carr = (m_all >= 9.5).sum()
n_g5 = ((m_all >= 9.0) & (m_all < 9.5)).sum()
n_g4_lo = (m_all < 9.0).sum()
print(f"  empirical marks: G4 = {n_g4_lo}, G5 (m=9.0) = {n_g5}, "
      f"Carrington-class (m=9.5) = {n_carr}")
print(f"  empirical fraction Carrington-class: {n_carr/len(m_all)*100:.2f}%")
print(f"  empirical recurrence (mean inter-event time): "
      f"{T_END_days/(365.25*n_carr):.1f} yr per Carrington-class event")

# Likelihood contribution of the Carrington event
# Locate Carrington in the sorted array
carr_idx = np.argmin(np.abs(t_all - (pd.Timestamp("1859-09-03") - T0).days))
carr_t = t_all[carr_idx]; carr_m = m_all[carr_idx]
S_e_all = S(t_all)
mu_e_all = mu0_h * (S_e_all / S_bar) ** gamma_h
# rate at Carrington from preceding history
R_carr = 0.0
for j in range(carr_idx):
    R_carr += np.exp(kappa_h*(m_all[j] - m0)) * np.exp(-beta_h*(carr_t - t_all[j]))
rate_at_carr = mu_e_all[carr_idx] + alpha_h * R_carr
print(f"\n[Carrington context]")
print(f"  intensity λ*(t) at Carrington 1859-09-03 = {rate_at_carr:.4f}/day = "
      f"{rate_at_carr*365.25:.2f}/yr")
print(f"  background contribution: {mu_e_all[carr_idx]:.4f}/day  "
      f"({mu_e_all[carr_idx]/rate_at_carr*100:.0f}% of total)")
print(f"  excitation contribution: {alpha_h*R_carr:.4f}/day  "
      f"({alpha_h*R_carr/rate_at_carr*100:.0f}% of total)")
print(f"  log-density of next event time conditional on history: "
      f"{np.log(rate_at_carr):.3f}")

# Per-event marginal log-likelihood at Carrington
# (compared to the mean across all events)
all_lls = []
R = 0.0
for i in range(len(t_all)):
    if i > 0:
        R = np.exp(-beta_h * (t_all[i]-t_all[i-1])) * (R + np.exp(kappa_h*(m_all[i-1]-m0)))
    r = mu_e_all[i] + alpha_h*R
    if r > 0:
        all_lls.append(np.log(r))
all_lls = np.array(all_lls)
print(f"  Carrington event-rank by log-density: "
      f"{(all_lls < all_lls[carr_idx]).sum() + 1} of {len(all_lls)} (higher rank = less surprising)")
print(f"  Carrington log-density quantile: {(all_lls < all_lls[carr_idx]).mean()*100:.0f}th percentile")
print(f"  (50th percentile would mean Carrington was a typical-looking event;")
print(f"   high percentile means the model thought it was *more* likely than average)")

# Predicted recurrence time for an Ak>=400 event
# Use: rate of m>=9.5 events ≈ N(m>=9.5)/T_total
rate_extreme = n_carr / (T_END_days/365.25)  # per year
print(f"\n  Implied recurrence rate for Carrington-class (m=9.5): "
      f"{1/rate_extreme:.0f} yr (empirical, 1 event in {T_END_days/365.25:.0f}yr)")

# ----------------------------------------------------------------------
# 7. Block bootstrap (B=200)
# ----------------------------------------------------------------------
print("\n[bootstrap] block bootstrap B=200, block=365d…")
BLOCK = 365.0
T_all = float(T_END_days)
NB = int(np.ceil(T_all/BLOCK))
B = 200
rng_bs = np.random.default_rng(20260523)
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
    if not new_t: continue
    tb = np.concatenate(new_t); mb = np.concatenate(new_m)
    ord_b = np.argsort(tb)
    tb = tb[ord_b]; mb = mb[ord_b]
    try:
        res_b = optimize.minimize(lambda p: -loglike(p, tb, mb, cur_T, 0.0, S_bar),
                                  best.x, method="Nelder-Mead",
                                  options={"xatol":1e-5,"fatol":1e-5,"maxiter":25000})
        if res_b.fun < 1e17 and np.all(np.isfinite(res_b.x)):
            params_bs[n_ok] = res_b.x
            n_ok += 1
    except Exception:
        pass
    if (b+1) % 50 == 0:
        print(f"  bootstrap {b+1}/{B}  ({n_ok} ok)  elapsed {time.time()-t_bs:.0f}s")
params_bs = params_bs[:n_ok]
print(f"\n[bootstrap] {n_ok}/{B} replicates ({time.time()-t_bs:.0f}s)")
names_p = ["μ0","γ","α","β","κ"]
hat = list(best.x)
ci = {}
print(f"\n  {'param':>4}  {'MLE':>10}  {'2.5%':>10}  {'97.5%':>10}")
for i, name in enumerate(names_p):
    q = np.quantile(params_bs[:, i], [0.025, 0.5, 0.975])
    ci[name] = (float(q[0]), float(q[2]))
    print(f"  {name:>4}  {hat[i]:10.5f}  {q[0]:10.5f}  {q[2]:10.5f}")

# Derived CIs
inv_beta_bs = 1/params_bs[:,3]
exp_kappa_bs = np.exp(params_bs[:,4])
print(f"  1/β   {1/best.x[3]:10.3f}  {np.quantile(inv_beta_bs,0.025):10.3f}  "
      f"{np.quantile(inv_beta_bs,0.975):10.3f}")
print(f"  exp(κ){np.exp(best.x[4]):10.3f}  {np.quantile(exp_kappa_bs,0.025):10.3f}  "
      f"{np.quantile(exp_kappa_bs,0.975):10.3f}")
np.save(os.path.join(DATA, "v11_bootstrap_params.npy"), params_bs)

# ----------------------------------------------------------------------
# 8. Save and plot
# ----------------------------------------------------------------------
summary = {
    "model": "v11 marked Hawkes — 1844-2025 (Helsinki + aa + Kp)",
    "total_events": int(len(t_all)),
    "T_years": float(T_END_days/365.25),
    "pre_1868_events": int(len(pre_1868_events)),
    "helsinki_threshold_Ak": int(best_thr),
    "carrington": {
        "date": "1859-09-03",
        "Ak_helsinki": float(carrington_row.ak_helsinki),
        "mark_assigned": float(carrington_row.mark),
        "lambda_at_event_per_day": float(rate_at_carr),
        "lambda_at_event_per_year": float(rate_at_carr*365.25),
        "background_share": float(mu_e_all[carr_idx]/rate_at_carr),
        "excitation_share": float(alpha_h*R_carr/rate_at_carr),
        "log_density_percentile": float((all_lls < all_lls[carr_idx]).mean()),
    },
    "mle_v11": {
        "mu0": float(mu0_h), "gamma": float(gamma_h), "alpha": float(alpha_h),
        "beta": float(beta_h), "kappa": float(kappa_h),
        "inv_beta": float(1/beta_h), "exp_kappa": float(np.exp(kappa_h)),
        "mu0_per_yr": float(mu0_h*365.25),
        "logL": float(ll_full), "AIC": float(2*5 - 2*ll_full),
        "eta_g4": float(alpha_h/beta_h),
        "eta_g5": float(alpha_h*np.exp(kappa_h)/beta_h),
        "eta_carr": float(alpha_h*np.exp(kappa_h*1.5)/beta_h),
    },
    "v7_reference": {  # what we got at 1868-2025
        "mu0_per_yr": 1.62, "gamma": 0.995,
        "inv_beta": 1.56, "exp_kappa": 2.46,
    },
    "bootstrap_ci_95": {k: list(v) for k,v in ci.items()},
}
with open(os.path.join(DATA, "v11_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print(f"\n[save] data/v11_summary.json")

# Figure 30: extended event series (tick plot)
fig, ax = plt.subplots(figsize=(14, 4.5))
for tt, mm in zip(t_all, m_all):
    if mm >= 9.5:
        c, h = "darkred", 1.0
    elif mm >= 9.0:
        c, h = "C3", 0.85
    else:
        c, h = "C0", 0.5
    ax.vlines(T0 + pd.to_timedelta(tt, unit="D"), 0, h, colors=c, lw=1.0, alpha=0.7)
# annotations
events_to_label = [
    ("1859-09-03", "Carrington 1859", 1.0, "darkred"),
    ("1872-02-04", "1872", 0.55, "C3"),
    ("1882-11-17", "Stewart 1882", 0.55, "C3"),
    ("1921-05-15", "NY Railroad 1921", 0.55, "C3"),
    ("1989-03-13", "Hydro-Québec 1989", 0.55, "C3"),
    ("2003-10-29", "Halloween 2003", 0.55, "C3"),
    ("2024-05-10", "Gannon 2024", 0.55, "C3"),
]
for date_str, label, h, c in events_to_label:
    try:
        d = pd.Timestamp(date_str)
        if T0 <= d <= T_END:
            ax.annotate(label, (d, h), xytext=(0, 8), textcoords="offset points",
                        fontsize=8, ha="center", color=c, alpha=0.9,
                        arrowprops=dict(arrowstyle="-", color=c, alpha=0.4, lw=0.5))
    except: pass
ax.set_xlim(T0, T_END); ax.set_ylim(0, 1.45); ax.set_yticks([])
ax.set_xlabel("Date")
ax.set_title(f"v11 — Geomagnetic G4+ event series, 1844-2025  "
             f"(Helsinki + aa + Kp,  N={len(t_all)} events over "
             f"{T_END_days/365.25:.0f} yr)")
plt.tight_layout()
plt.savefig(os.path.join(FIG, "30_v11_events_1844_2025.png"), dpi=140)
plt.close()
print("[plot] figures/30_v11_events_1844_2025.png")

# Figure 31: log-density of each event under v11; Carrington highlighted
fig, ax = plt.subplots(figsize=(11, 5))
dates_all = T0 + pd.to_timedelta(t_all, unit="D")
ax.scatter(dates_all, all_lls, s=20, c="C0", alpha=0.5, label="all events")
# Carrington
ax.scatter(dates_all[carr_idx], all_lls[carr_idx], s=200, c="darkred",
           marker="*", zorder=5, label=f"Carrington 1859 (log-λ = {all_lls[carr_idx]:.2f})")
ax.axhline(all_lls.mean(), color="0.5", ls="--", alpha=0.7,
           label=f"mean log-density = {all_lls.mean():.2f}")
ax.set_ylabel("log λ*(t) at each event\n(higher = model finds it more expected)")
ax.set_xlabel("Date")
ax.set_title(f"v11 — Event-by-event log conditional intensity\n"
             f"Carrington is at {(all_lls < all_lls[carr_idx]).mean()*100:.0f}th percentile of log-density")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "31_v11_carrington_logdensity.png"), dpi=140)
plt.close()
print("[plot] figures/31_v11_carrington_logdensity.png")

# Figure 32: Parameter overlay v7 vs v11
fig, axes = plt.subplots(1, 4, figsize=(15, 4))
labels_p = [("μ0 (per yr)", "mu0_per_yr"), ("γ", "gamma"),
            ("1/β (days)", "inv_beta"), ("exp(κ)", "exp_kappa")]
v7_vals = {"mu0_per_yr": 1.62, "gamma": 0.995, "inv_beta": 1.56, "exp_kappa": 2.46}
v11_vals = {"mu0_per_yr": mu0_h*365.25, "gamma": gamma_h,
            "inv_beta": 1/beta_h, "exp_kappa": np.exp(kappa_h)}
bs_arrays = {
    "mu0_per_yr": params_bs[:,0]*365.25,
    "gamma": params_bs[:,1],
    "inv_beta": 1/params_bs[:,3],
    "exp_kappa": np.exp(params_bs[:,4]),
}
for ax, (label, key) in zip(axes, labels_p):
    ax.hist(bs_arrays[key], bins=20, color="C0", alpha=0.6, edgecolor="white")
    ax.axvline(v11_vals[key], color="k", lw=2, label=f"v11 = {v11_vals[key]:.3f}")
    ax.axvline(v7_vals[key], color="C3", lw=2, ls="--", label=f"v7 = {v7_vals[key]:.3f}")
    ax.set_xlabel(label); ax.set_ylabel("count")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
plt.suptitle(f"v11 vs v7 parameters  (n_bootstrap = {n_ok})",
             fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "32_v11_v7_param_overlay.png"), dpi=140)
plt.close()
print("[plot] figures/32_v11_v7_param_overlay.png")

print("\n[done] v11 complete.")
