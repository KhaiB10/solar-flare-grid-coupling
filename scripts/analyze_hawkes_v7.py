#!/usr/bin/env python3
"""
v7: Pre-1932 extension via the aa geomagnetic index (1868-1932) bolted onto
the GFZ Kp record (1932-2025). 158-year extended marked Hawkes with bootstrap
confidence intervals.

Strategy:
  1. Daily aa_max (peak of 8 three-hour values) from NCEI 1868-2010.
  2. Calibrate aa_max thresholds in the 1932-2010 overlap so the resulting
     event rate matches the GFZ Kp G4+ and G5 rates over the same window
     (rate-matching, not classification — see FINDINGS_v7 for rationale).
  3. For pre-1932 events, assign a binary "G4" or "G5" mark using the two
     calibrated thresholds.  For 1932+ events, keep the original Kp_max mark.
  4. Refit the v6 marked Hawkes on the merged 1868-2025 series.
  5. Block-bootstrap 95% CIs for all parameters (block length 365 days to
     preserve within-cycle dependence).
"""

import os, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize, stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FIG  = os.path.join(ROOT, "figures")

# ----------------------------------------------------------------------
# 1. Load aa daily series + Kp daily series
# ----------------------------------------------------------------------
aa = pd.read_csv(os.path.join(DATA, "derived_aa_daily.csv"), parse_dates=["date"])
kp = pd.read_csv(os.path.join(DATA, "derived_daily_with_phase.csv"), parse_dates=["date"])
print(f"[load] aa: {len(aa):,} days, {aa.date.min().date()} → {aa.date.max().date()}")
print(f"[load] kp: {len(kp):,} days, {kp.date.min().date()} → {kp.date.max().date()}")

# ----------------------------------------------------------------------
# 2. Rate-matched aa thresholds from 1932-2010 overlap
# ----------------------------------------------------------------------
merge = aa.merge(kp[["date","Kp_max"]], on="date", how="inner")
n_g4 = int((merge.Kp_max >= 8).sum())
n_g5 = int((merge.Kp_max >= 9).sum())
aa_sorted = np.sort(merge.aa_max.values)[::-1]
THR_G4 = float(aa_sorted[n_g4 - 1])
THR_G5 = float(aa_sorted[n_g5 - 1])
agree = ((merge.aa_max >= THR_G4) & (merge.Kp_max >= 8)).sum()
print(f"[cal] Overlap: {len(merge):,} days, {n_g4} G4+ (Kp), {n_g5} G5 (Kp)")
print(f"[cal] Rate-matched aa_max thresholds: G4+ ≥ {THR_G4:.0f},  G5 ≥ {THR_G5:.0f}")
print(f"[cal] Day-level agreement on G4+ identification: {agree}/{n_g4} = {agree/n_g4*100:.1f}%")

# ----------------------------------------------------------------------
# 3. Build merged 1868-2025 event series
# ----------------------------------------------------------------------
# Pre-1932: use aa thresholds, mark = 9 if aa_max ≥ THR_G5 else 8
pre = aa[aa.date < pd.Timestamp("1932-01-01")].copy()
pre_events = pre[pre.aa_max >= THR_G4].copy()
pre_events["mark"] = np.where(pre_events.aa_max >= THR_G5, 9.0, 8.0)
pre_events["source"] = "aa"
pre_events = pre_events[["date", "mark", "source"]]

# 1932+: use Kp_max
post = kp[kp.Kp_max >= 8][["date", "Kp_max"]].copy().rename(columns={"Kp_max":"mark"})
post["source"] = "kp"

events = pd.concat([pre_events, post], ignore_index=True).sort_values("date").reset_index(drop=True)
print(f"\n[merge] Total events 1868-2025: {len(events)}")
print(f"        Pre-1932 (aa-derived):  {(events.source=='aa').sum()}  "
      f"({(events.source=='aa').sum() - (events[events.source=='aa'].mark>=9).sum()} G4, "
      f"{(events[events.source=='aa'].mark>=9).sum()} G5)")
print(f"        1932+ (Kp-derived):     {(events.source=='kp').sum()}  "
      f"({(events[events.source=='kp'].mark<9).sum()} G4-class, "
      f"{(events[events.source=='kp'].mark>=9).sum()} G5)")

# Time origin = 1868-01-01
T0 = pd.Timestamp("1868-01-01")
T_end = (aa.date.max() if aa.date.max() > kp.date.max() else kp.date.max())
T_end_days = (kp.date.max() - T0).days  # use Kp upper bound so series is consistent through 2025
print(f"        Observation window: {T_end_days:,} days = {T_end_days/365.25:.1f} years")
print(f"        Long-run G4+ rate: {len(events)/(T_end_days/365.25):.2f} events/yr  "
      f"(v6 used {246/94.0:.2f}/yr over 1932-2025)")

# Jitter ties and sort
rng = np.random.default_rng(20260523)
t_int = (events.date - T0).dt.days.values.astype(float)
jitter = rng.uniform(0.0, 1.0, size=len(t_int))
order = np.argsort(t_int + jitter)
t = (t_int + jitter)[order]
m = events.mark.values[order].astype(float)
source = events.source.values[order]
N = len(t)
T = float(T_end_days)
print(f"        N = {N}, T = {T:.0f} days")

# Save merged events
pd.DataFrame({"date": events.date.values, "mark": events.mark.values, "source": events.source.values}).to_csv(
    os.path.join(DATA, "derived_events_extended_1868_2025.csv"), index=False)

# ----------------------------------------------------------------------
# 4. SSN modulation (extended back to 1868)
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
mask_w = (ssn_t >= 0) & (ssn_t <= T)
S_bar = ssn_s[mask_w].mean()
print(f"[ssn] Window-mean smoothed SSN 1868-2025 = {S_bar:.2f}  "
      f"(v6 used 93.50 for 1932-2025)")

def S(td):
    return np.interp(td, ssn_t, ssn_s, left=ssn_s[0], right=ssn_s[-1])

S_events = S(t)
grid = np.arange(0, T + 1.0, 1.0)
S_grid = S(grid)

# ----------------------------------------------------------------------
# 5. Marked Hawkes log-likelihood (same form as v6)
# ----------------------------------------------------------------------
m0 = 8.0
def loglike(params, t, m, T, S_events, grid, S_grid, S_bar, m0):
    mu0, gamma, alpha, beta, kappa = params
    if mu0 <= 0 or alpha < 0 or beta <= 0: return -1e18
    mu_events = mu0 * (S_events / S_bar) ** gamma
    g = np.exp(kappa * (m - m0))
    R = 0.0; s_log = 0.0
    for i in range(len(t)):
        if i > 0:
            R = np.exp(-beta * (t[i] - t[i-1])) * (R + g[i-1])
        rate = mu_events[i] + alpha * R
        if rate <= 0: return -1e18
        s_log += np.log(rate)
    mu_grid = mu0 * (S_grid / S_bar) ** gamma
    s_int_mu = np.trapezoid(mu_grid, grid)
    s_comp = (alpha/beta) * np.sum(g * (1.0 - np.exp(-beta*(T - t))))
    return s_log - s_int_mu - s_comp

def neg_ll(p, *a): return -loglike(p, *a)

# ----------------------------------------------------------------------
# 6. MLE fit
# ----------------------------------------------------------------------
print("\n[fit] multi-start MLE on extended 1868-2025 series…")
starts = [
    (0.00549, 0.845, 0.111, 0.652, 1.00),   # v6 best
    (0.005,   0.7,   0.10,  0.50,  0.5),
    (0.006,   1.0,   0.15,  0.70,  1.2),
    (0.004,   1.2,   0.08,  0.30,  0.8),
    (0.007,   0.5,   0.18,  0.80,  0.3),
    (0.005,   0.9,   0.12,  0.60,  1.5),
]
best = None
t0 = time.time()
for trial, x0 in enumerate(starts):
    res = optimize.minimize(neg_ll, x0,
                            args=(t, m, T, S_events, grid, S_grid, S_bar, m0),
                            method="Nelder-Mead",
                            options={"xatol":1e-7,"fatol":1e-7,"maxiter":60000})
    if res.fun < 1e17:
        mu0r,gr,ar,br,kr = res.x
        print(f"  trial {trial+1}: μ0={mu0r:.5f} γ={gr:.3f} α={ar:.4f} β={br:.4f} κ={kr:+.3f}  -LL={res.fun:.2f}")
        if best is None or res.fun < best.fun:
            best = res
print(f"  fit time: {time.time()-t0:.1f}s")

mu0_h, gamma_h, alpha_h, beta_h, kappa_h = best.x
eta_g4 = alpha_h / beta_h
eta_g5 = alpha_h * np.exp(kappa_h) / beta_h
ll_v7 = -best.fun
aic_v7 = 2*5 - 2*ll_v7
print(f"\n[MLE v7]:")
print(f"  μ0 = {mu0_h:.5f} /day = {mu0_h*365.25:.3f}/yr at S=S_bar")
print(f"  γ  = {gamma_h:.4f}")
print(f"  α  = {alpha_h:.4f}")
print(f"  β  = {beta_h:.4f}  → 1/β = {1/beta_h:.2f} d")
print(f"  κ  = {kappa_h:+.4f}  → G5 productivity = exp(κ) = {np.exp(kappa_h):.3f}×")
print(f"  η(G4) = {eta_g4:.3f},  η(G5) = {eta_g5:.3f}")
print(f"  log-L = {ll_v7:.2f}, AIC = {aic_v7:.2f}")

# Comparison Poisson on extended data
lam_pois = N/T
ll_pois = N*np.log(lam_pois) - lam_pois*T
print(f"  Reference Poisson: λ = {lam_pois*365.25:.2f}/yr,  log-L = {ll_pois:.2f},  ΔAIC vs v7 = {aic_v7 - (2 - 2*ll_pois):.1f}")

# ----------------------------------------------------------------------
# 7. Time-rescaling GOF
# ----------------------------------------------------------------------
mu_grid_h = mu0_h * (S_grid/S_bar)**gamma_h
cum_mu = np.concatenate([[0.0], np.cumsum(0.5*(mu_grid_h[:-1]+mu_grid_h[1:])*np.diff(grid))])
def int_mu(tq): return np.interp(tq, grid, cum_mu)
g_arr = np.exp(kappa_h*(m - m0))
Lam = np.zeros(N); A = 0.0; G = 0.0
for i, ti in enumerate(t):
    if i == 0:
        Lam[i] = int_mu(ti)
    else:
        dt = ti - t[i-1]
        A = np.exp(-beta_h*dt)*(A + g_arr[i-1])
        G = G + g_arr[i-1]
        Lam[i] = int_mu(ti) + (alpha_h/beta_h)*(G - A)
tau = np.diff(Lam)
ks = stats.kstest(tau, "expon", args=(0, 1.0))
lag1 = stats.pearsonr(tau[:-1], tau[1:])
print(f"\n[GOF] τ: mean={tau.mean():.3f}, var={tau.var():.3f}")
print(f"     KS p={ks.pvalue:.3e},  lag-1 r={lag1.statistic:+.3f} p={lag1.pvalue:.2e}")

# ----------------------------------------------------------------------
# 8. Block bootstrap CIs for (μ0, γ, α, β, κ)
#    Block = 365 days. Resample blocks with replacement covering total T.
# ----------------------------------------------------------------------
print("\n[bootstrap] block bootstrap CIs (block=365d, B=200)…")
BLOCK = 365.0
NB = int(np.ceil(T/BLOCK))
B = 200
rng_bs = np.random.default_rng(20260523)
params_bs = np.zeros((B, 5))
t_bs_start = time.time()
n_ok = 0
for b in range(B):
    # Sample NB blocks with replacement, each is a (start, end) window in original time
    starts_idx = rng_bs.integers(0, NB, size=NB)
    # Build resampled event list by shifting events within each chosen block
    new_t, new_m = [], []
    cur_T = 0.0
    for k, si in enumerate(starts_idx):
        b_start = si * BLOCK
        b_end = min(b_start + BLOCK, T)
        sel = (t >= b_start) & (t < b_end)
        if sel.any():
            ts = t[sel] - b_start + cur_T
            ms = m[sel]
            new_t.append(ts); new_m.append(ms)
        cur_T += (b_end - b_start)
    if not new_t:
        continue
    tb = np.concatenate(new_t)
    mb = np.concatenate(new_m)
    order_b = np.argsort(tb)
    tb = tb[order_b]; mb = mb[order_b]
    Tb = cur_T
    Sb = S(tb)
    grid_b = np.arange(0, Tb+1.0, 1.0)
    Sg_b = S(grid_b)
    try:
        res_b = optimize.minimize(neg_ll, best.x,
                                  args=(tb, mb, Tb, Sb, grid_b, Sg_b, S_bar, m0),
                                  method="Nelder-Mead",
                                  options={"xatol":1e-5,"fatol":1e-5,"maxiter":20000})
        if res_b.fun < 1e17 and np.all(np.isfinite(res_b.x)):
            params_bs[n_ok] = res_b.x
            n_ok += 1
    except Exception:
        pass
    if (b+1) % 50 == 0:
        print(f"  bootstrap {b+1}/{B}  ({n_ok} ok)  elapsed {time.time()-t_bs_start:.0f}s")

params_bs = params_bs[:n_ok]
print(f"\n[bootstrap] {n_ok}/{B} replicates succeeded ({time.time()-t_bs_start:.0f}s total)")
names_p = ["μ0", "γ", "α", "β", "κ"]
hat = [mu0_h, gamma_h, alpha_h, beta_h, kappa_h]
print(f"\n  {'Param':>6}  {'MLE':>10}  {'2.5%':>10}  {'97.5%':>10}  {'SE':>10}")
ci = {}
for i, name in enumerate(names_p):
    q = np.quantile(params_bs[:, i], [0.025, 0.5, 0.975])
    se = params_bs[:, i].std()
    ci[name] = (q[0], q[2])
    print(f"  {name:>6}  {hat[i]:10.5f}  {q[0]:10.5f}  {q[2]:10.5f}  {se:10.5f}")

# Derived: η(G4), η(G5), G5 productivity multiplier
eta_g4_bs = params_bs[:, 2] / params_bs[:, 3]
eta_g5_bs = params_bs[:, 2] * np.exp(params_bs[:, 4]) / params_bs[:, 3]
mult_bs = np.exp(params_bs[:, 4])
mu_yr_bs = params_bs[:, 0] * 365.25
inv_beta_bs = 1.0 / params_bs[:, 3]

print(f"\n  Derived quantities (95% CI):")
print(f"    μ0 (events/yr at mean SSN): {mu0_h*365.25:.3f}   "
      f"CI [{np.quantile(mu_yr_bs, 0.025):.3f}, {np.quantile(mu_yr_bs, 0.975):.3f}]")
print(f"    1/β (excitation half-life, days): {1/beta_h:.2f}   "
      f"CI [{np.quantile(inv_beta_bs, 0.025):.2f}, {np.quantile(inv_beta_bs, 0.975):.2f}]")
print(f"    η(G4): {eta_g4:.3f}   CI [{np.quantile(eta_g4_bs, 0.025):.3f}, {np.quantile(eta_g4_bs, 0.975):.3f}]")
print(f"    η(G5): {eta_g5:.3f}   CI [{np.quantile(eta_g5_bs, 0.025):.3f}, {np.quantile(eta_g5_bs, 0.975):.3f}]")
print(f"    G5 productivity exp(κ): {np.exp(kappa_h):.2f}   "
      f"CI [{np.quantile(mult_bs, 0.025):.2f}, {np.quantile(mult_bs, 0.975):.2f}]")

# Save bootstrap
np.save(os.path.join(DATA, "v7_bootstrap_params.npy"), params_bs)

# ----------------------------------------------------------------------
# 9. Per-cycle stability: refit on each 11-year cycle subset
# ----------------------------------------------------------------------
print("\n[stability] solar-cycle leave-one-out: refit dropping each cycle")
# SC11 starts 1867-12, SC25 starts 2019-12 (approx). Use canonical starts (years).
cycle_starts = [1867+11/12, 1878+11/12, 1890+3/12, 1901+11/12, 1913+7/12,
                1923+7/12, 1933+8/12, 1944+1/12, 1954+3/12, 1964+9/12,
                1976+2/12, 1986+8/12, 1996+7/12, 2008+11/12, 2019+11/12]
cycle_ids = list(range(11, 26))
# Convert to days from T0
cs_days = [(yfrac - T0_year) * 365.25 for yfrac in cycle_starts] + [T + 1]

print(f"  {'leave-out':>10}  {'μ0':>8}  {'γ':>6}  {'1/β':>6}  {'η(G4)':>7}  {'η(G5)':>7}  {'exp(κ)':>7}")
for k in range(len(cycle_ids)):
    a, b_ = cs_days[k], cs_days[k+1]
    if a < 0: a = 0
    if b_ > T: b_ = T
    keep = (t < a) | (t >= b_)
    if keep.sum() < 50:
        continue
    tk = t[keep]; mk = m[keep]
    Sk = S(tk)
    grid_k = np.arange(0, T+1.0, 1.0)
    Sg_k = S(grid_k)
    res_k = optimize.minimize(neg_ll, best.x,
                              args=(tk, mk, T, Sk, grid_k, Sg_k, S_bar, m0),
                              method="Nelder-Mead",
                              options={"xatol":1e-5,"fatol":1e-5,"maxiter":20000})
    if res_k.fun < 1e17:
        p = res_k.x
        print(f"     SC{cycle_ids[k]:>4}  {p[0]*365.25:8.2f}  {p[1]:6.3f}  {1/p[3]:6.2f}  "
              f"{p[2]/p[3]:7.3f}  {p[2]*np.exp(p[4])/p[3]:7.3f}  {np.exp(p[4]):7.2f}")

# ----------------------------------------------------------------------
# 10. Figures
# ----------------------------------------------------------------------
# (17) Extended event series
fig, ax = plt.subplots(figsize=(12, 4))
years_e = 1868 + t/365.25
mask_aa = source == "aa"
mask_kp = source == "kp"
mask_g5 = m >= 9.0
ax.vlines(years_e[mask_aa & ~mask_g5], 0, 1, color="#888", lw=0.5, alpha=0.7, label=f"G4 (aa-derived, n={(mask_aa & ~mask_g5).sum()})")
ax.vlines(years_e[mask_aa & mask_g5],  0, 2, color="#b1361e", lw=1.0, alpha=0.9, label=f"G5 (aa-derived, n={(mask_aa & mask_g5).sum()})")
ax.vlines(years_e[mask_kp & ~mask_g5], 0, 1, color="#5b8dba", lw=0.5, alpha=0.7, label=f"G4 (Kp-derived, n={(mask_kp & ~mask_g5).sum()})")
ax.vlines(years_e[mask_kp & mask_g5],  0, 2, color="#d43", lw=1.0, alpha=0.9, label=f"G5 (Kp-derived, n={(mask_kp & mask_g5).sum()})")
ax.axvline(1932, color="black", ls="--", lw=0.8, alpha=0.6)
ax.text(1932, 2.2, "1932: aa→Kp join", fontsize=9, ha="center")
ax.set_xlim(1868, 2026)
ax.set_ylim(0, 2.6)
ax.set_yticks([])
ax.set_xlabel("Year")
ax.set_title(f"Extended G4+ event series 1868–2025 (N={N}, {T/365.25:.0f} years, "
             f"~{(2026-1868)/11.1:.1f} solar cycles)")
ax.legend(loc="upper right", fontsize=9)
ax.grid(alpha=0.3, axis="x")
plt.tight_layout()
plt.savefig(os.path.join(FIG, "17_v7_extended_events.png"), dpi=140)
plt.close()
print("\n[fig] 17_v7_extended_events.png")

# (18) Bootstrap histograms
fig, axes = plt.subplots(2, 3, figsize=(13, 7))
labels = ["μ0 (events/yr at S̄)", "γ", "1/β (days)", "η(G4)", "η(G5)", "G5 productivity exp(κ)"]
vals = [mu_yr_bs, params_bs[:,1], inv_beta_bs, eta_g4_bs, eta_g5_bs, mult_bs]
hats = [mu0_h*365.25, gamma_h, 1/beta_h, eta_g4, eta_g5, np.exp(kappa_h)]
for ax, lab, v, h_ in zip(axes.flat, labels, vals, hats):
    ax.hist(v, bins=30, color="#5b8dba", alpha=0.65, edgecolor="white")
    ax.axvline(h_, color="#b1361e", lw=2, label=f"MLE = {h_:.3f}")
    q = np.quantile(v, [0.025, 0.975])
    ax.axvline(q[0], color="black", ls="--", lw=1)
    ax.axvline(q[1], color="black", ls="--", lw=1)
    ax.set_title(f"{lab}\n95% CI = [{q[0]:.3f}, {q[1]:.3f}]")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
plt.suptitle(f"v7 block-bootstrap parameter uncertainty (B={n_ok} replicates, block=365d)", fontsize=12, y=1.00)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "18_v7_bootstrap_cis.png"), dpi=140)
plt.close()
print("[fig] 18_v7_bootstrap_cis.png")

# ----------------------------------------------------------------------
# 11. Summary file
# ----------------------------------------------------------------------
with open(os.path.join(DATA, "hawkes_v7_summary.txt"), "w") as f:
    f.write("v7: Pre-1932 extension via aa-index — marked Hawkes 1868-2025\n")
    f.write("="*60 + "\n\n")
    f.write(f"Total events: {N} G4+ over {T/365.25:.1f} years\n")
    f.write(f"  Pre-1932 (aa-derived): {(source=='aa').sum()}\n")
    f.write(f"  1932+ (Kp-derived):    {(source=='kp').sum()}\n")
    f.write(f"  G5 days total: {(m>=9).sum()}  ({(mask_aa & mask_g5).sum()} pre-1932 + {(mask_kp & mask_g5).sum()} 1932+)\n\n")
    f.write(f"Rate-matched aa_max thresholds: G4+ ≥ {THR_G4:.0f}, G5 ≥ {THR_G5:.0f}\n")
    f.write(f"Day-level agreement on G4+ in overlap: {agree/n_g4*100:.1f}%\n\n")
    f.write(f"S_bar (1868-2025) = {S_bar:.2f}  (v6 used 93.50 for 1932-2025)\n\n")
    f.write("MLE parameters:\n")
    f.write(f"  μ0 = {mu0_h:.5f} /day ({mu0_h*365.25:.3f}/yr at S=S_bar)\n")
    f.write(f"  γ  = {gamma_h:.4f}\n")
    f.write(f"  α  = {alpha_h:.4f}\n")
    f.write(f"  β  = {beta_h:.4f} /day,  1/β = {1/beta_h:.2f} d\n")
    f.write(f"  κ  = {kappa_h:+.4f}, G5 productivity exp(κ) = {np.exp(kappa_h):.3f}×\n\n")
    f.write(f"  η(G4) = {eta_g4:.4f}\n")
    f.write(f"  η(G5) = {eta_g5:.4f}\n\n")
    f.write(f"Log-likelihood: {ll_v7:.2f},  AIC = {aic_v7:.2f}\n")
    f.write(f"GOF: KS p = {ks.pvalue:.3e},  lag-1 autocorr = {lag1.statistic:+.3f} (p={lag1.pvalue:.2e})\n\n")
    f.write(f"Block-bootstrap 95% CIs (B={n_ok}):\n")
    for i, name in enumerate(names_p):
        q = np.quantile(params_bs[:, i], [0.025, 0.975])
        f.write(f"  {name}: [{q[0]:.5f}, {q[1]:.5f}]\n")
    f.write(f"\n  Derived 95% CIs:\n")
    f.write(f"  μ0/yr:    [{np.quantile(mu_yr_bs,0.025):.3f}, {np.quantile(mu_yr_bs,0.975):.3f}]\n")
    f.write(f"  1/β (d):  [{np.quantile(inv_beta_bs,0.025):.2f}, {np.quantile(inv_beta_bs,0.975):.2f}]\n")
    f.write(f"  η(G4):    [{np.quantile(eta_g4_bs,0.025):.3f}, {np.quantile(eta_g4_bs,0.975):.3f}]\n")
    f.write(f"  η(G5):    [{np.quantile(eta_g5_bs,0.025):.3f}, {np.quantile(eta_g5_bs,0.975):.3f}]\n")
    f.write(f"  exp(κ):   [{np.quantile(mult_bs,0.025):.2f}, {np.quantile(mult_bs,0.975):.2f}]\n")

print("\n[done] v7 extended-record marked Hawkes complete")
