#!/usr/bin/env python3
"""Re-plot SC26 count distribution with sensible x-axis clipping and clearer labels.
Rerun the simulation lightly to recover counts (or reload from a fresh pickle).
"""
import os, json, pickle, time
os.environ.setdefault("PYTENSOR_FLAGS", "cxx=,mode=FAST_RUN")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = "/home/user/workspace/solar-flare-grid-coupling"
DATA = os.path.join(ROOT, "data")
FIG  = os.path.join(ROOT, "figures")
SEED = 20260523
rng  = np.random.default_rng(SEED)

with open(os.path.join(DATA, "v15_idata.pkl"), "rb") as f:
    idata = pickle.load(f)
post = idata.posterior
mu_mu     = post["mu_mu"].values.ravel()
sigma_mu  = post["sigma_mu"].values.ravel()
mu_alpha  = post["mu_alpha"].values.ravel()
sigma_alpha = post["sigma_alpha"].values.ravel()
mu_beta   = post["mu_beta"].values.ravel()
sigma_beta = post["sigma_beta"].values.ravel()
mu_kappa  = post["mu_kappa"].values.ravel()
sigma_kappa = post["sigma_kappa"].values.ravel()
gamma_post = post["gamma"].values.ravel()
n_post = mu_mu.size

ev_df = pd.read_csv(os.path.join(DATA, "derived_events_extended_1844_2025.csv"),
                    parse_dates=["date"]).sort_values("date").reset_index(drop=True)
S_df  = pd.read_csv(os.path.join(DATA, "derived_S_daily_v12.csv"),
                    parse_dates=["date"]).sort_values("date").reset_index(drop=True)
S_bar = float(S_df["S_daily"].mean())
modern_marks = ev_df[ev_df.date >= pd.Timestamp("1996-08-01")].mark.values
m0 = 8.0
n_days = (pd.Timestamp("2042-12-31") - pd.Timestamp("2031-07-01")).days + 1

sc23_window = (S_df.date >= pd.Timestamp("1996-08-01")) & (S_df.date < pd.Timestamp("2008-12-01"))
sc23_S = S_df.loc[sc23_window, "S_daily"].values
template = np.interp(np.linspace(0, len(sc23_S)-1, n_days), np.arange(len(sc23_S)), sc23_S)
template01 = (template - template.min()) / (template.max() - template.min())
F107_FLOOR = 75.0

def simulate_sc26(mu0, gamma_g, alpha, beta, kappa, S_traj, t_end_days, rng, mark_pool):
    events_t = []; events_m = []; g_vals = []
    t = 0.0; MAX_EV = 1000
    while t < t_end_days and len(events_t) < MAX_EV:
        idx = min(int(t), len(S_traj) - 1); idx = max(idx, 0)
        S_t = S_traj[idx]; mu_t = mu0 * (S_t / S_bar) ** gamma_g
        if events_t:
            dt = t - np.array(events_t)
            exc = alpha * np.sum(np.array(g_vals) * np.exp(-beta * dt))
        else:
            exc = 0.0
        lam_upper = mu_t + exc + 1e-9
        dt_next = rng.exponential(1.0 / lam_upper); t_cand = t + dt_next
        if t_cand >= t_end_days: break
        idx = min(int(t_cand), len(S_traj) - 1); idx = max(idx, 0)
        S_tc = S_traj[idx]; mu_tc = mu0 * (S_tc / S_bar) ** gamma_g
        if events_t:
            dt = t_cand - np.array(events_t)
            exc_tc = alpha * np.sum(np.array(g_vals) * np.exp(-beta * dt))
        else:
            exc_tc = 0.0
        lam_tc = mu_tc + exc_tc
        if rng.random() < lam_tc / lam_upper:
            m_new = float(rng.choice(mark_pool))
            events_t.append(t_cand); events_m.append(m_new)
            g_vals.append(np.exp(kappa * (m_new - m0)))
        t = t_cand
    return list(zip(events_t, events_m))

B = 2000
sim_idx = rng.choice(n_post, size=B, replace=True)
scenarios = {
    "PHYSICAL":    ("F10.7 SC23-shape, peak ~N(155, 25²) sfu",  None,  "#9467bd"),
    "FLAT_QUIET":  ("F10.7 flat @ 118 sfu (Singh+2021 forecast)",  np.full(n_days, 118.0),  "#1f77b4"),
    "FLAT_AVG":    ("F10.7 flat @ 150 sfu (mid-modern)",            np.full(n_days, 150.0),  "#2ca02c"),
    "FLAT_LIKE25": ("F10.7 flat @ 180 sfu (SC25-like)",             np.full(n_days, 180.0),  "#d62728"),
}

print("re-running simulations (fast)…")
counts_all = {}
t0 = time.time()
for sc_name, (desc, S_flat, col) in scenarios.items():
    counts = np.zeros(B, dtype=int)
    for k, idx in enumerate(sim_idx):
        mu0_26   = float(np.exp(mu_mu[idx]     + sigma_mu[idx]    * rng.standard_normal()))
        alpha_26 = float(np.exp(mu_alpha[idx]  + sigma_alpha[idx] * rng.standard_normal()))
        beta_26  = float(np.exp(mu_beta[idx]   + sigma_beta[idx]  * rng.standard_normal()))
        kappa_26 = float(mu_kappa[idx]         + sigma_kappa[idx] * rng.standard_normal())
        gamma_26 = float(gamma_post[idx])
        if sc_name == "PHYSICAL":
            peak_amp = float(np.clip(rng.normal(155, 25), 80, 230))
            S_traj   = F107_FLOOR + (peak_amp - F107_FLOOR) * template01
        else:
            S_traj = S_flat
        evs = simulate_sc26(mu0_26, gamma_26, alpha_26, beta_26, kappa_26,
                            S_traj, float(n_days), rng, modern_marks)
        counts[k] = len(evs)
    counts_all[sc_name] = counts
    print(f"  {sc_name}: median={np.median(counts):.0f}  ({time.time()-t0:.1f}s)")

hist_counts = {22:28, 23:31, 24:3, 25:10}

# --- Plot A: bounded x-axis, side-by-side scenarios ---
fig, ax = plt.subplots(figsize=(12, 6))
X_MAX = 120
bins = np.arange(0, X_MAX + 4, 3)
order = ["PHYSICAL", "FLAT_QUIET", "FLAT_AVG", "FLAT_LIKE25"]
for sc in order:
    desc, _, col = scenarios[sc]
    c = counts_all[sc].copy()
    c_clip = np.minimum(c, X_MAX)
    med = np.median(c); q025, q975 = np.quantile(c, [0.025, 0.975])
    ax.hist(c_clip, bins=bins, alpha=0.45, color=col, density=True, edgecolor=col,
            label=f"{sc} — median={med:.0f}, 95% HDI [{q025:.0f}, {q975:.0f}]")
# Historical reference (with proper label placement)
ymax = ax.get_ylim()[1]
hist_y_positions = {22: 0.92, 23: 0.85, 24: 0.78, 25: 0.71}
for cnum, ccount in hist_counts.items():
    ax.axvline(ccount, color="k", ls=":", lw=1.2, alpha=0.75)
    label_suffix = " (partial)" if cnum == 25 else ""
    ax.text(ccount + 0.5, ymax * hist_y_positions[cnum],
            f"SC{cnum}={ccount}{label_suffix}",
            ha="left", va="top", fontsize=9, alpha=0.85,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="gray", alpha=0.7))
ax.set_xlim(0, X_MAX)
ax.set_xlabel("Total G4+ events during SC26 (2031-07 → 2042-12)")
ax.set_ylabel("posterior density")
ax.set_title("v16 — SC26 G4+ count forecast across F10.7 scenarios\n"
             "(historical SC22–25 counts shown as dotted vertical lines; X-axis clipped at 120 for clarity)")
ax.legend(loc="upper right", fontsize=10, framealpha=0.95)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "60_sc26_count_distribution.png"), dpi=150)
plt.close()
print("saved fig 60 (replot)")

# --- Plot B: ECDF view, cleaner for comparing scenarios ---
fig, ax = plt.subplots(figsize=(12, 5.5))
for sc in order:
    desc, _, col = scenarios[sc]
    c = np.sort(counts_all[sc])
    ecdf = np.arange(1, len(c)+1) / len(c)
    ax.step(c, ecdf, where="post", color=col, lw=2,
            label=f"{sc}: med={np.median(c):.0f}")
for cnum, ccount in hist_counts.items():
    ax.axvline(ccount, color="k", ls=":", lw=1, alpha=0.6)
    ax.text(ccount, 0.02, f"SC{cnum}", rotation=90, ha="right", va="bottom",
            fontsize=8, alpha=0.7)
ax.set_xlim(0, 100); ax.set_ylim(0, 1)
ax.set_xlabel("G4+ event count")
ax.set_ylabel("cumulative probability  P(SC26 count ≤ X)")
ax.set_title("v16 — SC26 G4+ count cumulative distribution by F10.7 scenario")
ax.legend(loc="lower right", fontsize=10)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "63_sc26_count_ecdf.png"), dpi=150)
plt.close()
print("saved fig 63")
print("done")
