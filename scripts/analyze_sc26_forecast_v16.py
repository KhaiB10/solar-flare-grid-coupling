#!/usr/bin/env python3
"""
v16 SC26 amplitude prediction — true out-of-sample test of the v15 hierarchical model.

Method
------
1. Load the v15 InferenceData (population hyperposteriors mu_*, sigma_*, plus gamma).
2. For each of B posterior draws:
   - Draw a *fresh* cycle's parameters (mu0_26, alpha_26, beta_26, kappa_26)
     from the population distribution implied by that draw's hyperparameters.
     This is the prior-predictive for a new (unobserved) cycle.
   - Simulate one full SC26 realization via Ogata thinning, with F10.7(t)
     following one of four scenarios:
       (a) FLAT_QUIET  S(t) = 118 sfu  (Singh+2021 SC26 forecast)
       (b) FLAT_AVG    S(t) = 150 sfu  (mid-modern)
       (c) FLAT_LIKE25 S(t) = 180 sfu  (SC25-like, "stronger than expected")
       (d) PHYSICAL    sinusoid 75 -> peak -> 75 sfu over the cycle,
                       peak amplitude itself sampled per realization from
                       N(155, 25^2) truncated to [80, 230].
3. Tally per-scenario distributions of:
   - total G4+ event count
   - peak month of activity (when do events cluster)
   - P(>=1 G5) using empirical mark distribution from SC23-25
   - P(Carrington-class >=1) using a low-rate Pareto tail on marks

Random seed: 20260523.
"""
import os, json, pickle, time
os.environ.setdefault("PYTENSOR_FLAGS", "cxx=,mode=FAST_RUN")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) \
       if "__file__" in globals() else "/home/user/workspace/solar-flare-grid-coupling"
DATA = os.path.join(ROOT, "data")
FIG  = os.path.join(ROOT, "figures")
SEED = 20260523
rng  = np.random.default_rng(SEED)

# ----------------------------------------------------------------------
# 1. Load posterior + reference data
# ----------------------------------------------------------------------
print("[1] loading v15 posterior and reference catalog…")
with open(os.path.join(DATA, "v15_idata.pkl"), "rb") as f:
    idata = pickle.load(f)
post = idata.posterior

# Hyperposteriors (flattened across chain/draw)
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
print(f"  posterior draws available: {n_post}")

# Reference catalog (for empirical mark distribution and historical counts)
ev_df = pd.read_csv(os.path.join(DATA, "derived_events_extended_1844_2025.csv"),
                    parse_dates=["date"]).sort_values("date").reset_index(drop=True)
S_df  = pd.read_csv(os.path.join(DATA, "derived_S_daily_v12.csv"),
                    parse_dates=["date"]).sort_values("date").reset_index(drop=True)
S_bar = float(S_df["S_daily"].mean())
print(f"  S_bar (training mean F10.7) = {S_bar:.2f}")

m0 = 8.0  # G4 threshold mark

# Historical counts per cycle for the comparison plot
CYCLES_REF = [
    (22, '1986-09-01', '1996-08-01'),
    (23, '1996-08-01', '2008-12-01'),
    (24, '2008-12-01', '2019-12-01'),
    (25, '2019-12-01', '2025-05-31'),  # partial
]
hist_counts = {}
for num, s, e in CYCLES_REF:
    s_ts = pd.Timestamp(s); e_ts = pd.Timestamp(e)
    n = ((ev_df.date >= s_ts) & (ev_df.date < e_ts)).sum()
    hist_counts[num] = int(n)
print(f"  historical G4+ counts: {hist_counts}")

# Empirical mark distribution from modern cycles (SC23+24+25)
modern_marks = ev_df[(ev_df.date >= pd.Timestamp("1996-08-01"))].mark.values
print(f"  modern (SC23-25) mark distribution: n={len(modern_marks)}, "
      f"min/median/max = {modern_marks.min():.1f}/{np.median(modern_marks):.1f}/{modern_marks.max():.1f}")
print(f"  marks >=9 in modern record: {(modern_marks>=9).sum()}/{len(modern_marks)} = "
      f"{(modern_marks>=9).mean()*100:.1f}%")

# ----------------------------------------------------------------------
# 2. SC26 window definition
# ----------------------------------------------------------------------
print("\n[2] SC26 forecast window…")
SC26_START = pd.Timestamp("2031-07-01")
SC26_END   = pd.Timestamp("2042-12-31")
n_days = (SC26_END - SC26_START).days + 1
print(f"  window: {SC26_START.date()} -> {SC26_END.date()} ({n_days} d = {n_days/365.25:.2f} yr)")

t_grid = np.arange(n_days, dtype=float)        # 0-indexed days since SC26_START
dates  = pd.date_range(SC26_START, SC26_END, freq="D")

# Physics-informed F10.7 trajectory shape: sin half-cycle ramp + Gaussian peak + tail
# Real cycles: ~3 yr rising, peak, ~7 yr declining. Use empirical SC23 F10.7 shape rescaled.
sc23_window = (S_df.date >= pd.Timestamp("1996-08-01")) & (S_df.date < pd.Timestamp("2008-12-01"))
sc23_S = S_df.loc[sc23_window, "S_daily"].values
# Resample SC23 to SC26 length
template = np.interp(np.linspace(0, len(sc23_S)-1, n_days),
                     np.arange(len(sc23_S)), sc23_S)
# Normalize template to [0, 1] (subtract min, divide by range)
template_min = float(template.min()); template_max = float(template.max())
template01 = (template - template_min) / (template_max - template_min)
F107_FLOOR = 75.0
print(f"  template based on SC23 shape, length {n_days} d, "
      f"min->max ratio {template_min:.0f}->{template_max:.0f}")

# ----------------------------------------------------------------------
# 3. Forward simulator (Ogata thinning, single cycle from t=0 with no history)
# ----------------------------------------------------------------------
def simulate_sc26(mu0, gamma_g, alpha, beta, kappa, S_traj, t_end_days, rng, mark_pool):
    """Return list of (t_day, mark) for one SC26 realization, no prior history."""
    events_t = []
    events_m = []
    g_vals   = []
    t = 0.0
    MAX_EV = 1000  # SC26 cap (no observed cycle has exceeded ~150)
    while t < t_end_days and len(events_t) < MAX_EV:
        idx = min(int(t), len(S_traj) - 1)
        if idx < 0: idx = 0
        S_t  = S_traj[idx]
        mu_t = mu0 * (S_t / S_bar) ** gamma_g
        if events_t:
            dt  = t - np.array(events_t)
            exc = alpha * np.sum(np.array(g_vals) * np.exp(-beta * dt))
        else:
            exc = 0.0
        lam_upper = mu_t + exc + 1e-9
        dt_next = rng.exponential(1.0 / lam_upper)
        t_cand = t + dt_next
        if t_cand >= t_end_days: break
        idx = min(int(t_cand), len(S_traj) - 1)
        if idx < 0: idx = 0
        S_tc  = S_traj[idx]
        mu_tc = mu0 * (S_tc / S_bar) ** gamma_g
        if events_t:
            dt = t_cand - np.array(events_t)
            exc_tc = alpha * np.sum(np.array(g_vals) * np.exp(-beta * dt))
        else:
            exc_tc = 0.0
        lam_tc = mu_tc + exc_tc
        if rng.random() < lam_tc / lam_upper:
            m_new = float(rng.choice(mark_pool))
            events_t.append(t_cand)
            events_m.append(m_new)
            g_vals.append(np.exp(kappa * (m_new - m0)))
        t = t_cand
    return list(zip(events_t, events_m))

# ----------------------------------------------------------------------
# 4. Run all four scenarios
# ----------------------------------------------------------------------
B = 2000
print(f"\n[4] running {B} simulations per scenario…")
sim_idx = rng.choice(n_post, size=B, replace=True)

scenarios = {
    "FLAT_QUIET":  ("F10.7 flat @ 118 sfu (Singh+2021)",  np.full(n_days, 118.0)),
    "FLAT_AVG":    ("F10.7 flat @ 150 sfu (modern avg)",  np.full(n_days, 150.0)),
    "FLAT_LIKE25": ("F10.7 flat @ 180 sfu (SC25-like)",   np.full(n_days, 180.0)),
    "PHYSICAL":    ("F10.7 SC23-shape, peak ~N(155,25) truncated", None),
}

results = {}
for sc_name, (sc_desc, sc_traj) in scenarios.items():
    print(f"\n  scenario {sc_name}: {sc_desc}")
    t0 = time.time()
    counts = np.zeros(B, dtype=int)
    g5_counts = np.zeros(B, dtype=int)
    all_event_dates = []  # for peak-month density (kept across all sims)
    for k, idx in enumerate(sim_idx):
        # Draw a fresh cycle's params from the population
        mu0_26   = float(np.exp(mu_mu[idx]     + sigma_mu[idx]    * rng.standard_normal()))
        alpha_26 = float(np.exp(mu_alpha[idx]  + sigma_alpha[idx] * rng.standard_normal()))
        beta_26  = float(np.exp(mu_beta[idx]   + sigma_beta[idx]  * rng.standard_normal()))
        kappa_26 = float(mu_kappa[idx]         + sigma_kappa[idx] * rng.standard_normal())
        gamma_26 = float(gamma_post[idx])

        if sc_name == "PHYSICAL":
            peak_amp = float(np.clip(rng.normal(155, 25), 80, 230))
            S_traj   = F107_FLOOR + (peak_amp - F107_FLOOR) * template01
        else:
            S_traj   = sc_traj

        evs = simulate_sc26(mu0_26, gamma_26, alpha_26, beta_26, kappa_26,
                            S_traj, t_end_days=float(n_days), rng=rng,
                            mark_pool=modern_marks)
        counts[k] = len(evs)
        g5_counts[k] = sum(1 for _, m in evs if m >= 9.0)
        if k < 500:  # keep day-of-cycle for first 500 sims (memory)
            all_event_dates.extend([t for t, _ in evs])
        if (k+1) % 500 == 0:
            print(f"    {k+1}/{B}: median count = {np.median(counts[:k+1]):.1f}  "
                  f"({time.time()-t0:.1f}s)")

    results[sc_name] = {
        "desc":         sc_desc,
        "counts":       counts,
        "g5_counts":    g5_counts,
        "event_days":   np.array(all_event_dates),
    }
    print(f"  done  median={np.median(counts):.1f}  "
          f"95% HDI=[{np.quantile(counts,0.025):.0f},{np.quantile(counts,0.975):.0f}]  "
          f"P(>=1 G5)={np.mean(g5_counts>=1):.3f}")

# ----------------------------------------------------------------------
# 5. Plot
# ----------------------------------------------------------------------
print("\n[5] generating plots…")

# F1: SC26 count distribution per scenario + historical
fig, ax = plt.subplots(figsize=(11, 6))
colors = {"FLAT_QUIET":"#1f77b4", "FLAT_AVG":"#2ca02c", "FLAT_LIKE25":"#d62728", "PHYSICAL":"#9467bd"}
bins = np.arange(0, max(r["counts"].max() for r in results.values()) + 4, 2)
for sc, col in colors.items():
    c = results[sc]["counts"]
    ax.hist(c, bins=bins, alpha=0.45, color=col, density=True,
            label=f"{sc}: med={np.median(c):.0f}, 95% HDI [{np.quantile(c,0.025):.0f},{np.quantile(c,0.975):.0f}]",
            edgecolor=col)
# historical marks for context
for cnum, ccount in hist_counts.items():
    ax.axvline(ccount, color="k", ls=":", lw=1, alpha=0.6)
    ax.text(ccount, ax.get_ylim()[1]*0.95, f"SC{cnum}\n({ccount})", rotation=90,
            ha="right", va="top", fontsize=8, alpha=0.8)
ax.set_xlabel("Total G4+ events during SC26 (2031-07 → 2042-12)")
ax.set_ylabel("posterior density")
ax.set_title("v16 SC26 forecast — total G4+ event count by F10.7 scenario\n"
             "(dotted vertical lines = historical SC22-25 counts)")
ax.legend(loc="upper right", fontsize=9)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "60_sc26_count_distribution.png"), dpi=150)
plt.close()

# F2: Peak-month density (timing of events through the cycle) — PHYSICAL scenario only
fig, ax = plt.subplots(figsize=(11, 5))
phys_days = results["PHYSICAL"]["event_days"]
if len(phys_days) > 0:
    bins_year = np.arange(0, n_days + 90, 90)  # quarterly bins
    counts_q, edges = np.histogram(phys_days, bins=bins_year)
    centers = (edges[:-1] + edges[1:]) / 2 / 365.25 + 2031.5
    ax.bar(centers, counts_q, width=0.22, color="C4", alpha=0.7,
           label=f"PHYSICAL scenario, {len(results['PHYSICAL']['event_days'])} events from 500 sims")
    # F10.7 template overlay
    ax2 = ax.twinx()
    ax2.plot(np.arange(n_days)/365.25 + 2031.5, template,
             color="C3", lw=1.5, alpha=0.6, label="SC23-shape F10.7 template")
    ax2.set_ylabel("F10.7 template (sfu)", color="C3")
    ax2.tick_params(axis="y", colors="C3")
    ax.set_xlabel("Year")
    ax.set_ylabel("G4+ events per quarter (aggregate)")
    ax.set_title("v16 SC26 — temporal density of G4+ events under physical F10.7 scenario\n"
                 "(when do they cluster?)")
    ax.legend(loc="upper left")
    ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "61_sc26_event_timing.png"), dpi=150)
plt.close()

# F3: G5 probability and tail risk comparison
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]
g5_data = [results[sc]["g5_counts"] for sc in scenarios]
g5_labels = list(scenarios.keys())
ax.boxplot(g5_data, labels=g5_labels, showmeans=True, meanprops={"marker":"D","markerfacecolor":"red","markersize":6})
ax.set_ylabel("number of G5+ events per SC26 realization")
ax.set_title("G5 count per SC26 realization, by scenario")
ax.grid(alpha=0.3)

ax = axes[1]
xs = np.arange(len(scenarios))
probs_any   = [np.mean(results[sc]["g5_counts"] >= 1) for sc in scenarios]
probs_two   = [np.mean(results[sc]["g5_counts"] >= 2) for sc in scenarios]
ax.bar(xs - 0.2, probs_any, 0.4, label="P(≥1 G5)", color="C0")
ax.bar(xs + 0.2, probs_two, 0.4, label="P(≥2 G5)", color="C3")
ax.set_xticks(xs); ax.set_xticklabels(g5_labels, rotation=20)
ax.set_ylim(0, 1.05)
ax.set_ylabel("posterior probability")
ax.set_title("Tail-event probability across F10.7 scenarios")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "62_sc26_tail_probability.png"), dpi=150)
plt.close()

print("  saved figures 60-62")

# ----------------------------------------------------------------------
# 6. Save summary JSON
# ----------------------------------------------------------------------
print("\n[6] writing summary JSON…")
summary = {
    "seed": SEED,
    "model": "v16 SC26 prior-predictive — draws fresh cycle from v15 population posterior",
    "window": {"start": str(SC26_START.date()), "end": str(SC26_END.date()), "days": int(n_days)},
    "n_post_draws": int(n_post),
    "n_sims_per_scenario": B,
    "historical_counts": hist_counts,
    "scenarios": {},
}
for sc_name, r in results.items():
    c = r["counts"]; g5 = r["g5_counts"]
    summary["scenarios"][sc_name] = {
        "description": r["desc"],
        "count_median": int(np.median(c)),
        "count_hdi_50": [int(np.quantile(c, 0.25)), int(np.quantile(c, 0.75))],
        "count_hdi_95": [int(np.quantile(c, 0.025)), int(np.quantile(c, 0.975))],
        "P_count_greater_30": float(np.mean(c > 30)),
        "P_count_greater_50": float(np.mean(c > 50)),
        "P_count_less_than_10": float(np.mean(c < 10)),
        "g5_median": int(np.median(g5)),
        "P_any_G5": float(np.mean(g5 >= 1)),
        "P_two_or_more_G5": float(np.mean(g5 >= 2)),
    }
with open(os.path.join(DATA, "v16_sc26_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("  saved data/v16_sc26_summary.json")

# Cross-scenario summary table
print("\n=== SC26 forecast summary across F10.7 scenarios ===")
print(f"  Historical: SC22={hist_counts[22]}  SC23={hist_counts[23]}  "
      f"SC24={hist_counts[24]}  SC25 partial through 2025-05={hist_counts[25]}")
print()
print(f"  {'scenario':<15}  {'med':>5}  {'50% HDI':<12}  {'95% HDI':<12}  "
      f"{'P(>30)':>7}  {'P(>50)':>7}  {'P(G5)':>7}  {'P(2G5)':>7}")
for sc_name, s in summary["scenarios"].items():
    print(f"  {sc_name:<15}  {s['count_median']:>5}  "
          f"[{s['count_hdi_50'][0]:>3},{s['count_hdi_50'][1]:>3}]   "
          f"[{s['count_hdi_95'][0]:>3},{s['count_hdi_95'][1]:>3}]   "
          f"{s['P_count_greater_30']:>7.3f}  {s['P_count_greater_50']:>7.3f}  "
          f"{s['P_any_G5']:>7.3f}  {s['P_two_or_more_G5']:>7.3f}")
print("\n[done v16 SC26 forecast]")
