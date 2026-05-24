#!/usr/bin/env python3
"""
Geomagnetic storm clustering at two timescales — pressure-test.

Key questions:
  Q1. Is G4+ wait-time distribution non-exponential? (KS test, AIC vs Poisson)
  Q2. Is there a separable two-component structure: short (CME-pulse) + 27-day (Bartels)?
  Q3. What is the "effective" decadal hazard if we model clustering correctly?
  Q4. How does the cluster-aware decadal disturbance-window count compare to Poisson?
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats, optimize

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FIG  = os.path.join(ROOT, "figures")

daily = pd.read_csv(os.path.join(DATA, "derived_daily_with_phase.csv"),
                    parse_dates=["date"])
daily = daily.sort_values("date").reset_index(drop=True)

# ---- G4+ event series ----
g4 = daily[daily.Kp_max >= 8].copy()
print(f"G4+ days 1932-2025: {len(g4)}")
waits = g4.date.diff().dt.days.dropna().values
print(f"Waiting times: {len(waits)}")
print(f"  Min: {waits.min()} d (same-day excluded since we work daily)")
print(f"  Mean: {waits.mean():.1f} d  Median: {np.median(waits):.0f} d")

# ----------------------------------------------------------------------
# Q1. Fit & compare distributions to waiting times
# ----------------------------------------------------------------------
print("\n=== Q1: Wait-time distribution comparison ===")
# Candidate distributions
results = {}

# Exponential (Poisson process baseline)
loc, scale = stats.expon.fit(waits, floc=0)
ll_exp = np.sum(stats.expon.logpdf(waits, loc, scale))
aic_exp = 2*1 - 2*ll_exp
results["exponential"] = (ll_exp, aic_exp, {"scale": scale})

# Gamma (allows clustering via shape < 1)
shape, loc, scale = stats.gamma.fit(waits, floc=0)
ll_gamma = np.sum(stats.gamma.logpdf(waits, shape, loc, scale))
aic_gamma = 2*2 - 2*ll_gamma
results["gamma"] = (ll_gamma, aic_gamma, {"shape": shape, "scale": scale})

# Weibull
c, loc, scale = stats.weibull_min.fit(waits, floc=0)
ll_wb = np.sum(stats.weibull_min.logpdf(waits, c, loc, scale))
aic_wb = 2*2 - 2*ll_wb
results["weibull"] = (ll_wb, aic_wb, {"c": c, "scale": scale})

# Lognormal (heavy tail)
s, loc, scale = stats.lognorm.fit(waits, floc=0)
ll_ln = np.sum(stats.lognorm.logpdf(waits, s, loc, scale))
aic_ln = 2*2 - 2*ll_ln
results["lognormal"] = (ll_ln, aic_ln, {"s": s, "scale": scale})

# Power-law tail (Pareto)
b, loc, scale = stats.pareto.fit(waits, floc=0)
ll_pl = np.sum(stats.pareto.logpdf(waits, b, loc, scale))
aic_pl = 2*2 - 2*ll_pl
results["pareto"] = (ll_pl, aic_pl, {"b": b, "scale": scale})

# Mixture: 2-component exponential (fast cluster + slow background)
def neg_ll_mix2(params, x):
    p, lam1, lam2 = params
    if not (0 < p < 1 and lam1 > 0 and lam2 > 0):
        return 1e10
    pdf = p * lam1*np.exp(-lam1*x) + (1-p) * lam2*np.exp(-lam2*x)
    pdf = np.maximum(pdf, 1e-300)
    return -np.sum(np.log(pdf))

x0 = [0.5, 1/10, 1/200]
res = optimize.minimize(neg_ll_mix2, x0, args=(waits,), method="Nelder-Mead",
                        options={"xatol": 1e-6, "fatol":1e-6, "maxiter":5000})
p, lam1, lam2 = res.x
ll_mix2 = -res.fun
aic_mix2 = 2*3 - 2*ll_mix2
mean1, mean2 = 1/lam1, 1/lam2
results["mixture_exp2"] = (ll_mix2, aic_mix2,
                           {"weight_fast": p, "mean_fast_d": mean1, "mean_slow_d": mean2})

# Print AIC ranking
print(f"{'Distribution':<18} {'LL':>10} {'AIC':>10} {'ΔAIC':>10}  params")
ranked = sorted(results.items(), key=lambda kv: kv[1][1])
best_aic = ranked[0][1][1]
for name, (ll, aic, par) in ranked:
    d = aic - best_aic
    p_str = ", ".join(f"{k}={v:.3f}" for k,v in par.items())
    print(f"  {name:<18} {ll:>10.1f} {aic:>10.1f} {d:>10.1f}  {p_str}")

# KS vs exponential (formal test)
rate = 1/scale  # not really used; recompute properly
exp_scale = stats.expon.fit(waits, floc=0)[1]
ks = stats.kstest(waits, "expon", args=(0, exp_scale))
print(f"\nKS test vs exponential: D={ks.statistic:.3f}, p={ks.pvalue:.2e}")

# ----------------------------------------------------------------------
# Q2. Two-timescale structure
# ----------------------------------------------------------------------
print("\n=== Q2: Two-timescale (CME-cluster + Bartels) structure ===")
print(f"Mixture-exponential best fit:")
print(f"  Fast component:  weight={p:.3f},  mean wait = {mean1:.1f} days")
print(f"  Slow component:  weight={1-p:.3f},  mean wait = {mean2:.1f} days")
print(f"  AIC improvement over single exponential: {aic_exp - aic_mix2:.1f}")
print(f"  (rule of thumb: ΔAIC > 10 = decisive)")

# Periodogram of inter-arrival check
within_5  = (waits <= 5).sum()
within_10 = (waits <= 10).sum()
within_27 = (waits <= 27).sum()
within_35 = ((waits >= 20) & (waits <= 35)).sum()
total = len(waits)
print(f"\nFraction of waits:")
print(f"  ≤  5 days: {within_5/total:.3f}  ({within_5})")
print(f"  ≤ 10 days: {within_10/total:.3f}  ({within_10})")
print(f"  20–35 days (Bartels band): {within_35/total:.3f}  ({within_35})")
print(f"  ≤ 27 days (one Bartels):   {within_27/total:.3f}  ({within_27})")

# ----------------------------------------------------------------------
# Q3 & Q4. Cluster-aware decadal hazard: count "storm events" rather than
# "storm days" — collapse any G4+ days within 7 days into one event
# ----------------------------------------------------------------------
print("\n=== Q3: Cluster-aware event count ===")
# Define an "event" as a maximal run of G4+ days where consecutive G4+
# days are within EVENT_GAP days of each other.
EVENT_GAP = 7  # storms within a week are one CME-driven event
g4_dates = g4.date.values
events = []
cur = [g4_dates[0]]
for d in g4_dates[1:]:
    if (d - cur[-1])/np.timedelta64(1,"D") <= EVENT_GAP:
        cur.append(d)
    else:
        events.append(cur); cur = [d]
events.append(cur)
print(f"  Raw G4+ days:    {len(g4_dates)}")
print(f"  G4+ clusters (≤{EVENT_GAP}d apart):  {len(events)}")
print(f"  Cluster size dist: mean={np.mean([len(e) for e in events]):.2f}, "
      f"max={max(len(e) for e in events)}")

cluster_size_dist = pd.Series([len(e) for e in events]).value_counts().sort_index()
print(f"  Cluster size counts:\n{cluster_size_dist.to_string()}")

# Decadal hazard reframed
n_years = 94
events_per_year = len(events) / n_years
days_per_year   = len(g4_dates) / n_years
print(f"\n  G4+ DAYS per year:   {days_per_year:.2f}")
print(f"  G4+ EVENTS per year: {events_per_year:.2f}  (a {days_per_year/events_per_year:.2f}× collapse)")

# ----------------------------------------------------------------------
# Q4. What does cluster-aware Monte Carlo say about decadal RECOVERY-WINDOW count?
# ----------------------------------------------------------------------
# Naive Poisson day-based: lambda_days * 10 = expected G4+ days per decade
# Cluster-aware:         lambda_events * 10 = expected independent CME-driven events
print("\n=== Q4: Decadal storm-event Monte Carlo ===")
rng = np.random.default_rng(20260523)
N = 50_000
poisson_decade_days  = rng.poisson(days_per_year*10,  size=N)
poisson_decade_evts  = rng.poisson(events_per_year*10, size=N)

print(f"  Poisson (DAYS counted) decadal G4+:")
print(f"    mean={poisson_decade_days.mean():.1f},  ≥10 in: {(poisson_decade_days>=10).mean():.3f}")
print(f"  Cluster-aware (EVENTS counted) decadal G4+:")
print(f"    mean={poisson_decade_evts.mean():.1f},  ≥10 in: {(poisson_decade_evts>=10).mean():.3f}")
print(f"  Effective reduction in 'independent shocks' per decade: "
      f"{poisson_decade_days.mean()/poisson_decade_evts.mean():.2f}×")

# ----------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

# (a) Wait-time histogram + fitted exponential + mixture
ax = axes[0]
hist_bins = np.linspace(0, 200, 41)
ax.hist(waits[waits<=200], bins=hist_bins, density=True,
        color="#5b8dba", edgecolor="white", alpha=0.85,
        label=f"Observed (n={(waits<=200).sum()})")
xx = np.linspace(0.1, 200, 400)
ax.plot(xx, stats.expon.pdf(xx, 0, exp_scale), color="#888", lw=2,
        label=f"Exponential  λ⁻¹={exp_scale:.0f} d")
mix_pdf = p*lam1*np.exp(-lam1*xx) + (1-p)*lam2*np.exp(-lam2*xx)
ax.plot(xx, mix_pdf, color="#b1361e", lw=2.5,
        label=f"2-exp mix: fast={mean1:.1f}d ({p:.0%}), slow={mean2:.0f}d")
ax.set_xlabel("Days between successive G4+ storm days")
ax.set_ylabel("Density")
ax.set_title("G4+ wait-time distribution, 1932–2025")
ax.legend()
ax.grid(alpha=0.3)

# (b) Empirical survival vs exponential (log scale)
ax = axes[1]
sorted_w = np.sort(waits)
sf = 1 - np.arange(1, len(sorted_w)+1)/(len(sorted_w)+1)
ax.semilogy(sorted_w, sf, "o", ms=3, color="#5b8dba", alpha=0.7,
            label="Observed")
ax.semilogy(xx, stats.expon.sf(xx, 0, exp_scale), color="#888", lw=2,
            label=f"Exponential (Poisson process)")
ax.semilogy(xx, p*np.exp(-lam1*xx) + (1-p)*np.exp(-lam2*xx),
            color="#b1361e", lw=2.5, label="2-exp mixture")
ax.set_xlabel("Days t")
ax.set_ylabel("P(wait > t)")
ax.set_title("Survival function — log scale")
ax.legend()
ax.grid(alpha=0.3, which="both")
ax.set_xlim(0, 500)
ax.set_ylim(1e-3, 1.1)

plt.tight_layout()
plt.savefig(os.path.join(FIG, "08_clustering_waittime.png"), dpi=140)
plt.close()
print("\n[fig] 08_clustering_waittime.png")

# (c) Cluster-collapsed timeline
fig, ax = plt.subplots(figsize=(12, 3.5))
for e in events:
    ax.scatter(e[0], len(e), color="#b1361e" if len(e)>=3 else "#f0a04b",
               s=30+10*len(e), edgecolor="black", linewidth=0.4, alpha=0.85)
ax.set_xlabel("Year")
ax.set_ylabel("G4+ days in cluster")
ax.set_title(f"G4+ storm clusters, 1932–2025  "
             f"(n_clusters={len(events)}, raw days={len(g4_dates)})")
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "09_clusters_timeline.png"), dpi=140)
plt.close()
print("[fig] 09_clusters_timeline.png")

# Save derived
pd.DataFrame({
    "cluster_start": [e[0] for e in events],
    "cluster_end":   [e[-1] for e in events],
    "n_days":        [len(e) for e in events],
}).to_csv(os.path.join(DATA, "derived_g4_clusters.csv"), index=False)

# Save summary
with open(os.path.join(DATA, "clustering_summary.txt"), "w") as f:
    f.write("Geomagnetic G4+ wait-time clustering analysis\n")
    f.write("=============================================\n\n")
    f.write(f"G4+ days 1932-2025: {len(g4)}\n")
    f.write(f"Wait-times: n={len(waits)}, mean={waits.mean():.1f}d, median={np.median(waits):.0f}d\n\n")
    f.write("AIC ranking (lower = better):\n")
    for name, (ll, aic, par) in ranked:
        f.write(f"  {name:<18} AIC={aic:.1f}  ΔAIC={aic-best_aic:.1f}\n")
    f.write(f"\nKS vs exponential: D={ks.statistic:.3f}, p={ks.pvalue:.2e}\n\n")
    f.write(f"2-exponential mixture:\n")
    f.write(f"  Fast component: weight={p:.3f}, mean={mean1:.1f}d (CME-cluster timescale)\n")
    f.write(f"  Slow component: weight={1-p:.3f}, mean={mean2:.0f}d (background+Bartels)\n\n")
    f.write(f"Cluster collapse (gap={EVENT_GAP}d):\n")
    f.write(f"  Raw G4+ days: {len(g4_dates)}\n")
    f.write(f"  G4+ clusters: {len(events)}\n")
    f.write(f"  Collapse factor: {len(g4_dates)/len(events):.2f}x\n")

print("[done]")
