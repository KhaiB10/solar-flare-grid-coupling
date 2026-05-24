#!/usr/bin/env python3
"""
Solar Flare → Grid Coupling Analysis
=====================================
Defensive publication analysis for Diatom Sky R&D.

Inputs
------
data/Kp_ap_since_1932.txt   : GFZ Potsdam Kp/ap 3-hour index, 1932-present
data/known_gic_grid_events.csv : Curated table of documented GIC grid events

Outputs
-------
figures/*.png : Plots
data/derived_*.csv : Derived tables
FINDINGS.md (written by separate step)

Methodology
-----------
1. Parse Kp/ap 3-hour record (1932-2025, ~94 years).
2. Aggregate to daily and storm-day stats.
3. Fit a generalized Pareto tail to ap daily maxima above the 95th percentile.
4. Monte Carlo: sample N=10,000 future 10-year windows under the fitted
   tail to estimate the probability of at least one "Carrington-equivalent"
   storm (ap >= 400, equivalent to Kp ~ 9o) per decade.
5. Compare modeled storm-day frequencies to documented grid impact events
   to estimate the conditional probability P(grid impact | storm intensity).

Caveats (see FINDINGS.md): This is an open-data exploratory analysis, NOT
an operational risk model. We use the ap index (linear) rather than Kp
(quasi-log) for tail fitting because Pareto tails require a linear scale.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FIG  = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

# ----------------------------------------------------------------------
# 1. Load Kp/ap index
# ----------------------------------------------------------------------
# Format (from GFZ header): YYYY MM DD hh.h hh._m days days_m Kp ap D
cols = ["year","month","day","hh","hh_m","days","days_m","Kp","ap","D"]
df = pd.read_csv(
    os.path.join(DATA, "Kp_ap_since_1932.txt"),
    sep=r"\s+",
    comment="#",
    names=cols,
    engine="python",
)
df["date"] = pd.to_datetime(dict(year=df.year, month=df.month, day=df.day))
df = df[(df.year >= 1932) & (df.year <= 2025)].copy()

# Drop placeholder rows where Kp = -1 (definitive flag not set)
df = df[df.Kp >= 0]

print(f"[load] {len(df):,} three-hour Kp records from "
      f"{df.date.min().date()} to {df.date.max().date()}")

# ----------------------------------------------------------------------
# 2. Daily aggregates
# ----------------------------------------------------------------------
daily = df.groupby("date").agg(
    Kp_max=("Kp", "max"),
    Kp_sum=("Kp", "sum"),
    ap_max=("ap", "max"),
    ap_mean=("ap", "mean"),
).reset_index()
daily["year"] = daily.date.dt.year

# G-scale classification (NOAA): G1 Kp=5, G2=6, G3=7, G4=8, G5=9
def g_scale(kp):
    if kp >= 9: return 5
    if kp >= 8: return 4
    if kp >= 7: return 3
    if kp >= 6: return 2
    if kp >= 5: return 1
    return 0

daily["G"] = daily.Kp_max.apply(g_scale)
daily.to_csv(os.path.join(DATA, "derived_daily.csv"), index=False)

storms_per_year = daily.groupby(["year","G"]).size().unstack(fill_value=0)
storms_per_year.to_csv(os.path.join(DATA, "derived_storms_per_year.csv"))
print("\n[storms per year, last 10 yrs]")
print(storms_per_year.tail(10))

# ----------------------------------------------------------------------
# 3. Solar-cycle envelope
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(11, 4.5))
yearly = daily.groupby("year").agg(
    severe=("G", lambda x: (x >= 3).sum()),
    extreme=("G", lambda x: (x >= 4).sum()),
).reset_index()
ax.bar(yearly.year, yearly.severe, color="#f0a04b", label="G3+ days (Strong)")
ax.bar(yearly.year, yearly.extreme, color="#b1361e", label="G4+ days (Severe/Extreme)")
ax.set_xlabel("Year")
ax.set_ylabel("Storm days per year")
ax.set_title("Geomagnetic Storm Days per Year, 1932–2025  (Kp ≥ 7 = G3; Kp ≥ 8 = G4+)")
ax.legend(loc="upper right")
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "01_storm_days_per_year.png"), dpi=140)
plt.close()
print("[fig] 01_storm_days_per_year.png")

# ----------------------------------------------------------------------
# 4. ap-tail Generalized Pareto fit
# ----------------------------------------------------------------------
# ap_max daily values, threshold = 95th percentile
ap_vals = daily.ap_max.dropna().values
threshold = np.quantile(ap_vals, 0.95)
excess = ap_vals[ap_vals > threshold] - threshold
shape, loc, scale = stats.genpareto.fit(excess, floc=0)
print(f"\n[GPD] threshold ap={threshold:.1f}  shape={shape:.3f}  scale={scale:.2f}")
print(f"[GPD] exceedances: {len(excess):,} of {len(ap_vals):,} days "
      f"({100*len(excess)/len(ap_vals):.2f}%)")

# Empirical vs fitted tail
fig, ax = plt.subplots(figsize=(8, 5))
sorted_excess = np.sort(excess)
empirical_sf = 1.0 - np.arange(1, len(sorted_excess)+1) / (len(sorted_excess)+1)
ax.semilogy(sorted_excess + threshold, empirical_sf, "o", ms=3,
            alpha=0.5, label="Empirical (1932–2025)")
xx = np.linspace(0, sorted_excess.max()*1.1, 200)
ax.semilogy(xx + threshold, stats.genpareto.sf(xx, shape, 0, scale),
            "-", lw=2, color="#b1361e", label=f"GPD fit (ξ={shape:.2f})")
ax.axvline(400, ls="--", color="gray", label="ap=400 (Carrington-class)")
ax.set_xlabel("Daily max ap index")
ax.set_ylabel("Survival probability (log)")
ax.set_title("Tail of daily ap distribution — POT fit above 95th percentile")
ax.legend()
ax.grid(alpha=0.3, which="both")
plt.tight_layout()
plt.savefig(os.path.join(FIG, "02_ap_tail_fit.png"), dpi=140)
plt.close()
print("[fig] 02_ap_tail_fit.png")

# ----------------------------------------------------------------------
# 5. Monte Carlo — decadal Carrington probability
# ----------------------------------------------------------------------
# Rate of exceedances above threshold per year
n_years = daily.year.nunique()
lam = len(excess) / n_years   # Poisson rate
print(f"\n[MC] exceedance rate λ = {lam:.2f} / year")

rng = np.random.default_rng(20260523)
N_TRIALS = 20_000
DECADE = 10

# For each decade trial: draw Poisson(λ * 10) exceedance counts, then
# GPD-distributed magnitudes; record whether any sample exceeds (400 - threshold).
carrington_excess = 400.0 - threshold
extreme_excess    = 600.0 - threshold   # ~ a 1.5× Carrington event

p_carr = np.zeros(N_TRIALS, dtype=bool)
p_ext  = np.zeros(N_TRIALS, dtype=bool)
worst  = np.zeros(N_TRIALS)

for i in range(N_TRIALS):
    n = rng.poisson(lam * DECADE)
    if n == 0:
        worst[i] = threshold
        continue
    samples = stats.genpareto.rvs(shape, 0, scale, size=n, random_state=rng)
    worst[i] = samples.max() + threshold
    if samples.max() >= carrington_excess:
        p_carr[i] = True
    if samples.max() >= extreme_excess:
        p_ext[i]  = True

print(f"[MC] P(≥1 ap≥400 in a decade) = {p_carr.mean():.3f}  "
      f"({p_carr.sum()}/{N_TRIALS})")
print(f"[MC] P(≥1 ap≥600 in a decade) = {p_ext.mean():.4f}")
print(f"[MC] median decadal worst ap = {np.median(worst):.0f}")
print(f"[MC] 95th-pct decadal worst ap = {np.quantile(worst, 0.95):.0f}")
print(f"[MC] 99th-pct decadal worst ap = {np.quantile(worst, 0.99):.0f}")

# Plot
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(worst, bins=60, color="#5b8dba", edgecolor="white")
ax.axvline(400, ls="--", color="#b1361e", label="ap=400 (Carrington-class)")
ax.axvline(np.median(worst), ls=":", color="black",
           label=f"Median worst = {np.median(worst):.0f}")
ax.set_xlabel("Worst daily ap in a 10-year window")
ax.set_ylabel("Monte Carlo trials")
ax.set_title(f"Monte Carlo: distribution of worst decadal ap  "
             f"(N={N_TRIALS:,} trials, λ={lam:.1f}/yr)")
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "03_monte_carlo_decadal.png"), dpi=140)
plt.close()
print("[fig] 03_monte_carlo_decadal.png")

# ----------------------------------------------------------------------
# 6. Documented grid events overlay
# ----------------------------------------------------------------------
events = pd.read_csv(os.path.join(DATA, "known_gic_grid_events.csv"))
events["date"] = pd.to_datetime(events["date"])
events = events.merge(daily[["date","ap_max","Kp_max"]], on="date", how="left")
events.to_csv(os.path.join(DATA, "derived_events_with_ap.csv"), index=False)
print("\n[events] documented GIC grid events with ap/Kp:")
print(events[["date","event_name","Kp_max","ap_max","min_dst_nT"]].to_string(index=False))

# Save run summary
with open(os.path.join(DATA, "run_summary.txt"), "w") as f:
    f.write(f"Records: {len(df):,}\n")
    f.write(f"Daily aggregates: {len(daily):,}\n")
    f.write(f"Years covered: {n_years}\n")
    f.write(f"GPD threshold ap: {threshold:.1f}\n")
    f.write(f"GPD shape (xi): {shape:.4f}\n")
    f.write(f"GPD scale: {scale:.4f}\n")
    f.write(f"Exceedance rate (per year): {lam:.4f}\n")
    f.write(f"P(>=1 Carrington-class in decade): {p_carr.mean():.4f}\n")
    f.write(f"P(>=1 1.5x-Carrington in decade): {p_ext.mean():.4f}\n")
    f.write(f"Median decadal worst ap: {np.median(worst):.0f}\n")
    f.write(f"95th-pct decadal worst ap: {np.quantile(worst, 0.95):.0f}\n")
    f.write(f"99th-pct decadal worst ap: {np.quantile(worst, 0.99):.0f}\n")

print("\n[done] all outputs written")
