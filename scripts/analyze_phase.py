#!/usr/bin/env python3
"""
Solar-cycle-phase-conditioned hazard extension
==============================================
Builds on analyze.py.  Adds:
  - SILSO monthly sunspot record (1749-present)
  - Per-year solar-cycle phase tagging (min / rising / max / declining)
  - Phase-conditional Poisson rates + GPD parameters
  - Monte Carlo for decadal hazard under realistic Cycle 25/26 phase mixes
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.signal import find_peaks

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FIG  = os.path.join(ROOT, "figures")

# ----------------------------------------------------------------------
# 1. SILSO monthly sunspot number
# ----------------------------------------------------------------------
ssn_cols = ["year","month","frac_year","SN","SN_std","n_obs","provisional"]
ssn = pd.read_csv(
    os.path.join(DATA, "SN_m_tot_V2.0.txt"),
    sep=r"\s+",
    names=ssn_cols,
    engine="python",
    on_bad_lines="skip",
)
ssn = ssn[(ssn.year >= 1932) & (ssn.year <= 2025)].copy()
ssn["date"] = pd.to_datetime(dict(year=ssn.year, month=ssn.month, day=15))
ssn["SN"] = pd.to_numeric(ssn.SN, errors="coerce")
ssn = ssn.dropna(subset=["SN"])
print(f"[ssn] {len(ssn):,} monthly records {ssn.date.min().date()} → {ssn.date.max().date()}")

# 13-month smoothed sunspot number (standard SIDC definition)
ssn["SN_smooth"] = ssn.SN.rolling(window=13, center=True, min_periods=7).mean()

# ----------------------------------------------------------------------
# 2. Detect solar-cycle minima and maxima -> phase labels
# ----------------------------------------------------------------------
sn_smooth = ssn.SN_smooth.values
# Maxima: prominent peaks at least 7 years apart
max_idx, _ = find_peaks(sn_smooth, distance=7*12, prominence=20)
# Minima: invert and find peaks
min_idx, _ = find_peaks(-sn_smooth, distance=7*12, prominence=5)

cycle_max_dates = ssn.date.iloc[max_idx].tolist()
cycle_min_dates = ssn.date.iloc[min_idx].tolist()
print(f"[ssn] detected {len(cycle_max_dates)} maxima, {len(cycle_min_dates)} minima")
print("[ssn] maxima:", [d.strftime("%Y-%m") for d in cycle_max_dates])
print("[ssn] minima:", [d.strftime("%Y-%m") for d in cycle_min_dates])

def phase_for(date):
    """Tag a date with its position in the 11-year cycle."""
    # find surrounding min and max
    prev_min = max([d for d in cycle_min_dates if d <= date], default=None)
    next_min = min([d for d in cycle_min_dates if d >  date], default=None)
    nearest_max = min(cycle_max_dates, key=lambda d: abs((d-date).days))
    if prev_min is None or next_min is None:
        return "unknown"
    # rising = prev_min to max; declining = max to next_min
    # min phase = within ±18 months of a minimum
    # max phase = within ±18 months of a maximum
    if abs((nearest_max - date).days) <= 18*30:
        return "max"
    if min(abs((prev_min - date).days), abs((next_min - date).days)) <= 18*30:
        return "min"
    if date < nearest_max:
        return "rising"
    return "declining"

# ----------------------------------------------------------------------
# 3. Re-load Kp/ap data and merge phase
# ----------------------------------------------------------------------
cols = ["year","month","day","hh","hh_m","days","days_m","Kp","ap","D"]
kp = pd.read_csv(
    os.path.join(DATA, "Kp_ap_since_1932.txt"),
    sep=r"\s+", comment="#", names=cols, engine="python",
)
kp["date"] = pd.to_datetime(dict(year=kp.year, month=kp.month, day=kp.day))
kp = kp[(kp.year >= 1932) & (kp.year <= 2025) & (kp.Kp >= 0)].copy()
daily = kp.groupby("date").agg(
    Kp_max=("Kp","max"),
    ap_max=("ap","max"),
).reset_index()
daily["year"] = daily.date.dt.year

print(f"[kp] {len(daily):,} daily aggregates")
print("[merge] tagging phases (this takes ~30s)…")
daily["phase"] = daily.date.apply(phase_for)
print(daily.phase.value_counts())

daily.to_csv(os.path.join(DATA, "derived_daily_with_phase.csv"), index=False)

# ----------------------------------------------------------------------
# 4. Phase-conditional Poisson rate + GPD
# ----------------------------------------------------------------------
threshold = np.quantile(daily.ap_max, 0.95)
print(f"\n[fit] common 95th-pct threshold ap = {threshold:.1f}")

phase_stats = {}
for phase in ["min","rising","max","declining"]:
    sub = daily[daily.phase == phase]
    n_years = sub.date.dt.year.nunique()
    fraction_of_record = len(sub) / len(daily)
    exc = sub.ap_max[sub.ap_max > threshold].values - threshold
    lam = len(exc) / (fraction_of_record * 94)   # rate per year-of-this-phase
    if len(exc) >= 30:
        shape, _, scale = stats.genpareto.fit(exc, floc=0)
    else:
        shape, scale = -0.02, 65.0
    phase_stats[phase] = dict(
        n_days=len(sub),
        n_exc=len(exc),
        lam_per_year=lam,
        gpd_shape=shape,
        gpd_scale=scale,
        fraction_of_record=fraction_of_record,
    )
    print(f"  {phase:10s}  days={len(sub):>6}  exc={len(exc):>4}  "
          f"λ={lam:5.2f}/yr-of-phase  ξ={shape:+.3f}  σ={scale:5.2f}")

pd.DataFrame(phase_stats).T.to_csv(os.path.join(DATA, "derived_phase_stats.csv"))

# ----------------------------------------------------------------------
# 5. Monte Carlo — decadal worst storm under a representative phase mix
# ----------------------------------------------------------------------
# A typical 10-year window spans roughly: 2 yr min + 3 yr rising + 2 yr max + 3 yr declining
PHASE_MIX = {"min":2.0, "rising":3.0, "max":2.0, "declining":3.0}

rng = np.random.default_rng(20260523)
N_TRIALS = 30_000
carrington_excess = 400.0 - threshold

p_carr_phase  = np.zeros(N_TRIALS, dtype=bool)   # one realistic decade
p_carr_uncond = np.zeros(N_TRIALS, dtype=bool)   # using global rate (re-derived)

# also build phase-only decades for contrast
p_carr_all_max     = np.zeros(N_TRIALS, dtype=bool)   # 10 years entirely at max
p_carr_all_min     = np.zeros(N_TRIALS, dtype=bool)   # 10 years entirely at min

# global rate
all_exc = (daily.ap_max[daily.ap_max > threshold] - threshold).values
g_shape, _, g_scale = stats.genpareto.fit(all_exc, floc=0)
g_lam = len(all_exc) / 94.0

def draw_max(n, shape, scale):
    if n == 0: return 0.0
    return stats.genpareto.rvs(shape, 0, scale, size=n, random_state=rng).max()

for i in range(N_TRIALS):
    # realistic mix
    worst = 0.0
    for ph, yrs in PHASE_MIX.items():
        ps = phase_stats[ph]
        n  = rng.poisson(ps["lam_per_year"] * yrs)
        m  = draw_max(n, ps["gpd_shape"], ps["gpd_scale"])
        if m > worst: worst = m
    p_carr_phase[i] = worst >= carrington_excess
    # global rate
    n = rng.poisson(g_lam * 10)
    p_carr_uncond[i] = draw_max(n, g_shape, g_scale) >= carrington_excess
    # all-max decade
    ps = phase_stats["max"]
    n = rng.poisson(ps["lam_per_year"] * 10)
    p_carr_all_max[i] = draw_max(n, ps["gpd_shape"], ps["gpd_scale"]) >= carrington_excess
    # all-min decade
    ps = phase_stats["min"]
    n = rng.poisson(ps["lam_per_year"] * 10)
    p_carr_all_min[i] = draw_max(n, ps["gpd_shape"], ps["gpd_scale"]) >= carrington_excess

print("\n[MC] decadal P(>=1 Carrington-class) by scenario:")
print(f"  realistic mix (2/3/2/3 yrs):  {p_carr_phase.mean():.3f}")
print(f"  unconditional (global rate):  {p_carr_uncond.mean():.3f}")
print(f"  decade entirely at solar max: {p_carr_all_max.mean():.3f}")
print(f"  decade entirely at solar min: {p_carr_all_min.mean():.3f}")

# ----------------------------------------------------------------------
# 6. Figures
# ----------------------------------------------------------------------
# (a) phase-conditional storm-day rate
phase_order = ["min","rising","max","declining"]
g3p_rate = []
g4p_rate = []
for ph in phase_order:
    sub = daily[daily.phase == ph]
    yrs = ps["fraction_of_record"] = len(sub) / 365.25
    g3p_rate.append((sub.Kp_max >= 7).sum() / yrs)
    g4p_rate.append((sub.Kp_max >= 8).sum() / yrs)

fig, ax = plt.subplots(figsize=(8, 4.5))
x = np.arange(len(phase_order))
ax.bar(x-0.2, g3p_rate, width=0.4, color="#f0a04b", label="G3+ days / year")
ax.bar(x+0.2, g4p_rate, width=0.4, color="#b1361e", label="G4+ days / year")
ax.set_xticks(x); ax.set_xticklabels(phase_order)
ax.set_ylabel("Storm days per year in this phase")
ax.set_title("Storm-day rate by solar-cycle phase  (1932–2025)")
ax.legend(); ax.grid(alpha=0.3, axis="y")
plt.tight_layout()
plt.savefig(os.path.join(FIG, "04_phase_storm_rate.png"), dpi=140)
plt.close()
print("[fig] 04_phase_storm_rate.png")

# (b) decadal hazard bar chart
fig, ax = plt.subplots(figsize=(8, 4.5))
scenarios = ["Solar min\n(10yr)", "Realistic mix\n(2/3/2/3)",
             "Unconditional\n(global rate)", "Solar max\n(10yr)"]
probs = [p_carr_all_min.mean(), p_carr_phase.mean(),
         p_carr_uncond.mean(),  p_carr_all_max.mean()]
colors = ["#5b8dba","#4f8c4f","#888888","#b1361e"]
bars = ax.bar(scenarios, probs, color=colors)
for b, p in zip(bars, probs):
    ax.text(b.get_x()+b.get_width()/2, p+0.01, f"{p:.1%}",
            ha="center", fontsize=11, fontweight="bold")
ax.set_ylim(0, 1)
ax.set_ylabel("P(>=1 ap >= 400 in 10 years)")
ax.set_title("Decadal Carrington-class hazard by solar-cycle phase mix")
ax.grid(alpha=0.3, axis="y")
plt.tight_layout()
plt.savefig(os.path.join(FIG, "05_phase_hazard_compare.png"), dpi=140)
plt.close()
print("[fig] 05_phase_hazard_compare.png")

# (c) SSN smoothed + storm overlay
fig, ax = plt.subplots(figsize=(11, 4.5))
ax.plot(ssn.date, ssn.SN_smooth, color="#5b8dba", lw=1.4, label="13-mo smoothed SSN")
ax.fill_between(ssn.date, 0, ssn.SN_smooth, color="#5b8dba", alpha=0.15)
g4_days = daily[daily.Kp_max >= 8]
ax.scatter(g4_days.date, [-20]*len(g4_days), marker="|",
           color="#b1361e", s=80, label=f"G4+ days (n={len(g4_days)})")
ax.set_ylabel("Smoothed sunspot number")
ax.set_xlabel("Year")
ax.set_title("Solar cycle (SSN) and G4+ geomagnetic storm days, 1932–2025")
ax.legend(loc="upper right")
ax.grid(alpha=0.3)
ax.set_ylim(-40, 280)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "06_ssn_with_storms.png"), dpi=140)
plt.close()
print("[fig] 06_ssn_with_storms.png")

# Summary file append
with open(os.path.join(DATA, "run_summary.txt"), "a") as f:
    f.write("\n\n=== PHASE-CONDITIONED EXTENSION ===\n")
    f.write(f"Realistic-mix decadal P(Carrington): {p_carr_phase.mean():.4f}\n")
    f.write(f"All-max decadal P(Carrington):       {p_carr_all_max.mean():.4f}\n")
    f.write(f"All-min decadal P(Carrington):       {p_carr_all_min.mean():.4f}\n")
    f.write(f"Unconditional (re-derived):          {p_carr_uncond.mean():.4f}\n")
    for ph, s in phase_stats.items():
        f.write(f"  {ph}: λ={s['lam_per_year']:.2f}/yr  "
                f"ξ={s['gpd_shape']:+.3f}  σ={s['gpd_scale']:.2f}\n")

print("\n[done] phase-conditioned extension complete")
