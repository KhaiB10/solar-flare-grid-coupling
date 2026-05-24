#!/usr/bin/env python3
"""
Exploratory: hunt for novel patterns in 94 years of Kp/ap + SSN.

Hypotheses to test:
  H1. Semiannual / equinoctial effect — strength + statistical significance
  H2. Lunar / synodic month modulation (recently claimed in literature, contested)
  H3. Phase asymmetry within cycle — declining > rising at matched SSN level
  H4. "Storm clustering" — are G4+ days self-exciting (Hawkes-like)?
  H5. Cycle-to-cycle change in tail thickness — is the tail getting heavier?
  H6. Weekly / 27-day Bartels rotation signal in G4+ days
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FIG  = os.path.join(ROOT, "figures")

# Reload daily aggregates with phase
daily = pd.read_csv(os.path.join(DATA, "derived_daily_with_phase.csv"),
                    parse_dates=["date"])
daily["doy"] = daily.date.dt.dayofyear
daily["month"] = daily.date.dt.month
daily["year"] = daily.date.dt.year
daily["dow"] = daily.date.dt.dayofweek

print(f"Loaded {len(daily):,} daily records")
print(f"G4+ days: {(daily.Kp_max >= 8).sum()}")
print(f"G3+ days: {(daily.Kp_max >= 7).sum()}")

# ----------------------------------------------------------------------
# H1. Semiannual / equinoctial
# ----------------------------------------------------------------------
print("\n=== H1: Semiannual / equinoctial ===")
monthly_rate = daily.groupby("month").agg(
    n_days=("Kp_max", "size"),
    g3p=("Kp_max", lambda x: (x>=7).sum()),
    g4p=("Kp_max", lambda x: (x>=8).sum()),
    mean_ap=("ap_max","mean"),
).reset_index()
monthly_rate["g3p_per_year"] = monthly_rate.g3p / (monthly_rate.n_days/30.4) * 12
monthly_rate["g4p_per_year"] = monthly_rate.g4p / (monthly_rate.n_days/30.4) * 12
print(monthly_rate[["month","g3p","g4p","mean_ap"]])

# Test: are March/September higher than December/June?
equinox = daily[daily.month.isin([3,4,9,10])]
solstice = daily[daily.month.isin([6,7,12,1])]
eq_rate = (equinox.Kp_max>=7).mean()
sol_rate = (solstice.Kp_max>=7).mean()
# binomial test
from scipy.stats import binomtest
total_g3 = (daily.Kp_max>=7).sum()
n_eq_g3 = (equinox.Kp_max>=7).sum()
expected_p = len(equinox) / len(daily)
bt = binomtest(n_eq_g3, total_g3, expected_p, alternative="greater")
print(f"  Equinox G3+ rate: {eq_rate:.4f}, Solstice: {sol_rate:.4f}, ratio: {eq_rate/sol_rate:.2f}")
print(f"  Binomial test (equinox enrichment): p = {bt.pvalue:.2e}")

# ----------------------------------------------------------------------
# H2. Lunar / synodic modulation
# ----------------------------------------------------------------------
print("\n=== H2: Lunar synodic month (29.53 d) modulation ===")
# Days since new moon 1932-01-01 reference. Use a known new moon: 2000-01-06 18:14 UTC
ref = pd.Timestamp("2000-01-06")
synodic = 29.530588
daily["lunar_phase"] = ((daily.date - ref).dt.days % synodic) / synodic  # 0=new, 0.5=full

# Bin into 8 phases
bins = np.linspace(0, 1, 9)
daily["lunar_bin"] = np.digitize(daily.lunar_phase, bins) - 1
lunar = daily.groupby("lunar_bin").agg(
    n=("Kp_max","size"),
    g3p=("Kp_max", lambda x: (x>=7).sum()),
    g4p=("Kp_max", lambda x: (x>=8).sum()),
    mean_ap=("ap_max","mean"),
).reset_index()
print(lunar)
# Chi-square uniformity
expected = lunar.n * (daily.Kp_max>=7).sum() / len(daily)
chi2, p = stats.chisquare(lunar.g3p, expected)
print(f"  Chi-square uniformity of G3+ across lunar phases: χ²={chi2:.2f}, p={p:.3f}")

# ----------------------------------------------------------------------
# H3. Rising vs declining phase asymmetry at matched SSN
# ----------------------------------------------------------------------
print("\n=== H3: Rising vs declining at matched SSN ===")
ssn_m = pd.read_csv(os.path.join(DATA,"SN_m_tot_V2.0.txt"),
                    sep=r"\s+", header=None,
                    names=["year","month","frac","SN","SN_std","n","prov"],
                    engine="python", on_bad_lines="skip")
ssn_m = ssn_m[(ssn_m.year>=1932) & (ssn_m.year<=2025)]
ssn_m["date"] = pd.to_datetime(dict(year=ssn_m.year, month=ssn_m.month, day=15))
ssn_m["SN_smooth"] = ssn_m.SN.rolling(13, center=True, min_periods=7).mean()
# Monthly merge
daily["ym"] = daily.date.dt.to_period("M").dt.to_timestamp() + pd.Timedelta(days=14)
ssn_m["ym"] = ssn_m.date
merged = daily.merge(ssn_m[["ym","SN_smooth"]], on="ym", how="left")

# Bin by SSN
ssn_bins = [0, 25, 50, 75, 100, 150, 300]
merged["ssn_bin"] = pd.cut(merged.SN_smooth, ssn_bins)
phase_ssn = merged.groupby(["ssn_bin","phase"], observed=True).agg(
    n=("Kp_max","size"),
    g3p_rate=("Kp_max", lambda x: (x>=7).mean()*365.25),
).reset_index()
print(phase_ssn[phase_ssn.phase.isin(["rising","declining"])].to_string())

# ----------------------------------------------------------------------
# H4. Storm clustering — waiting time between G4+ days
# ----------------------------------------------------------------------
print("\n=== H4: Waiting times between G4+ days ===")
g4 = daily[daily.Kp_max >= 8].sort_values("date").reset_index(drop=True)
waits = g4.date.diff().dt.days.dropna().values
print(f"  N waiting times: {len(waits)}")
print(f"  Mean: {waits.mean():.1f} d, Median: {np.median(waits):.0f} d")
print(f"  Fraction within 7 days: {(waits<=7).mean():.3f}")
print(f"  Fraction within 27 days (one Bartels): {(waits<=27).mean():.3f}")
# Expected under Poisson with same rate
rate_per_day = len(g4) / len(daily)
exp_within_7  = 1 - np.exp(-rate_per_day * 7)
exp_within_27 = 1 - np.exp(-rate_per_day * 27)
print(f"  Poisson expected within 7 d: {exp_within_7:.3f}")
print(f"  Poisson expected within 27 d: {exp_within_27:.3f}")
# KS test against exponential
ks = stats.kstest(waits, "expon", args=(0, 1/rate_per_day))
print(f"  KS vs exponential: D={ks.statistic:.3f}, p={ks.pvalue:.2e}")
# Variance ratio (Fano factor) — Poisson should be 1
print(f"  Variance/mean of waits = {waits.var()/waits.mean():.2f}  (1.0 = Poisson)")

# Also do this for the 27-day Bartels rotation: do G4+ days recur near 27, 54, 81 d gaps?
print("  Wait-time histogram (5-day bins, first 100 d):")
hist, edges = np.histogram(waits[waits<=100], bins=np.arange(0,101,5))
for h, e in zip(hist, edges[:-1]):
    bar = "#" * int(h)
    print(f"    {e:3d}–{e+5:3d} d: {h:3d} {bar}")

# ----------------------------------------------------------------------
# H5. Cycle-to-cycle tail thickness drift
# ----------------------------------------------------------------------
print("\n=== H5: Per-cycle GPD shape (tail thickness) ===")
# Cycles 17 (~1933) → 25 (~2019+)
cycle_boundaries = [
    (1933, 1944, "C17"),
    (1944, 1954, "C18"),
    (1954, 1964, "C19"),
    (1964, 1976, "C20"),
    (1976, 1986, "C21"),
    (1986, 1996, "C22"),
    (1996, 2008, "C23"),
    (2008, 2019, "C24"),
    (2019, 2026, "C25"),
]
thresh = np.quantile(daily.ap_max, 0.95)
cycle_fits = []
for y0, y1, name in cycle_boundaries:
    sub = daily[(daily.year>=y0) & (daily.year<y1)]
    exc = (sub.ap_max[sub.ap_max>thresh] - thresh).values
    if len(exc) < 30:
        cycle_fits.append((name, y0, y1, len(exc), None, None, (sub.Kp_max>=8).sum()))
        continue
    shape, _, scale = stats.genpareto.fit(exc, floc=0)
    cycle_fits.append((name, y0, y1, len(exc), shape, scale, (sub.Kp_max>=8).sum()))
    print(f"  {name} ({y0}-{y1}): n_exc={len(exc):>4}  ξ={shape:+.3f}  σ={scale:5.2f}  G4+={(sub.Kp_max>=8).sum()}")

# ----------------------------------------------------------------------
# H6. Day-of-week / Bartels 27-d periodicity in G4+ days
# ----------------------------------------------------------------------
print("\n=== H6: Bartels 27-day periodicity ===")
# Power spectrum of daily G4+ indicator
indicator = (daily.Kp_max>=7).astype(float).values
N = len(indicator)
mean = indicator.mean()
fft = np.fft.rfft(indicator - mean)
freqs = np.fft.rfftfreq(N, d=1.0)   # cycles per day
power = np.abs(fft)**2
periods = 1/freqs[1:]
# Look in 20-35 day range
mask = (periods>=20) & (periods<=35)
peak_period = periods[mask][np.argmax(power[1:][mask])]
print(f"  Strongest period in 20-35 d band: {peak_period:.1f} d")

# Plot the spectrum near 27d
fig, ax = plt.subplots(figsize=(9,4))
mask2 = (periods>=10) & (periods<=200)
ax.semilogx(periods[mask2], power[1:][mask2], color="#4f8c4f")
ax.axvline(27.0, ls="--", color="#b1361e", label="27-day Bartels rotation")
ax.axvline(182.5, ls=":", color="#888", label="182.5-day semiannual")
ax.set_xlabel("Period (days)")
ax.set_ylabel("Spectral power (G3+ indicator)")
ax.set_title("Power spectrum of daily G3+ storm-day indicator, 1932-2025")
ax.legend()
ax.grid(alpha=0.3, which="both")
plt.tight_layout()
plt.savefig(os.path.join(FIG, "07_power_spectrum.png"), dpi=140)
plt.close()
print("  Figure: 07_power_spectrum.png")

# Save key derived data
pd.DataFrame(cycle_fits, columns=["cycle","y0","y1","n_exc","gpd_shape","gpd_scale","g4plus"]
            ).to_csv(os.path.join(DATA, "derived_per_cycle_fits.csv"), index=False)

print("\n=== EXPLORATION COMPLETE ===")
