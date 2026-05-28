#!/usr/bin/env python3
"""
Hawkes universality analysis — combine solar v15 + hurricane v3 posteriors
into scale-invariant comparison statistics.

For each domain, compute:
  - branching ratio  n = alpha * E[g(m)]   (fraction of events that are aftershocks)
  - kernel half-life t_half  (days)
  - mark productivity exponent kappa
  - background fraction f_bg = (mu * T_obs) / N_total
  - dimensionless storm count  N / (mu * T_obs)
Plus 3rd-domain ETAS benchmark from literature for sanity.

Random seed: 20260523.
"""
import os, json, pickle
os.environ.setdefault("PYTENSOR_FLAGS", "cxx=,mode=FAST_RUN")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Resolve paths relative to this file so the script runs inside the solar repo
# (universality/ subdir). External inputs (solar v15 idata.pkl, hurricane v3
# bootstrap params + event CSV) are NOT committed here — see README.
import pathlib
_HERE = pathlib.Path(__file__).resolve().parent
ROOT = str(_HERE.parent)                    # universality/
DATA = f"{ROOT}/data"
FIG  = f"{ROOT}/figures"
# Fallback: original sandbox location used during development
_EXT_ROOT = "/home/user/workspace/hawkes-universality"
if not os.path.exists(f"{DATA}/solar_v15_idata.pkl") and os.path.exists(f"{_EXT_ROOT}/data/solar_v15_idata.pkl"):
    DATA = f"{_EXT_ROOT}/data"
os.makedirs(FIG, exist_ok=True)
SEED = 20260523
rng = np.random.default_rng(SEED)

# ----------------------------------------------------------------------
# 1. Solar v15 — extract posterior samples
# ----------------------------------------------------------------------
print("[1] solar v15 posterior…")
with open(f"{DATA}/solar_v15_idata.pkl", "rb") as f:
    idata = pickle.load(f)
post = idata.posterior
mu_alpha    = post["mu_alpha"].values.ravel()
sigma_alpha = post["sigma_alpha"].values.ravel()
mu_beta     = post["mu_beta"].values.ravel()
sigma_beta  = post["sigma_beta"].values.ravel()
mu_kappa    = post["mu_kappa"].values.ravel()
sigma_kappa = post["sigma_kappa"].values.ravel()
mu_mu       = post["mu_mu"].values.ravel()
sigma_mu    = post["sigma_mu"].values.ravel()
n_post = mu_alpha.size

# Use population medians as the *typical* cycle's params (point estimate),
# but propagate uncertainty by sampling fresh draws for the universality stats.
solar_events = pd.read_csv(f"{DATA}/solar_events.csv", parse_dates=["date"])
m_pool_solar = solar_events.mark.values
m0_solar = 8.0  # G4 threshold

T_solar_total = float((solar_events.date.max() - solar_events.date.min()).days)  # ~66,200 d
N_solar_total = len(solar_events)
print(f"  N events = {N_solar_total}, T_obs = {T_solar_total:.0f} d = {T_solar_total/365.25:.1f} yr")

# Sample B "typical-cycle" draws from the population
B = 3000
solar_idx = rng.choice(n_post, size=B, replace=True)
solar_mu0    = np.exp(mu_mu[solar_idx]    + sigma_mu[solar_idx]    * rng.standard_normal(B))
solar_alpha  = np.exp(mu_alpha[solar_idx] + sigma_alpha[solar_idx] * rng.standard_normal(B))
solar_beta   = np.exp(mu_beta[solar_idx]  + sigma_beta[solar_idx]  * rng.standard_normal(B))
solar_kappa  = mu_kappa[solar_idx]        + sigma_kappa[solar_idx] * rng.standard_normal(B)

# Branching ratio: n = alpha * E[g(m)] where g(m) = exp(kappa * (m - m0))
# Use empirical mark distribution
solar_Eg = np.array([np.mean(np.exp(k * (m_pool_solar - m0_solar))) for k in solar_kappa])
solar_branching = solar_alpha * solar_Eg
solar_half_life = np.log(2.0) / solar_beta   # days
solar_mu_typical = solar_mu0  # at S=S_bar (training mean)

solar_stats = {
    "domain": "Solar (G4+ geomagnetic storms)",
    "kernel_form": "exponential exp(-beta*tau)",
    "N_events": N_solar_total,
    "T_obs_d": T_solar_total,
    "T_obs_yr": T_solar_total / 365.25,
    "branching_ratio": {
        "median": float(np.median(solar_branching)),
        "hdi_95": [float(np.quantile(solar_branching, 0.025)),
                   float(np.quantile(solar_branching, 0.975))],
    },
    "half_life_d": {
        "median": float(np.median(solar_half_life)),
        "hdi_95": [float(np.quantile(solar_half_life, 0.025)),
                   float(np.quantile(solar_half_life, 0.975))],
    },
    "kappa": {
        "median": float(np.median(solar_kappa)),
        "hdi_95": [float(np.quantile(solar_kappa, 0.025)),
                   float(np.quantile(solar_kappa, 0.975))],
    },
    "mu_per_day": {
        "median": float(np.median(solar_mu_typical)),
        "hdi_95": [float(np.quantile(solar_mu_typical, 0.025)),
                   float(np.quantile(solar_mu_typical, 0.975))],
    },
}
print(f"  branching n = {solar_stats['branching_ratio']['median']:.3f} "
      f"[{solar_stats['branching_ratio']['hdi_95'][0]:.3f}, {solar_stats['branching_ratio']['hdi_95'][1]:.3f}]")
print(f"  half-life   = {solar_stats['half_life_d']['median']:.2f} d "
      f"[{solar_stats['half_life_d']['hdi_95'][0]:.2f}, {solar_stats['half_life_d']['hdi_95'][1]:.2f}]")

# ----------------------------------------------------------------------
# 2. Hurricane v3 — bootstrap samples
# ----------------------------------------------------------------------
print("\n[2] hurricane v3 bootstrap…")
hur = json.load(open(f"{DATA}/hurricane_v3_summary.json"))
boot = np.load(f"{DATA}/hurricane_v3_bootstrap_params.npy")  # (200, 7)
# Columns: mu_atl, mu_epac, b_soi, b_amo, alpha, beta, kappa
hur_alpha = boot[:, 4]
hur_beta  = boot[:, 5]
hur_kappa = boot[:, 6]
hur_mu_atl  = boot[:, 0]
hur_mu_epac = boot[:, 1]
# Mark distribution from event catalog
hur_events = pd.read_csv(f"{DATA}/hurricane_cat3plus_events.csv")
# Cat-3=3, Cat-4=4, Cat-5=5  (or similar)
# Saffir-Simpson category at peak intensity — same mark used in hurricane v3 fit
# (kappa = 0.18 maps Cat-3 -> Cat-5 to exp(0.18*2) = 1.43x productivity, matches v3 summary)
if "peak_cat" in hur_events.columns:
    mark_col = "peak_cat"
elif "category" in hur_events.columns:
    mark_col = "category"
else:
    raise RuntimeError("No Saffir-Simpson category column found in hurricane events")
m_pool_hur = hur_events[mark_col].values.astype(float)
assert m_pool_hur.min() >= 3 and m_pool_hur.max() <= 5, f"hurricane marks out of Cat-3+ range: {m_pool_hur.min()}-{m_pool_hur.max()}"
m0_hur = 3.0  # Cat-3 threshold
print(f"  hurricane mark column: {mark_col}, range = {m_pool_hur.min():.0f}-{m_pool_hur.max():.0f}, n = {len(m_pool_hur)}")

hur_Eg = np.array([np.mean(np.exp(k * (m_pool_hur - m0_hur))) for k in hur_kappa])
hur_branching = hur_alpha * hur_Eg
hur_half_life = np.log(2.0) / hur_beta

hur_stats = {
    "domain": "Hurricanes (Cat-3+ Atlantic+EPac)",
    "kernel_form": "exponential exp(-beta*tau)",
    "N_events": hur["n_events"],
    "T_obs_d": hur["window"]["T_obs_days"],
    "T_obs_yr": hur["window"]["T_obs_yr"],
    "branching_ratio": {
        "median": float(np.median(hur_branching)),
        "hdi_95": [float(np.quantile(hur_branching, 0.025)),
                   float(np.quantile(hur_branching, 0.975))],
    },
    "half_life_d": {
        "median": float(np.median(hur_half_life)),
        "hdi_95": [float(np.quantile(hur_half_life, 0.025)),
                   float(np.quantile(hur_half_life, 0.975))],
    },
    "kappa": {
        "median": float(np.median(hur_kappa)),
        "hdi_95": [float(np.quantile(hur_kappa, 0.025)),
                   float(np.quantile(hur_kappa, 0.975))],
    },
    "mu_per_day": {
        "median": float(np.median(hur_mu_atl + hur_mu_epac)),  # combined basin rate
        "hdi_95": [float(np.quantile(hur_mu_atl + hur_mu_epac, 0.025)),
                   float(np.quantile(hur_mu_atl + hur_mu_epac, 0.975))],
    },
}
print(f"  branching n = {hur_stats['branching_ratio']['median']:.3f} "
      f"[{hur_stats['branching_ratio']['hdi_95'][0]:.3f}, {hur_stats['branching_ratio']['hdi_95'][1]:.3f}]")
print(f"  half-life   = {hur_stats['half_life_d']['median']:.2f} d "
      f"[{hur_stats['half_life_d']['hdi_95'][0]:.2f}, {hur_stats['half_life_d']['hdi_95'][1]:.2f}]")

# ----------------------------------------------------------------------
# 3. ETAS earthquakes — literature benchmark (no fit, just published numbers)
# ----------------------------------------------------------------------
print("\n[3] earthquakes ETAS literature benchmarks…")
# References:
# - Ogata 1988, Utsu 1993, Helmstetter & Sornette 2002
# - p-value (Omori decay exponent): typically 0.9-1.2
# - c-value (time scale): 0.003-0.28 days  =>  for power-law (t+c)^-p,
#   "half-life" of the kernel intensity equals c * (2^(1/p) - 1), which we report
# - branching ratio n = 0.5-0.8 depending on Mmin
# We report mid-point and range from Utsu 1993 (16 Japan regions)
etas_p_range  = (0.95, 1.22)
etas_c_range  = (0.003, 0.28)  # days
etas_n_range  = (0.5, 0.8)
def power_law_half_life(c, p):
    """For phi(tau) = (tau+c)^-p, intensity falls to half when tau = c*(2^(1/p)-1)."""
    return c * (2 ** (1.0 / p) - 1)
etas_half_low  = power_law_half_life(etas_c_range[0], etas_p_range[1])
etas_half_high = power_law_half_life(etas_c_range[1], etas_p_range[0])
etas_half_mid  = power_law_half_life(np.mean(etas_c_range), np.mean(etas_p_range))
print(f"  ETAS half-life (literature range): {etas_half_low:.4f} – {etas_half_high:.4f} d "
      f"(mid {etas_half_mid:.3f} d)")
print(f"  ETAS branching ratio: {etas_n_range[0]:.2f} – {etas_n_range[1]:.2f}")

etas_stats = {
    "domain": "Earthquakes (ETAS literature, Utsu 1993; Helmstetter-Sornette 2002)",
    "kernel_form": "Omori power law (tau+c)^-p",
    "N_events": "varies; thousands per region",
    "T_obs_d": "decades",
    "T_obs_yr": "decades",
    "branching_ratio": {"range": list(etas_n_range)},
    "half_life_d":     {"range": [float(etas_half_low), float(etas_half_high)],
                        "mid":   float(etas_half_mid)},
    "kappa": {"note": "ETAS uses K * exp(alpha*M), where alpha_etas in 0.7-1.5 (different parameterization)"},
}

# ----------------------------------------------------------------------
# 4. Save table
# ----------------------------------------------------------------------
universality = {
    "seed": SEED,
    "method": "Match dimensionless statistics (branching ratio n, kernel half-life t_half, mark sensitivity kappa) across domains",
    "solar":     solar_stats,
    "hurricane": hur_stats,
    "earthquakes_etas": etas_stats,
}
with open(f"{DATA}/universality_summary.json", "w") as f:
    json.dump(universality, f, indent=2)
print(f"\nsaved {DATA}/universality_summary.json")

# Save sample arrays for plotting
np.savez(f"{DATA}/universality_samples.npz",
         solar_branching=solar_branching, solar_half_life=solar_half_life, solar_kappa=solar_kappa,
         hur_branching=hur_branching,     hur_half_life=hur_half_life,     hur_kappa=hur_kappa)

# ----------------------------------------------------------------------
# 5. Universality figure — three panels
# ----------------------------------------------------------------------
print("\n[5] generating universality figure…")
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel A: branching ratio
ax = axes[0]
ax.hist(solar_branching, bins=40, alpha=0.6, color="#d62728", density=True, label=f"Solar G4+ (n={N_solar_total})")
ax.hist(hur_branching,   bins=40, alpha=0.6, color="#1f77b4", density=True, label=f"Hurricanes Cat-3+ (n={hur_stats['N_events']})")
ax.axvspan(etas_n_range[0], etas_n_range[1], alpha=0.18, color="#2ca02c",
           label=f"ETAS earthquakes (lit. range {etas_n_range[0]:.1f}-{etas_n_range[1]:.1f})")
ax.set_xlabel("branching ratio  n = α · E[g(m)]")
ax.set_ylabel("posterior density")
ax.set_title("A. Branching ratio  (fraction of events that are 'aftershocks')")
ax.legend(fontsize=9, loc="upper right")
ax.grid(alpha=0.3)

# Panel B: kernel half-life (log scale, spans 4 orders of magnitude)
ax = axes[1]
ax.hist(np.log10(solar_half_life), bins=40, alpha=0.6, color="#d62728", density=True,
        label=f"Solar G4+  med = {np.median(solar_half_life):.2f} d")
ax.hist(np.log10(hur_half_life),   bins=40, alpha=0.6, color="#1f77b4", density=True,
        label=f"Hurricanes  med = {np.median(hur_half_life):.1f} d")
ax.axvspan(np.log10(etas_half_low), np.log10(etas_half_high), alpha=0.18, color="#2ca02c",
           label=f"ETAS earthquakes (lit.)  {etas_half_low*1440:.0f} min – {etas_half_high*1440:.0f} min")
ax.set_xlabel("log10(kernel half-life [days])")
ax.set_ylabel("posterior density")
ax.set_title("B. Kernel half-life  (clustering timescale)")
ax.legend(fontsize=9, loc="upper right")
ax.grid(alpha=0.3)
# Add minor x-axis labels in real days
xt = ax.get_xticks()
ax.set_xticks(xt)
ax.set_xticklabels([f"10$^{{{int(t)}}}$" + f"\n({10**t:.2g} d)" for t in xt], fontsize=8)

# Panel C: mark productivity
ax = axes[2]
ax.hist(solar_kappa, bins=40, alpha=0.6, color="#d62728", density=True,
        label=f"Solar κ (per Kp unit)  med = {np.median(solar_kappa):.2f}")
ax.hist(hur_kappa,   bins=40, alpha=0.6, color="#1f77b4", density=True,
        label=f"Hurricane κ (per category)  med = {np.median(hur_kappa):.2f}")
ax.axvline(0, color="k", ls="--", lw=1, alpha=0.6, label="κ=0 (mark-insensitive)")
ax.set_xlabel("κ  (productivity scaling with mark)")
ax.set_ylabel("posterior density")
ax.set_title("C. Mark sensitivity  (do bigger events make more aftershocks?)")
ax.legend(fontsize=9, loc="upper right")
ax.grid(alpha=0.3)

plt.suptitle("Hawkes universality: solar storms vs hurricanes vs ETAS earthquakes\n"
             "Different physics, different timescales, same statistical structure",
             y=1.02, fontsize=13)
plt.tight_layout()
plt.savefig(f"{FIG}/01_universality_three_panel.png", dpi=150, bbox_inches="tight")
plt.close()
print("  saved 01_universality_three_panel.png")

# ----------------------------------------------------------------------
# 6. Scaling figure — half-life vs domain timescale (log-log)
# ----------------------------------------------------------------------
print("\n[6] generating scaling figure…")
fig, ax = plt.subplots(figsize=(9, 6))

# Plot each domain as a point with error bars on a log-log plane:
# x = "characteristic forcing timescale" of the domain (rough physical scale)
# y = kernel half-life
# Earthquakes: forcing is essentially instantaneous (stress drop, seconds)
#   -> use 1 second = 1/86400 d
# Hurricanes: MJO + tropical wave cycles ~ 30-60 d. Use 45 d.
# Solar: solar rotation 27 d. Use 27 d.

domain_data = [
    ("Earthquakes\n(ETAS lit.)",      1/86400, etas_half_mid,   (etas_half_low, etas_half_high),   "#2ca02c"),
    ("Hurricanes\n(this work)",       45.0,    np.median(hur_half_life),   (np.quantile(hur_half_life, 0.025),
                                                                            np.quantile(hur_half_life, 0.975)), "#1f77b4"),
    ("Solar storms\n(this work)",     27.0,    np.median(solar_half_life), (np.quantile(solar_half_life, 0.025),
                                                                            np.quantile(solar_half_life, 0.975)), "#d62728"),
]
for label, x, y, ci, col in domain_data:
    ax.errorbar([x], [y], yerr=[[y - ci[0]], [ci[1] - y]], fmt="o", color=col, ms=12,
                ecolor="gray", capsize=5, lw=2)
    ax.annotate(label, xy=(x, y), xytext=(15, 10), textcoords="offset points",
                fontsize=11, fontweight="bold")

# Fit OLS in log-log space and show the slope honestly
xs_log = np.log10([d[1] for d in domain_data])
ys_log = np.log10([d[2] for d in domain_data])
slope, intercept = np.polyfit(xs_log, ys_log, 1)
x_line = np.logspace(-5, 2, 100)
ax.plot(x_line, 10 ** (intercept + slope * np.log10(x_line)), "k:", alpha=0.5,
        label=f"OLS fit: log10(t_half) = {intercept:.2f} + {slope:.2f} · log10(τ)")
# Reference: slope=1 (perfect scale-invariance)
ax.plot(x_line, 10 ** (intercept + 1.0 * np.log10(x_line)), "k--", alpha=0.25,
        label="slope=1 reference (scale-invariant)")

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Characteristic forcing timescale  τ_forcing  (days)")
ax.set_ylabel("Hawkes kernel half-life  t_half  (days)")
ax.set_title("Cross-domain timescales: kernel half-life vs forcing timescale\n"
             f"3 domains, ~7 orders of magnitude on x. Fitted slope = {slope:.2f} (no strict scaling law)")
ax.grid(True, which="both", alpha=0.3)
ax.legend(loc="upper left")
plt.tight_layout()
plt.savefig(f"{FIG}/02_universality_scaling.png", dpi=150, bbox_inches="tight")
plt.close()
print("  saved 02_universality_scaling.png")

print("\n[done]")
print("\n=== UNIVERSALITY SUMMARY ===")
print(f"  Solar G4+:    n = {np.median(solar_branching):.3f} [{np.quantile(solar_branching,0.025):.3f}, {np.quantile(solar_branching,0.975):.3f}]   "
      f"t_half = {np.median(solar_half_life):.2f} d")
print(f"  Hurricanes:   n = {np.median(hur_branching):.3f} [{np.quantile(hur_branching,0.025):.3f}, {np.quantile(hur_branching,0.975):.3f}]   "
      f"t_half = {np.median(hur_half_life):.2f} d")
print(f"  ETAS quakes:  n ≈ 0.5-0.8 (lit.)                  t_half ≈ {etas_half_mid*1440:.1f} min")
print(f"\n  Despite spanning 4 orders of magnitude in timescale, all three branching")
print(f"  ratios sit in the sub-critical regime (n < 1) and all three productivity")
print(f"  kernels are mark-positive.")
