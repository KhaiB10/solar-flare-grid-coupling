"""Quantify the log-log scaling: t_half vs tau_forcing across domains."""
import json
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent  # universality/
summary = json.load(open(ROOT / "data" / "universality_summary.json"))

# Characteristic forcing timescales (days)
# Earthquakes: S-wave traversal of a fault patch ~ seconds = 1.16e-5 d
# Solar: Carrington rotation 27.3 d
# Hurricanes: MJO ~ 45 d (also tropical wave 3-5 d) -- use MJO
tau_forcing = {
    "earthquakes": 10.0 / 86400,           # 10 s in days
    "solar":       27.3,                   # Carrington rotation
    "hurricane":   45.0,                   # MJO oscillation
}

t_half = {
    "earthquakes": summary["earthquakes_etas"]["half_life_d"]["mid"],
    "solar":       summary["solar"]["half_life_d"]["median"],
    "hurricane":   summary["hurricane"]["half_life_d"]["median"],
}

domains = ["earthquakes", "solar", "hurricane"]
x = np.array([np.log10(tau_forcing[d]) for d in domains])
y = np.array([np.log10(t_half[d]) for d in domains])

# Simple OLS in log-log space: log10(t_half) = a + b*log10(tau)
b, a = np.polyfit(x, y, 1)
y_hat = a + b * x
ss_res = np.sum((y - y_hat) ** 2)
ss_tot = np.sum((y - y.mean()) ** 2)
r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

# Implied ratio
ratios = {d: t_half[d] / tau_forcing[d] for d in domains}

print("=== Scaling test: log10(t_half) = a + b*log10(tau_forcing) ===\n")
for d in domains:
    print(f"  {d:14s}  tau={tau_forcing[d]:.3e} d   t_half={t_half[d]:.3e} d   ratio={ratios[d]:.3f}")

print(f"\n  slope  b = {b:.3f}   (b=1 means perfect scale-invariance)")
print(f"  inter. a = {a:.3f}   (10^a = {10**a:.3e}  prefactor)")
print(f"  R^2     = {r2:.4f}")
print(f"  ratio range across 3 domains: {min(ratios.values()):.3f} - {max(ratios.values()):.3f}")

out = {
    "tau_forcing_d": tau_forcing,
    "t_half_d": t_half,
    "log_log_slope": float(b),
    "log_log_intercept": float(a),
    "log_log_prefactor": float(10 ** a),
    "R2": float(r2),
    "ratio_t_half_over_tau": ratios,
    "interpretation": (
        f"Slope = {b:.2f} (b=1 = perfect scale-invariance). "
        f"R^2 = {r2:.3f} across 3 domains spanning 7 orders of magnitude on x. "
        f"All three domains have t_half between {min(ratios.values())*100:.1f}% and "
        f"{max(ratios.values())*100:.1f}% of their characteristic forcing timescale."
    ),
}
json.dump(out, open(ROOT / "data" / "scaling_test.json", "w"), indent=2)
print(f"\nsaved {ROOT/'data'/'scaling_test.json'}")
