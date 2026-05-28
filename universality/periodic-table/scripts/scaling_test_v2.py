"""
Scaling test for periodic_table_v2.csv

Fits log10(t_half) = a + b * log10(n) by ordinary least squares,
plus a robust Theil-Sen estimator. Reports slope/intercept/R^2 and
saves a JSON snapshot for the FINDINGS doc.

Strict scale-invariance would predict slope ~= -1 (t_half ~ 1/n).
v1 result: slope ~0.24-0.29 (effectively no scaling).

Random seed: 20260523.
"""
import csv
import json
import math
import numpy as np
from pathlib import Path
from scipy import stats

np.random.seed(20260523)

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "data" / "periodic_table_v2.csv"
OUT_JSON = ROOT / "data" / "scaling_test_v2.json"


def load_pairs():
    rows = []
    with open(CSV_PATH) as f:
        for r in csv.DictReader(f):
            try:
                n = float(r["n_branching"]) if r["n_branching"] else None
            except ValueError:
                n = None
            try:
                t = float(r["t_half_s"]) if r["t_half_s"] else None
            except ValueError:
                t = None
            if n is not None and t is not None and n > 0 and t > 0:
                rows.append({
                    "row_id": r["row_id"],
                    "domain": r["domain"],
                    "subdomain": r["subdomain"],
                    "n": n,
                    "t_half_s": t,
                    "log_n": math.log10(n),
                    "log_t": math.log10(t),
                })
    return rows


def main():
    pairs = load_pairs()
    print(f"Loaded {len(pairs)} (n, t_half) pairs from v2 table")
    if len(pairs) < 3:
        print("Not enough points for a regression")
        return

    x = np.array([p["log_n"] for p in pairs])
    y = np.array([p["log_t"] for p in pairs])

    # OLS
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    # Theil-Sen (robust)
    ts_slope, ts_intercept, ts_lo, ts_hi = stats.theilslopes(y, x, 0.95)

    # Domain breakdown
    by_domain = {}
    for p in pairs:
        by_domain.setdefault(p["domain"], []).append(p)

    print(f"\n=== OLS ===")
    print(f"slope     = {slope:+.4f}  (strict scaling would predict -1.0)")
    print(f"intercept = {intercept:+.4f}")
    print(f"R^2       = {r_value**2:.4f}")
    print(f"p-value   = {p_value:.4f}")
    print(f"std_err   = {std_err:.4f}")

    print(f"\n=== Theil-Sen (robust) ===")
    print(f"slope     = {ts_slope:+.4f}  (95% CI: [{ts_lo:+.4f}, {ts_hi:+.4f}])")
    print(f"intercept = {ts_intercept:+.4f}")

    print(f"\n=== Domain breakdown (rows with both n & t_half) ===")
    for d, plist in sorted(by_domain.items()):
        ns = [p["n"] for p in plist]
        ts = [p["t_half_s"] for p in plist]
        print(f"  {d:20s}  N={len(plist):2d}  n=[{min(ns):.3f}, {max(ns):.3f}]  t_half=[{min(ts):.1f}, {max(ts):.1f}] s")

    print(f"\n=== All pairs ===")
    for p in pairs:
        print(f"  {p['row_id']:8s} {p['domain']:18s} {p['subdomain'][:30]:30s} n={p['n']:.3f}  t_half={p['t_half_s']:.2f}s")

    result = {
        "n_pairs": len(pairs),
        "ols": {
            "slope": float(slope),
            "intercept": float(intercept),
            "r_squared": float(r_value**2),
            "p_value": float(p_value),
            "std_err": float(std_err),
        },
        "theil_sen": {
            "slope": float(ts_slope),
            "intercept": float(ts_intercept),
            "slope_95ci_lo": float(ts_lo),
            "slope_95ci_hi": float(ts_hi),
        },
        "by_domain_counts": {d: len(plist) for d, plist in by_domain.items()},
        "interpretation": (
            "Strict scale-invariance under universal Hawkes would predict slope=-1 "
            "(faster systems = higher branching). Observed slope is the empirical scaling. "
            "v1 OLS slope was ~0.24-0.29 (effectively flat, falsifying strict scale-invariance)."
        ),
    }
    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
