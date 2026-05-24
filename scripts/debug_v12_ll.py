#!/usr/bin/env python3
"""Quick sanity check of v12 loglike at sensible starting points."""
import os, sys, numpy as np, pandas as pd

DATA = "data"
T0 = pd.Timestamp("1844-07-01")
T_END = pd.Timestamp("2025-05-31")
full_dates = pd.date_range(T0, T_END, freq="D")
n_full = len(full_dates)

# SSN smoothed proxy
ssn = pd.read_csv(os.path.join(DATA, "SN_m_tot_V2.0.txt"), sep=r"\s+", header=None,
                  names=["year","month","yfrac","sn","sd","n","prov"], engine="python")
ssn = ssn[["yfrac","sn"]]
ssn.loc[ssn.sn < 0, "sn"] = np.nan
ssn["sn_smooth"] = ssn.sn.rolling(window=13, center=True, min_periods=7).mean()
ssn = ssn.dropna().reset_index(drop=True)
yfrac_full = full_dates.year + (full_dates.dayofyear - 1)/365.25
ssn_at_full = np.interp(yfrac_full, ssn.yfrac.values, ssn.sn_smooth.values)

# GFZ F10.7
gfz = pd.read_csv(os.path.join(DATA, "Kp_ap_Ap_SN_F107_since_1932.txt"),
                  sep=r"\s+", comment="#", header=None, engine="python")
gfz["date"] = pd.to_datetime(dict(year=gfz[0], month=gfz[1], day=gfz[2]))
gfz = gfz[["date", 25]].rename(columns={25:"F107"})
gfz.loc[gfz.F107 < 0, "F107"] = np.nan
gfz_aligned = gfz.set_index("date").reindex(full_dates).F107.values

# splice
slope, intercept = 0.6666, 61.547
S_proxy = intercept + slope * ssn_at_full
S_daily = np.where(full_dates < pd.Timestamp("1947-02-14"), S_proxy, gfz_aligned)
S_daily = np.where(np.isnan(S_daily), S_proxy, S_daily)
S_daily = np.minimum(S_daily, 300.0)
S_daily = pd.Series(S_daily).rolling(5, center=True, min_periods=1).median().values
S_bar = float(np.nanmean(S_daily))
print(f"S_bar={S_bar:.3f}, min={S_daily.min():.1f}, max={S_daily.max():.1f}")

# events
events = pd.read_csv(os.path.join(DATA, "derived_events_extended_1844_2025.csv"), parse_dates=["date"])
events = events.sort_values("date").reset_index(drop=True)
t = (events.date - T0).dt.days.values.astype(float)
m = events.mark.values.astype(float)
T_obs = float((T_END - T0).days)
print(f"N={len(t)}, T_obs={T_obs:.0f}d")

S_grid_t = np.arange(n_full, dtype=float)
m0 = 8.0

def ll(p):
    mu0, gamma, alpha, beta, kappa = p
    if mu0 <= 0 or alpha < 0 or beta <= 0: return -np.inf
    S_e = np.interp(t, S_grid_t, S_daily)
    mu_e = mu0 * (S_e / S_bar) ** gamma
    g = np.exp(kappa * (m - m0))
    R = 0; s_log = 0
    for i in range(len(t)):
        if i > 0: R = np.exp(-beta*(t[i]-t[i-1]))*(R + g[i-1])
        rate = mu_e[i] + alpha*R
        if rate <= 0: return -np.inf
        s_log += np.log(rate)
    mu_grid = mu0 * (S_daily / S_bar)**gamma
    s_int_mu = np.trapezoid(mu_grid, S_grid_t)
    s_comp = (alpha/beta) * np.sum(g * (1.0 - np.exp(-beta*(T_obs - t))))
    return s_log - s_int_mu - s_comp

# Evaluate at v11 fit (it was using smoothed SSN background; F10.7-equiv parameters
# should be very close since splice is monotone affine)
print("\nLogL at candidate starting points:")
for label, p in [
    ("v11 fit (warm)",   (0.00499, 1.012, 0.0962, 0.582, 1.092)),
    ("start 1",          (0.00500, 1.00, 0.10, 0.60, 1.00)),
    ("start 2 (gamma=1.2)",(0.00500, 1.20, 0.10, 0.60, 1.00)),
    ("start 3 (gamma=0.8)",(0.00400, 0.80, 0.12, 0.55, 1.10)),
    ("low mu, hi gamma", (0.00010, 2.50, 0.10, 0.60, 1.00)),
    ("no-trigger Poisson",(0.0066, 1.00, 0.0,  1.00, 0.0)),
]:
    print(f"  {label:30s}: {ll(p):+10.2f}")
