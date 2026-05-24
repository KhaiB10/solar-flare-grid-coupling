#!/usr/bin/env python3
"""
v9 post-mortem: where does the reliability gap come from?
We confirmed δ ≈ +0.22 (insignificant). But v8 showed obs-frequency
≈ 95-100% in the 30-65% prediction bins. If smooth SSN scaling doesn't
fix it, *what* does?

This script localizes the gap by binning the test-window forecast days
by calendar period and computing local Brier/observed-rate within each
2024 quarter vs the rest.
"""
import os, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FIG  = os.path.join(ROOT, "figures")

# Reload v9 summary
with open(os.path.join(DATA, "v9_summary.json")) as f:
    v9 = json.load(f)

# Recompute the rolling forecast bookkeeping using the saved params
events = pd.read_csv(os.path.join(DATA, "derived_events_extended_1868_2025.csv"),
                     parse_dates=["date"])
kp = pd.read_csv(os.path.join(DATA, "derived_daily_with_phase.csv"), parse_dates=["date"])
T0 = pd.Timestamp("1868-01-01"); SPLIT = pd.Timestamp("2016-01-01")
T_END_days = (kp.date.max() - T0).days
SPLIT_days = (SPLIT - T0).days
rng = np.random.default_rng(20260523)
t_int = (events.date - T0).dt.days.values.astype(float)
jitter = rng.uniform(0.0, 1.0, size=len(t_int))
order = np.argsort(t_int + jitter)
t_all = (t_int + jitter)[order]; m_all = events.mark.values[order].astype(float)
train_mask = t_all < SPLIT_days
t_tr = t_all[train_mask]; m_tr = m_all[train_mask]
t_te = t_all[~train_mask]; m_te = m_all[~train_mask]
T_tr = float(SPLIT_days); T_te = float(T_END_days - SPLIT_days)
T_all_days = T_END_days

ssn = pd.read_csv(os.path.join(DATA, "SN_m_tot_V2.0.txt"),
                  sep=r"\s+", header=None,
                  names=["year","month","yfrac","sn","sd","n","prov"], engine="python")
ssn = ssn[["yfrac","sn"]].copy(); ssn.loc[ssn.sn<0,"sn"] = np.nan
ssn["sn_smooth"] = ssn.sn.rolling(window=13, center=True, min_periods=7).mean()
ssn = ssn.dropna().reset_index(drop=True)
T0_year = T0.year + (T0.month-1)/12.0
ssn_t = (ssn.yfrac.values - T0_year) * 365.25
ssn_s = ssn.sn_smooth.values
def S(td): return np.interp(td, ssn_t, ssn_s, left=ssn_s[0], right=ssn_s[-1])
S_bar_train = ssn_s[(ssn_t >= 0) & (ssn_t <= SPLIT_days)].mean()

p = v9["train_mle_1868_2015"]
mu0, gamma, alpha0, beta, kappa, delta = p["mu0"], p["gamma"], p["alpha0"], p["beta"], p["kappa"], p["delta"]
m0_const = 8.0

g_pre = np.exp(kappa*(m_tr - m0_const))
a_pre = alpha0 * (S(t_tr)/S_bar_train)**delta
g_te = np.exp(kappa*(m_te - m0_const))
a_te = alpha0 * (S(t_te)/S_bar_train)**delta

grid_te = np.arange(T_tr, T_tr + T_te + 1.0, 1.0)
mu_grid = mu0*(S(grid_te)/S_bar_train)**gamma
cum_mu = np.concatenate([[0.0], np.cumsum(0.5*(mu_grid[:-1]+mu_grid[1:])*np.diff(grid_te))])
def int_mu(tq): return np.interp(tq, grid_te, cum_mu)
def Lambda(tq):
    base = int_mu(tq)
    pre = float(np.sum(a_pre*g_pre*(np.exp(-beta*(T_tr - t_tr))-np.exp(-beta*(tq - t_tr)))))
    mask = t_te < tq
    te = float(np.sum(a_te[mask]*g_te[mask]*(1.0 - np.exp(-beta*(tq - t_te[mask])))))
    return base+pre+te

WIN = 30.0
forecast_days = np.arange(T_tr, T_tr + T_te - WIN + 1, 1.0)
p_pred = np.zeros(len(forecast_days)); obs = np.zeros(len(forecast_days), dtype=int)
for k, td in enumerate(forecast_days):
    p_pred[k] = 1.0 - np.exp(-(Lambda(td+WIN)-Lambda(td)))
    obs[k] = int(((t_te >= td) & (t_te < td + WIN)).any())

dates_fc = T0 + pd.to_timedelta(forecast_days, unit="D")
df = pd.DataFrame({"date": dates_fc, "p_pred": p_pred, "obs": obs})

# Slice by "is the next 30d in the 2024 cluster window?"
gannon_window  = (df.date >= "2024-04-15") & (df.date <= "2024-06-30")
oct2024_window = (df.date >= "2024-09-15") & (df.date <= "2024-11-15")
def stats_for(mask, name):
    d = df[mask]
    if len(d) == 0: return
    pm, om = d.p_pred.mean(), d.obs.mean()
    br = ((d.p_pred - d.obs)**2).mean()
    clim = om*(1-om)
    bss = 1 - br/clim if clim > 0 else float("nan")
    print(f"  {name:30s}  n={len(d):4d}  pred={pm:.3f}  obs={om:.3f}  Brier={br:.3f}  BSS={bss:+.3f}")

print("\n[postmortem] Forecast performance sliced by calendar window")
stats_for(np.ones(len(df), dtype=bool), "all test forecast days")
stats_for(gannon_window,                "May 2024 Gannon window (Apr15-Jun30)")
stats_for(oct2024_window,               "Oct 2024 cluster window (Sep15-Nov15)")
stats_for(~(gannon_window|oct2024_window), "all OTHER days")

# Pull both windows' predicted-vs-observed daily trace
fig, axes = plt.subplots(2, 1, figsize=(13, 7), sharex=False)
for ax, win, title in zip(axes, [gannon_window, oct2024_window],
                          ["May 2024 Gannon window", "Oct 2024 cluster window"]):
    d = df[win]
    ax.plot(d.date, d.p_pred, color="C2", lw=2, label="v9 30-day forecast")
    # mark observed events in this window
    ev_dates = [T0 + pd.to_timedelta(t, unit="D") for t in t_te]
    ev_marks = m_te
    for ed, em in zip(ev_dates, ev_marks):
        if d.date.min() <= ed <= d.date.max() + pd.Timedelta(days=30):
            c = "darkred" if em >= 9 else "C3"
            ax.axvline(ed, color=c, alpha=0.5, lw=1.5)
    ax.set_title(title + f"   |   pred mean = {d.p_pred.mean():.2f},  obs ≥1 freq = {d.obs.mean():.2f}")
    ax.set_ylabel("P(≥1 in 30d)")
    ax.set_ylim(0, 1.05); ax.grid(alpha=0.3); ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(FIG, "27_v9_postmortem_2024.png"), dpi=140)
plt.close()
print("\n[plot] figures/27_v9_postmortem_2024.png  saved")
