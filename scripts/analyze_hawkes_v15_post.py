#!/usr/bin/env python3
"""
v15 post-sampling pipeline: loads the pickled InferenceData from analyze_hawkes_v15.py
and runs diagnostics, per-cycle summaries, SC25 forecast simulation, plots, and the
summary JSON. Separated so we don't have to re-sample on cosmetic errors.

Random seed: 20260523.
"""
import os, sys, time, json, pickle
os.environ.setdefault("PYTENSOR_FLAGS", "cxx=,mode=FAST_RUN")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import arviz as az

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) if "__file__" in globals() \
       else "/home/user/workspace/solar-flare-grid-coupling"
DATA = os.path.join(ROOT, "data")
FIG  = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

SEED = 20260523

# ----------------------------------------------------------------------
# Re-derive the data setup (mirrors analyze_hawkes_v15.py up to event/cycle assignment)
# ----------------------------------------------------------------------
print("[1] re-loading data + cycle assignment…")
S_df  = pd.read_csv(os.path.join(DATA, "derived_S_daily_v12.csv"), parse_dates=["date"]).sort_values("date").reset_index(drop=True)
ev_df = pd.read_csv(os.path.join(DATA, "derived_events_extended_1844_2025.csv"), parse_dates=["date"]).sort_values("date").reset_index(drop=True)
T0 = S_df.date.iloc[0]; T_END = S_df.date.iloc[-1]
S_df["t"]  = (S_df.date - T0).dt.total_seconds() / 86400
n_full = len(S_df); S_daily = S_df["S_daily"].values; S_bar = float(S_daily.mean())
ev = ev_df[(ev_df.date >= T0) & (ev_df.date <= T_END)].sort_values("date").reset_index(drop=True)
ev["t"] = (ev.date - T0).dt.total_seconds() / 86400
t_all = ev.t.values; m_all = ev.mark.values; N = len(t_all); m0 = 8.0
T_obs = float((T_END - T0).days)

CYCLES = [
    (8,  '1833-11-01', '1843-07-01'),  (9,  '1843-07-01', '1855-12-01'),
    (10, '1855-12-01', '1867-03-01'),  (11, '1867-03-01', '1878-12-01'),
    (12, '1878-12-01', '1890-03-01'),  (13, '1890-03-01', '1902-01-01'),
    (14, '1902-01-01', '1913-07-01'),  (15, '1913-07-01', '1923-08-01'),
    (16, '1923-08-01', '1933-09-01'),  (17, '1933-09-01', '1944-02-01'),
    (18, '1944-02-01', '1954-04-01'),  (19, '1954-04-01', '1964-10-01'),
    (20, '1964-10-01', '1976-03-01'),  (21, '1976-03-01', '1986-09-01'),
    (22, '1986-09-01', '1996-08-01'),  (23, '1996-08-01', '2008-12-01'),
    (24, '2008-12-01', '2019-12-01'),  (25, '2019-12-01', '2031-07-01'),
]
def days_of(s): return float((pd.Timestamp(s) - T0).days)
cycle_starts = np.array([days_of(c[1]) for c in CYCLES])
cycle_ends   = np.array([days_of(c[2]) for c in CYCLES])
cycle_nums   = np.array([c[0] for c in CYCLES])
def ci(t):
    for i, (s, e) in enumerate(zip(cycle_starts, cycle_ends)):
        if s <= t < e: return i
    return -1
c_of_event = np.array([ci(t) for t in t_all], dtype=int)
used_cycles = sorted(set(c_of_event.tolist()))
n_used = len(used_cycles)
remap = {o: n for n, o in enumerate(used_cycles)}
c_idx = np.array([remap[c] for c in c_of_event])
print(f"  N={N}, n_used={n_used}, cycles {[int(cycle_nums[i]) for i in used_cycles]}")

# v12 warm-start (used by FINDINGS to compare)
v12_mu0 = 0.00423; v12_gamma = 2.18; v12_alpha = 0.0951; v12_beta = 0.583; v12_kappa = 1.08

# ----------------------------------------------------------------------
# Load posterior
# ----------------------------------------------------------------------
print("[2] loading pickled InferenceData…")
with open(os.path.join(DATA, "v15_idata.pkl"), "rb") as f:
    idata = pickle.load(f)
post = idata.posterior
print("  posterior variables:", list(post.data_vars))
n_chain = post.dims["chain"]; n_draw = post.dims["draw"]
print(f"  chains={n_chain}, draws={n_draw}")

# ----------------------------------------------------------------------
# 5. Convergence diagnostics
# ----------------------------------------------------------------------
print("\n[5] convergence diagnostics…")
top_vars = ["mu_mu","sigma_mu","mu_alpha","sigma_alpha",
            "mu_beta","sigma_beta","mu_kappa","sigma_kappa","gamma"]
summ = az.summary(idata, var_names=top_vars, ci_prob=0.95)
print(summ.to_string())
max_rhat = float(summ["r_hat"].max())
min_ess  = float(summ["ess_bulk"].min())
n_div    = int(idata.sample_stats["diverging"].sum())
print(f"\n  max R-hat       = {max_rhat:.3f}  (target < 1.01)")
print(f"  min ESS_bulk    = {min_ess:.0f}     (target > 400)")
print(f"  divergent steps = {n_div}")

# ----------------------------------------------------------------------
# 6. Per-cycle summaries + SC24 outlier test
# ----------------------------------------------------------------------
print("\n[6] per-cycle posterior summaries…")
per_cycle_summ = az.summary(idata, var_names=["mu0_c","alpha_c","beta_c","kappa_c"], ci_prob=0.95)
print(per_cycle_summ.to_string())
per_cycle_summ.to_csv(os.path.join(DATA, "v15_per_cycle_summary.csv"))

sc24_pos = next((k for k, ci_idx in enumerate(used_cycles) if cycle_nums[ci_idx] == 24), None)
sc25_pos = next((k for k, ci_idx in enumerate(used_cycles) if cycle_nums[ci_idx] == 25), None)
print(f"  SC24 dense-index: {sc24_pos}; SC25 dense-index: {sc25_pos}")

mu0_post   = post["mu0_c"].values
alpha_post = post["alpha_c"].values
beta_post  = post["beta_c"].values
kappa_post = post["kappa_c"].values
mu_mu_post    = post["mu_mu"].values
sigma_mu_post = post["sigma_mu"].values

sc24_z_med = None; sc24_z_p025 = None; sc24_z_p975 = None; sc24_lt_pop = None
if sc24_pos is not None:
    mu0_sc24    = mu0_post[..., sc24_pos].ravel()
    mu0_pop_med = np.exp(mu_mu_post.ravel())
    sc24_lt_pop = float(np.mean(mu0_sc24 < mu0_pop_med))
    sc24_z      = (np.log(mu0_sc24) - mu_mu_post.ravel()) / sigma_mu_post.ravel()
    sc24_z_med  = float(np.median(sc24_z))
    sc24_z_p025, sc24_z_p975 = [float(x) for x in np.quantile(sc24_z, [0.025, 0.975])]
    print(f"\n  SC24 mu0 posterior z-score (cycles from population mean):")
    print(f"    median = {sc24_z_med:.2f}, 95% HDI = [{sc24_z_p025:.2f}, {sc24_z_p975:.2f}]")
    print(f"    P(SC24 mu0 < population median) = {sc24_lt_pop:.3f}")
    if abs(sc24_z_med) > 1.5:
        print(f"    SC24 is a statistical outlier (|z| > 1.5)")
    else:
        print(f"    SC24 is not a strong outlier in the hierarchical fit")

# ----------------------------------------------------------------------
# 7. SC25 posterior predictive simulation through 2030
# ----------------------------------------------------------------------
print("\n[7] SC25 posterior predictive simulation through 2030…")
T_FORECAST_END = pd.Timestamp("2030-12-31")
t_fc_end       = float((T_FORECAST_END - T0).days)

sc25_observed_idx = np.where(c_of_event == cycle_nums.tolist().index(25))[0]
sc25_obs_t = t_all[sc25_observed_idx]
sc25_obs_m = m_all[sc25_observed_idx]
print(f"  SC25 observed so far: {len(sc25_obs_t)} events through {T_END.date()}")

t_start_forecast = T_obs + 1.0
forecast_days    = np.arange(int(t_start_forecast), int(t_fc_end) + 1)
print(f"  forecast window: {len(forecast_days)} days = {len(forecast_days)/365.25:.2f} years")

recent_S    = float(np.mean(S_daily[-365:]))
S_forecast  = np.full(len(forecast_days), recent_S)
print(f"  forecast S(t) = recent 365-d mean = {recent_S:.2f} sfu (flat extrapolation)")

mu0_sc25_post   = mu0_post[..., sc25_pos].ravel()
alpha_sc25_post = alpha_post[..., sc25_pos].ravel()
beta_sc25_post  = beta_post[..., sc25_pos].ravel()
kappa_sc25_post = kappa_post[..., sc25_pos].ravel()
gamma_post      = post["gamma"].values.ravel()
n_post_draws    = mu0_sc25_post.size

# Mark distribution: empirical from SC23+24+25
sc_recent_idx = np.where(np.isin(c_of_event, [cycle_nums.tolist().index(23),
                                               cycle_nums.tolist().index(24),
                                               cycle_nums.tolist().index(25)]))[0]
sc_recent_marks = m_all[sc_recent_idx]
print(f"  empirical mark distribution from SC23+24+25: n={len(sc_recent_marks)}, "
      f"min/median/max = {sc_recent_marks.min():.1f}/{np.median(sc_recent_marks):.1f}/{sc_recent_marks.max():.1f}")

B_sims     = 2000
rng_sim    = np.random.default_rng(SEED + 7)
sim_indices = rng_sim.choice(n_post_draws, size=B_sims, replace=False)

def simulate_hawkes_forward(mu0, gamma_g, alpha, beta, kappa,
                             obs_t, obs_m, t_start, t_end, S_forecast_vals, rng):
    events_t = list(obs_t)
    events_m = list(obs_m)
    g_vals   = [np.exp(kappa * (mm - m0)) for mm in obs_m]
    t = t_start
    MAX_EV = 500
    n_gen  = 0
    while t < t_end and n_gen < MAX_EV:
        idx = min(int(t - t_start), len(S_forecast_vals) - 1)
        if idx < 0: idx = 0
        S_t  = S_forecast_vals[idx]
        mu_t = mu0 * (S_t / S_bar) ** gamma_g
        if events_t:
            dt  = t - np.array(events_t)
            exc = alpha * np.sum(np.array(g_vals) * np.exp(-beta * dt))
        else:
            exc = 0.0
        lam_upper = mu_t + exc + 1e-9
        dt_next   = rng.exponential(1.0 / lam_upper)
        t_cand    = t + dt_next
        if t_cand >= t_end:
            break
        idx = min(int(t_cand - t_start), len(S_forecast_vals) - 1)
        if idx < 0: idx = 0
        S_tc  = S_forecast_vals[idx]
        mu_tc = mu0 * (S_tc / S_bar) ** gamma_g
        dt    = t_cand - np.array(events_t)
        exc_tc = alpha * np.sum(np.array(g_vals) * np.exp(-beta * dt))
        lam_tc = mu_tc + exc_tc
        if rng.random() < lam_tc / lam_upper:
            m_new = float(rng.choice(sc_recent_marks))
            events_t.append(t_cand)
            events_m.append(m_new)
            g_vals.append(np.exp(kappa * (m_new - m0)))
            n_gen += 1
        t = t_cand
    return [(events_t[i], events_m[i]) for i in range(len(obs_t), len(events_t))]

print("  running posterior predictive simulations…")
t_sim_start = time.time()
forward_counts = np.zeros(B_sims, dtype=int)
all_forward_events = []
for k, idx in enumerate(sim_indices):
    fwd = simulate_hawkes_forward(
        mu0=mu0_sc25_post[idx], gamma_g=gamma_post[idx],
        alpha=alpha_sc25_post[idx], beta=beta_sc25_post[idx],
        kappa=kappa_sc25_post[idx],
        obs_t=sc25_obs_t, obs_m=sc25_obs_m,
        t_start=t_start_forecast, t_end=t_fc_end,
        S_forecast_vals=S_forecast, rng=rng_sim,
    )
    forward_counts[k] = len(fwd)
    if k < 200:
        all_forward_events.append(fwd)
    if (k+1) % 500 == 0:
        print(f"    sim {k+1}/{B_sims}: median count = {np.median(forward_counts[:k+1]):.1f}  "
              f"({(time.time()-t_sim_start):.1f}s elapsed)")

print(f"\n  G4+ count forecast 2025-06 → 2030-12 (forward only):")
print(f"    median   = {np.median(forward_counts):.1f}")
print(f"    50% HDI  = [{np.quantile(forward_counts, 0.25):.1f}, {np.quantile(forward_counts, 0.75):.1f}]")
print(f"    95% HDI  = [{np.quantile(forward_counts, 0.025):.1f}, {np.quantile(forward_counts, 0.975):.1f}]")
print(f"    P(>20)   = {np.mean(forward_counts > 20):.3f}")
print(f"    P(>10)   = {np.mean(forward_counts > 10):.3f}")

# G5 estimate (from kept-200 traces)
sc25_g5_count = np.array([sum(1 for _, mm in evs if mm >= 9) for evs in all_forward_events])
p_any_g5 = float(np.mean(sc25_g5_count >= 1))
print(f"    P(≥1 G5 by 2030) ≈ {p_any_g5:.3f}  (from {len(all_forward_events)} kept sims)")

# ----------------------------------------------------------------------
# 8. Plots
# ----------------------------------------------------------------------
print("\n[8] generating plots…")

cycle_labels = [f"SC{cycle_nums[i]}" for i in used_cycles]
x = np.arange(n_used)

# F1: per-cycle mu0
fig, ax = plt.subplots(figsize=(11, 5))
mu0_med = np.array([np.median(mu0_post[..., k]) for k in range(n_used)])
mu0_lo  = np.array([np.quantile(mu0_post[..., k], 0.025) for k in range(n_used)])
mu0_hi  = np.array([np.quantile(mu0_post[..., k], 0.975) for k in range(n_used)])
colors = ["C3" if cycle_nums[used_cycles[k]] == 24 else
          "C2" if cycle_nums[used_cycles[k]] == 25 else "C0" for k in range(n_used)]
ax.errorbar(x, mu0_med, yerr=[mu0_med - mu0_lo, mu0_hi - mu0_med],
            fmt="o", ecolor="gray", capsize=4)
for xi, mi, col in zip(x, mu0_med, colors):
    ax.plot(xi, mi, "o", color=col, ms=8, zorder=3)
pop_med = float(np.median(np.exp(mu_mu_post.ravel())))
ax.axhline(pop_med, color="k", ls="--", lw=1, label=f"population median = {pop_med:.4f}")
ax.axhline(v12_mu0, color="C1", ls=":", lw=1, label=f"v12 pooled = {v12_mu0:.4f}")
ax.set_xticks(x); ax.set_xticklabels(cycle_labels, rotation=45)
ax.set_ylabel("μ₀  (G4+ events per day at S̄)")
ax.set_title("v15: per-cycle background rate μ₀ — posterior median + 95% HDI\n(red = SC24, green = SC25)")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "50_v15_per_cycle_mu0.png"), dpi=150)
plt.close()

# F2: per-cycle half-life 1/β
fig, ax = plt.subplots(figsize=(11, 5))
beta_med = np.array([np.median(beta_post[..., k]) for k in range(n_used)])
hl_med = 1.0 / beta_med
hl_lo  = 1.0 / np.array([np.quantile(beta_post[..., k], 0.975) for k in range(n_used)])
hl_hi  = 1.0 / np.array([np.quantile(beta_post[..., k], 0.025) for k in range(n_used)])
ax.errorbar(x, hl_med, yerr=[hl_med - hl_lo, hl_hi - hl_med], fmt="o", color="C2", capsize=4)
ax.axhline(1/v12_beta, color="C3", ls="--", label=f"v12 pooled = {1/v12_beta:.2f} d")
ax.set_xticks(x); ax.set_xticklabels(cycle_labels, rotation=45)
ax.set_ylabel("kernel half-life 1/β  (days)")
ax.set_title("v15: per-cycle excitation half-life — posterior median + 95% HDI")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "51_v15_per_cycle_halflife.png"), dpi=150)
plt.close()

# F3: SC25 forecast distribution
fig, ax = plt.subplots(figsize=(10, 5))
max_c = int(forward_counts.max())
ax.hist(forward_counts, bins=np.arange(max_c + 2) - 0.5,
        color="C3", alpha=0.7, edgecolor="white")
med = float(np.median(forward_counts))
q025, q25, q75, q975 = np.quantile(forward_counts, [0.025, 0.25, 0.75, 0.975])
ax.axvline(med, color="k", lw=2, label=f"median = {med:.0f}")
ax.axvline(q025, color="k", ls=":", label="95% HDI")
ax.axvline(q975, color="k", ls=":")
ax.axvspan(q25, q75, alpha=0.15, color="C3", label=f"50% HDI [{q25:.0f}, {q75:.0f}]")
ax.set_xlabel(f"number of G4+ events  {pd.Timedelta(days=1) + T_END:%Y-%m-%d} → 2030-12-31  (forward only)")
ax.set_ylabel("posterior probability (counts)")
ax.set_title(f"v15 SC25 forward forecast — {B_sims} posterior predictive draws\n"
             f"conditioned on {len(sc25_obs_t)} observed SC25 events")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "52_v15_sc25_forecast.png"), dpi=150)
plt.close()

# F4: hyperposteriors
fig, axes = plt.subplots(2, 4, figsize=(15, 7))
for ax_i, (name, label) in zip(axes.flat,
                                [("mu_mu","log μ₀ pop mean"),
                                 ("sigma_mu","log μ₀ pop sigma"),
                                 ("mu_alpha","log α pop mean"),
                                 ("sigma_alpha","log α pop sigma"),
                                 ("mu_beta","log β pop mean"),
                                 ("sigma_beta","log β pop sigma"),
                                 ("mu_kappa","κ pop mean"),
                                 ("sigma_kappa","κ pop sigma")]):
    vals = post[name].values.ravel()
    ax_i.hist(vals, bins=50, color="C0", alpha=0.7)
    ax_i.axvline(np.median(vals), color="k", lw=2)
    ax_i.set_title(f"{label}\nmed={np.median(vals):.3f},  95% HDI=[{np.quantile(vals,0.025):.3f}, {np.quantile(vals,0.975):.3f}]",
                   fontsize=9)
    ax_i.grid(alpha=0.3)
plt.suptitle("v15 hierarchical hyperposteriors", y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "53_v15_hyperposteriors.png"), dpi=150)
plt.close()
print("  saved figures 50-53")

# ----------------------------------------------------------------------
# 9. Summary JSON
# ----------------------------------------------------------------------
print("\n[9] writing summary…")
summary = {
    "seed": SEED,
    "model": "Hierarchical Bayes marked Hawkes (exp kernel), per-cycle (mu0, alpha, beta, kappa), pooled gamma",
    "n_events": int(N),
    "n_cycles": int(n_used),
    "cycles_modeled": [int(cycle_nums[i]) for i in used_cycles],
    "n_params_total": int(4 * n_used + 1 + 8),
    "sampling": {
        "tune": 1500, "draws": 1500, "chains": int(n_chain),
        "max_rhat": max_rhat, "min_ess_bulk": min_ess, "divergent": n_div,
    },
    "hyperposteriors": {},
    "per_cycle": {},
    "sc25_forecast": {
        "window_start": str((T_END + pd.Timedelta(days=1)).date()),
        "window_end":   str(T_FORECAST_END.date()),
        "S_forecast_assumed_sfu": recent_S,
        "n_sc25_observed_to_date": int(len(sc25_obs_t)),
        "median_forward_G4plus": int(np.median(forward_counts)),
        "hdi_50_forward_G4plus": [int(np.quantile(forward_counts, 0.25)),
                                   int(np.quantile(forward_counts, 0.75))],
        "hdi_95_forward_G4plus": [int(np.quantile(forward_counts, 0.025)),
                                   int(np.quantile(forward_counts, 0.975))],
        "P_more_than_20": float(np.mean(forward_counts > 20)),
        "P_more_than_10": float(np.mean(forward_counts > 10)),
        "P_any_G5_by_2030_approx": p_any_g5,
        "n_sims": B_sims,
    },
}
for name in ["mu_mu","sigma_mu","mu_alpha","sigma_alpha","mu_beta","sigma_beta",
             "mu_kappa","sigma_kappa","gamma"]:
    vals = post[name].values.ravel()
    summary["hyperposteriors"][name] = {
        "median": float(np.median(vals)),
        "hdi_95": [float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975))],
    }
for k, ci_idx in enumerate(used_cycles):
    summary["per_cycle"][f"SC{cycle_nums[ci_idx]}"] = {
        "n_events": int(np.sum(c_of_event == ci_idx)),
        "mu0_median": float(np.median(mu0_post[..., k])),
        "mu0_hdi95":  [float(np.quantile(mu0_post[..., k], 0.025)),
                       float(np.quantile(mu0_post[..., k], 0.975))],
        "alpha_median": float(np.median(alpha_post[..., k])),
        "beta_median":  float(np.median(beta_post[..., k])),
        "halflife_days_median": float(1.0 / np.median(beta_post[..., k])),
        "kappa_median": float(np.median(kappa_post[..., k])),
    }
if sc24_pos is not None:
    summary["sc24_anomaly_test"] = {
        "z_score_median": sc24_z_med,
        "z_score_hdi95":  [sc24_z_p025, sc24_z_p975],
        "P_mu0_below_population": sc24_lt_pop,
        "is_outlier_abs_z_over_1p5": bool(abs(sc24_z_med) > 1.5),
    }

with open(os.path.join(DATA, "v15_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("  saved data/v15_summary.json")
print("\n[done v15 post-processing]")
