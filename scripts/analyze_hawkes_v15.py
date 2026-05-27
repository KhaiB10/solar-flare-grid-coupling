#!/usr/bin/env python3
"""
v15: Hierarchical Bayesian marked Hawkes model with per-solar-cycle random
effects, fit with PyMC NUTS.

Structural change from v12: instead of one (mu0, alpha, beta, kappa) for the
whole 1844-2025 catalog, each solar cycle c gets its own draw from a
population distribution:

    mu0^(c)   ~ LogNormal(mu_mu, sigma_mu)
    alpha^(c) ~ LogNormal(mu_a,  sigma_a)
    beta^(c)  ~ LogNormal(mu_b,  sigma_b)
    kappa^(c) ~ Normal   (mu_k,  sigma_k)

gamma is pooled globally (it's a property of the F10.7 -> rate scaling, not
per-cycle physics). The exponential kernel is used (not Omori) for tractable
likelihood recursion inside PyTensor.

Likelihood per cycle:
    lambda^*_(c)(t) = mu0^(c) * (S(t)/S_bar)^gamma
                    + alpha^(c) * sum_{ti<t, cycle_i=c} g_i * exp(-beta^(c)*(t-ti))

Excitation is intra-cycle (storms in cycle c don't trigger storms in cycle c+1).
This is physical: the magnetosphere fully relaxes over the 11-year cycle gap.

Random seed: 20260523.
"""

import os, sys, time, json
# Disable the PyTensor C backend (linker bug with hidden symbols on this system’s libpython).
# The Python backend is fast enough for our model size (≤2 ms / gradient eval).
os.environ.setdefault("PYTENSOR_FLAGS", "cxx=,mode=FAST_RUN")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pytensor
import pytensor.tensor as pt
import pymc as pm
import arviz as az
from scipy import optimize

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FIG  = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

SEED = 20260523
rng = np.random.default_rng(SEED)

print(f"PyMC {pm.__version__}, PyTensor {pytensor.__version__}")

# ----------------------------------------------------------------------
# 1. Load data and define solar cycles
# ----------------------------------------------------------------------
print("\n[1] loading event catalog and daily F10.7 background…")
S_df = pd.read_csv(os.path.join(DATA, "derived_S_daily_v12.csv"), parse_dates=["date"])
S_df = S_df.sort_values("date").reset_index(drop=True)
T0 = S_df.date.iloc[0]; T_END = S_df.date.iloc[-1]
S_daily = S_df["S_daily"].values.astype(np.float64)
n_full = len(S_daily)
S_bar = float(np.mean(S_daily))
print(f"  S_daily {n_full} days, {T0.date()}→{T_END.date()}, S̄={S_bar:.2f} sfu")

ev = pd.read_csv(os.path.join(DATA, "derived_events_extended_1844_2025.csv"),
                 parse_dates=["date"])
ev = ev[(ev.date >= T0) & (ev.date <= T_END)].sort_values("date").reset_index(drop=True)
t_all = ((ev.date - T0).dt.days).values.astype(np.float64)
m_all = ev["mark"].values.astype(np.float64)
T_obs = float((T_END - T0).days)
N = len(t_all)
print(f"  {N} events in window, T_obs = {T_obs:.0f} d = {T_obs/365.25:.2f} yr")

# Solar cycle definitions (SIDC/NOAA standard minima)
CYCLES = [
    (8,  '1833-11-01', '1843-07-01'),
    (9,  '1843-07-01', '1855-12-01'),
    (10, '1855-12-01', '1867-03-01'),
    (11, '1867-03-01', '1878-12-01'),
    (12, '1878-12-01', '1890-03-01'),
    (13, '1890-03-01', '1902-01-01'),
    (14, '1902-01-01', '1913-07-01'),
    (15, '1913-07-01', '1923-08-01'),
    (16, '1923-08-01', '1933-09-01'),
    (17, '1933-09-01', '1944-02-01'),
    (18, '1944-02-01', '1954-04-01'),
    (19, '1954-04-01', '1964-10-01'),
    (20, '1964-10-01', '1976-03-01'),
    (21, '1976-03-01', '1986-09-01'),
    (22, '1986-09-01', '1996-08-01'),
    (23, '1996-08-01', '2008-12-01'),
    (24, '2008-12-01', '2019-12-01'),
    (25, '2019-12-01', '2031-07-01'),
]
# Window of each cycle in days-since-T0
def days_of(s): return float((pd.Timestamp(s) - T0).days)
cycle_starts = np.array([days_of(c[1]) for c in CYCLES])
cycle_ends   = np.array([days_of(c[2]) for c in CYCLES])
cycle_nums   = np.array([c[0] for c in CYCLES])
n_cycles = len(CYCLES)

# Assign each event to its cycle index
def cycle_idx(t):
    for i, (st, en) in enumerate(zip(cycle_starts, cycle_ends)):
        if st <= t < en:
            return i
    return -1

c_of_event = np.array([cycle_idx(t) for t in t_all], dtype=np.int32)
assert (c_of_event >= 0).all(), "Some events fell outside defined cycles"
print(f"\n  Solar cycles covered: SC{cycle_nums[0]}-SC{cycle_nums[-1]} ({n_cycles} cycles)")
print(f"  Events per cycle:")
unique, counts = np.unique(c_of_event, return_counts=True)
for i, c in zip(unique, counts):
    yr_s = pd.Timestamp(CYCLES[i][1]).year; yr_e = pd.Timestamp(CYCLES[i][2]).year
    print(f"    SC{cycle_nums[i]} ({yr_s}-{yr_e}): {c} events")

# Drop cycles with no in-window events from the fit (here all 17 have ≥3)
used_cycles = sorted(set(c_of_event.tolist()))
n_used = len(used_cycles)
print(f"\n  Cycles in fit: {n_used}")

# Restrict the observation window per cycle to its [start, end] intersected with [0, T_obs]
cycle_T_start = np.clip(cycle_starts[used_cycles], 0, T_obs)
cycle_T_end   = np.clip(cycle_ends[used_cycles],   0, T_obs)
cycle_T_dur   = cycle_T_end - cycle_T_start
print(f"  per-cycle observed durations (yr): "
      + ", ".join(f"SC{cycle_nums[i]}={d/365.25:.1f}" for i, d in zip(used_cycles, cycle_T_dur)))

# Remap c_of_event to dense indices in [0, n_used)
remap = {orig:new for new, orig in enumerate(used_cycles)}
c_idx = np.array([remap[c] for c in c_of_event], dtype=np.int32)

# Pre-compute per-event quantities for likelihood
m0 = 8.0
S_grid_t = np.arange(n_full, dtype=np.float64)

# Per-event background driver value S(t_i)
S_at_event = np.interp(t_all, S_grid_t, S_daily)

# For efficient excitation recursion per cycle, we need events sorted within each cycle
# (they already are, since the global event list is sorted)
# Precompute per-event Δt to the prior event in the SAME cycle (or -1 if first in cycle)
prev_in_cycle = -np.ones(N, dtype=np.int64)
for c in range(n_used):
    idx_in_c = np.where(c_idx == c)[0]
    for k in range(1, len(idx_in_c)):
        prev_in_cycle[idx_in_c[k]] = idx_in_c[k-1]
dt_to_prev = np.where(prev_in_cycle >= 0,
                      t_all - t_all[np.where(prev_in_cycle >= 0, prev_in_cycle, 0)],
                      0.0)

# Per-cycle background integral via cumulative S^gamma — easier to recompute inside likelihood
# We'll compute the background integrand as mu0_c * (S/S_bar)^gamma and integrate
# trapezoidally restricted to each cycle's [T_start, T_end].

# Pre-compute index ranges on the daily grid per cycle
cyc_lo = np.clip(np.floor(cycle_T_start).astype(int), 0, n_full)
cyc_hi = np.clip(np.floor(cycle_T_end).astype(int),   0, n_full)
# For each cycle, the daily S values restricted to the cycle window
print("\n[2] preparing PyTensor likelihood…")

# ----------------------------------------------------------------------
# PyTensor likelihood: vectorized across events using scan within each cycle.
# We'll write this as a numpy-style logp using pm.Potential to inject the
# Hawkes log-likelihood.
# ----------------------------------------------------------------------

# Convert to PyTensor shared variables
t_all_pt   = pt.as_tensor_variable(t_all)
m_all_pt   = pt.as_tensor_variable(m_all)
c_idx_pt   = pt.as_tensor_variable(c_idx, dtype="int32")
S_event_pt = pt.as_tensor_variable(S_at_event)
S_daily_pt = pt.as_tensor_variable(S_daily)
S_bar_pt   = pt.constant(S_bar)
S_grid_pt  = pt.as_tensor_variable(S_grid_t)

# Per-event excitation R_i requires recursion within each cycle. We implement
# this by precomputing in pure numpy at MLE time and inside PyTensor using a
# manual O(N) Python loop unrolled via aesara.scan with a "reset on cycle change" flag.
# However, scan with mutable state is finicky. We use a simpler approach:
# precompute (R_i) closed-form within each cycle as a function of per-event
# parameters using a Python-level scan applied symbolically.

# Precompute pairwise mask M_ij = 1 if j<i and j,i are in the same cycle, else 0.
# This lets us vectorize the Hawkes recursion: R_i = sum_{j<i, same cycle} g_j * exp(-beta_c * (t_i - t_j)).
# N=434, so this is a 434x434 matrix — trivially small.
ti = t_all[:, None]   # (N,1)
tj = t_all[None, :]   # (1,N)
dt_pair_np  = ti - tj                              # (N,N) Δt for each pair (>0 when j<i)
same_cyc_np = (c_idx[:, None] == c_idx[None, :])   # (N,N)
strict_lower = np.tri(N, N, k=-1, dtype=bool)      # j < i
mask_np      = (same_cyc_np & strict_lower).astype(np.float64)
# CRITICAL: clamp Δt to non-negative inside the exp domain. Upper-triangle entries
# have Δt<0 which would overflow exp(-β·Δt); even though we mask them out, 0×inf=NaN.
dt_pair_safe_np = np.where(mask_np > 0, dt_pair_np, 0.0)
dt_pair_pt   = pt.as_tensor_variable(dt_pair_safe_np)
mask_pt      = pt.as_tensor_variable(mask_np)

# Cycle-window cumulative-sum bounds (numpy constants for vectorized indexing)
cyc_lo_np = cyc_lo.astype(np.int64)
cyc_hi_np = cyc_hi.astype(np.int64)


def hawkes_loglik(mu0_vec, gamma, alpha_vec, beta_vec, kappa_vec):
    """
    mu0_vec, alpha_vec, beta_vec, kappa_vec: shape (n_used,) PyTensor tensors.
    gamma: scalar.
    Vectorized pairwise implementation (no scan) — gradients work cleanly.
    """
    # Per-event background rate
    log_S_ratio = pt.log(S_event_pt / S_bar_pt)
    mu_event    = mu0_vec[c_idx_pt] * pt.exp(gamma * log_S_ratio)

    # Per-event productivity g_i = exp(kappa_c * (m_i - m0))
    g_event     = pt.exp(kappa_vec[c_idx_pt] * (m_all_pt - m0))
    alpha_event = alpha_vec[c_idx_pt]      # (N,)
    beta_event  = beta_vec[c_idx_pt]       # (N,)

    # Pairwise self-excitation R_i = Σ_{j<i, same cycle} g_j * exp(-β_c · (t_i - t_j))
    # We use β at row i (β_i) since both i and j share the same cycle inside the mask.
    decay_pair = pt.exp(-beta_event[:, None] * dt_pair_pt)   # (N,N)
    contrib    = mask_pt * decay_pair * g_event[None, :]      # (N,N)
    R_event    = pt.sum(contrib, axis=1)                       # (N,)

    rate     = mu_event + alpha_event * R_event
    log_term = pt.sum(pt.log(rate + 1e-300))

    # Background integral per cycle via cumsum on the daily grid
    S_pow      = pt.exp(gamma * pt.log(S_daily_pt / S_bar_pt))   # (n_full,)
    cum_S_pow  = pt.concatenate([pt.zeros((1,), dtype="float64"),
                                  pt.cumsum(S_pow)])             # (n_full+1,) — shifted
    # sum over [lo, hi) = cum_S_pow[hi] - cum_S_pow[lo]
    cyc_lo_t   = pt.as_tensor_variable(cyc_lo_np)
    cyc_hi_t   = pt.as_tensor_variable(cyc_hi_np)
    bg_window  = cum_S_pow[cyc_hi_t] - cum_S_pow[cyc_lo_t]        # (n_used,)
    bg_total   = pt.sum(mu0_vec * bg_window)

    # Excitation compensator: per event α_c · g_i · (1 - exp(-β_c · (T_end_c - t_i))) / β_c
    T_end_per_event = pt.as_tensor_variable(cycle_T_end[c_idx])
    tau_remain      = T_end_per_event - t_all_pt
    exc_per_event   = alpha_event * g_event * (1.0 - pt.exp(-beta_event * tau_remain)) / beta_event
    exc_total       = pt.sum(exc_per_event)

    return log_term - bg_total - exc_total


# ----------------------------------------------------------------------
# 3. Build PyMC hierarchical model
# ----------------------------------------------------------------------
print("\n[3] building hierarchical model…")

# Warm-start values from v14 MLE (exp-kernel rescaled if needed)
# v12 (exp kernel) values: mu0=0.00423, gamma=2.18, alpha=0.0951, beta=0.583, kappa=1.08
v12_mu0 = 0.00423; v12_gamma = 2.18; v12_alpha = 0.0951; v12_beta = 0.583; v12_kappa = 1.08

with pm.Model() as model_v15:
    # Population (hyper) priors
    # Log-mu_mu centered at log(v12_mu0), wide sigma
    mu_mu      = pm.Normal("mu_mu",      mu=np.log(v12_mu0), sigma=1.0)
    sigma_mu   = pm.HalfNormal("sigma_mu", sigma=1.0)
    mu_alpha   = pm.Normal("mu_alpha",   mu=np.log(v12_alpha), sigma=1.0)
    sigma_alpha= pm.HalfNormal("sigma_alpha", sigma=1.0)
    mu_beta    = pm.Normal("mu_beta",    mu=np.log(v12_beta), sigma=0.5)
    sigma_beta = pm.HalfNormal("sigma_beta", sigma=0.5)
    mu_kappa   = pm.Normal("mu_kappa",   mu=v12_kappa, sigma=1.0)
    sigma_kappa= pm.HalfNormal("sigma_kappa", sigma=1.0)

    # Global parameter
    gamma = pm.Normal("gamma", mu=v12_gamma, sigma=1.0)

    # Per-cycle non-centered parameterization (helps NUTS sampling)
    z_mu     = pm.Normal("z_mu",    mu=0, sigma=1, shape=n_used)
    z_alpha  = pm.Normal("z_alpha", mu=0, sigma=1, shape=n_used)
    z_beta   = pm.Normal("z_beta",  mu=0, sigma=1, shape=n_used)
    z_kappa  = pm.Normal("z_kappa", mu=0, sigma=1, shape=n_used)

    log_mu0_c    = pm.Deterministic("log_mu0_c",    mu_mu    + sigma_mu    * z_mu)
    log_alpha_c  = pm.Deterministic("log_alpha_c",  mu_alpha + sigma_alpha * z_alpha)
    log_beta_c   = pm.Deterministic("log_beta_c",   mu_beta  + sigma_beta  * z_beta)
    kappa_c      = pm.Deterministic("kappa_c",      mu_kappa + sigma_kappa * z_kappa)

    mu0_c   = pm.Deterministic("mu0_c",   pt.exp(log_mu0_c))
    alpha_c = pm.Deterministic("alpha_c", pt.exp(log_alpha_c))
    beta_c  = pm.Deterministic("beta_c",  pt.exp(log_beta_c))

    # Hawkes likelihood via potential
    ll = hawkes_loglik(mu0_c, gamma, alpha_c, beta_c, kappa_c)
    pm.Potential("hawkes_loglik", ll)

print("  model built")
print("  parameters:", len(model_v15.free_RVs), "RVs")
for v in model_v15.free_RVs:
    print(f"    {v.name}: shape {v.type.shape}")

# ----------------------------------------------------------------------
# 4. Sample with NUTS
# ----------------------------------------------------------------------
print("\n[4] sampling with NUTS (4 chains × 1500 tune + 1500 draws)…")
t_sample = time.time()
with model_v15:
    idata = pm.sample(
        draws=1500, tune=1500, chains=4, cores=1,   # serial chains — Python backend
        target_accept=0.95,
        random_seed=SEED,
        progressbar=False,
        return_inferencedata=True,
    )
print(f"  sampling took {(time.time()-t_sample)/60:.1f} minutes")

# Save inference data (always pickle first — cheap insurance against backend bugs)
import pickle as _pkl
with open(os.path.join(DATA, "v15_idata.pkl"), "wb") as _f:
    _pkl.dump(idata, _f)
print("  trace saved to data/v15_idata.pkl")
try:
    idata.to_netcdf(os.path.join(DATA, "v15_idata.nc"), engine="h5netcdf")
    print("  trace also saved as data/v15_idata.nc")
except Exception as _e:
    print(f"  netcdf save skipped: {_e!r}")

# ----------------------------------------------------------------------
# 5. Diagnostics
# ----------------------------------------------------------------------
print("\n[5] convergence diagnostics…")
summ = az.summary(idata, var_names=["mu_mu","sigma_mu","mu_alpha","sigma_alpha",
                                     "mu_beta","sigma_beta","mu_kappa","sigma_kappa",
                                     "gamma"], hdi_prob=0.95)
print(summ.to_string())
max_rhat = float(summ["r_hat"].max())
min_ess = float(summ["ess_bulk"].min())
n_div = int(idata.sample_stats["diverging"].sum())
print(f"\n  max R-hat = {max_rhat:.3f}  (target < 1.01)")
print(f"  min ESS_bulk = {min_ess:.0f}  (target > 400)")
print(f"  divergent transitions = {n_div}")

# ----------------------------------------------------------------------
# 6. Per-cycle posterior summaries
# ----------------------------------------------------------------------
print("\n[6] per-cycle posterior summaries…")
per_cycle_summ = az.summary(idata, var_names=["mu0_c","alpha_c","beta_c","kappa_c"],
                            hdi_prob=0.95)
print(per_cycle_summ.to_string())
per_cycle_summ.to_csv(os.path.join(DATA, "v15_per_cycle_summary.csv"))

# SC24 outlier test: posterior P(mu0_SC24 < mu_population)
sc24_idx_used = used_cycles.index(cycle_nums.tolist().index(24)) if 24 in cycle_nums.tolist() else None
# Robust: find by mapping cycle_nums[used_cycles[k]] == 24
sc24_pos = next((k for k, ci in enumerate(used_cycles) if cycle_nums[ci] == 24), None)
sc25_pos = next((k for k, ci in enumerate(used_cycles) if cycle_nums[ci] == 25), None)
print(f"  SC24 position in used_cycles: {sc24_pos}; SC25 position: {sc25_pos}")

# Extract posterior arrays
post = idata.posterior
mu0_post   = post["mu0_c"].values   # shape (chain, draw, n_used)
alpha_post = post["alpha_c"].values
beta_post  = post["beta_c"].values
kappa_post = post["kappa_c"].values
mu_mu_post = post["mu_mu"].values   # population log-mean
sigma_mu_post = post["sigma_mu"].values

if sc24_pos is not None:
    mu0_sc24 = mu0_post[..., sc24_pos].ravel()
    mu0_pop_med = np.exp(mu_mu_post.ravel())  # population median
    sc24_lt_pop = float(np.mean(mu0_sc24 < mu0_pop_med))
    sc24_z = (np.log(mu0_sc24) - mu_mu_post.ravel()) / sigma_mu_post.ravel()
    sc24_z_med = float(np.median(sc24_z))
    sc24_z_p025, sc24_z_p975 = np.quantile(sc24_z, [0.025, 0.975])
    print(f"\n  SC24 mu0 posterior z-score (cycles from population mean):")
    print(f"    median = {sc24_z_med:.2f}, 95% HDI = [{sc24_z_p025:.2f}, {sc24_z_p975:.2f}]")
    print(f"    P(SC24 mu0 < population median) = {sc24_lt_pop:.3f}")
    if abs(sc24_z_med) > 1.5:
        print(f"    → SC24 is a statistical outlier (|z| > 1.5)")
    else:
        print(f"    → SC24 is not a strong outlier in the hierarchical fit")

# ----------------------------------------------------------------------
# 7. SC25 posterior predictive simulation 2024-2030
# ----------------------------------------------------------------------
print("\n[7] SC25 posterior predictive simulation through 2030…")
T_FORECAST_END = pd.Timestamp("2030-12-31")
t_fc_end = float((T_FORECAST_END - T0).days)
# SC25 observed events so far (in catalog) - use as conditioning
sc25_observed_idx = np.where(c_of_event == cycle_nums.tolist().index(25))[0]
sc25_obs_events = [(t_all[i], m_all[i]) for i in sc25_observed_idx]
print(f"  SC25 observed so far: {len(sc25_obs_events)} events through {T_END.date()}")
# Forecast window: from T_END+1 day to T_FORECAST_END
t_start_forecast = T_obs + 1.0  # day after T_END
forecast_days = np.arange(int(t_start_forecast), int(t_fc_end) + 1)
print(f"  forecast window: {len(forecast_days)} days = {len(forecast_days)/365.25:.2f} years")

# For S(t) in forecast window: we don't have F10.7 data past T_END.
# Use a climatological projection from SC25's currently-observed S(t) trajectory,
# extrapolated by SC25 mean F10.7. Simple approach: use the mean of last 365 days of S_daily as a flat extrapolation.
recent_S = float(np.mean(S_daily[-365:]))
S_forecast = np.full(len(forecast_days), recent_S)
print(f"  forecast S(t) = recent 365-d mean = {recent_S:.2f} sfu (flat extrapolation)")

# Simulate from the posterior: for each posterior draw of (mu0_SC25, alpha_SC25,
# beta_SC25, kappa_SC25, gamma), simulate a Hawkes process on [T_END+1, T_FORECAST_END]
# conditioned on the SC25 events observed so far.
n_post_draws = mu0_post.shape[0] * mu0_post.shape[1]
mu0_sc25_post   = mu0_post[..., sc25_pos].ravel()
alpha_sc25_post = alpha_post[..., sc25_pos].ravel()
beta_sc25_post  = beta_post[..., sc25_pos].ravel()
kappa_sc25_post = kappa_post[..., sc25_pos].ravel()
gamma_post      = post["gamma"].values.ravel()

# Empirical mark distribution from observed SC25 events (fallback: full catalog SC23+SC24+SC25 marks)
sc_recent_idx = np.where(np.isin(c_of_event, [cycle_nums.tolist().index(23),
                                                cycle_nums.tolist().index(24),
                                                cycle_nums.tolist().index(25)]))[0]
sc_recent_marks = m_all[sc_recent_idx]
print(f"  empirical mark distribution from SC23+24+25: n={len(sc_recent_marks)}, "
      f"min/median/max = {sc_recent_marks.min():.1f}/{np.median(sc_recent_marks):.1f}/{sc_recent_marks.max():.1f}")

# Simulation: thinning algorithm
B_sims = 2000  # subsample posterior to 2000 sims
rng_sim = np.random.default_rng(SEED + 7)
sim_indices = rng_sim.choice(n_post_draws, size=B_sims, replace=False)

sc25_obs_t = np.array([e[0] for e in sc25_obs_events])
sc25_obs_m = np.array([e[1] for e in sc25_obs_events])

def simulate_hawkes_forward(mu0, gamma_g, alpha, beta, kappa,
                             obs_t, obs_m, t_start, t_end, S_forecast_vals,
                             forecast_day_arr, rng):
    """Ogata thinning simulation, conditioned on observed past events."""
    events_t = list(obs_t)
    events_m = list(obs_m)
    g_vals = [np.exp(kappa * (m - m0)) for m in obs_m]
    t = t_start
    # Cap to prevent runaway: max possible events
    MAX_EV = 500
    while t < t_end and len(events_t) < len(obs_t) + MAX_EV:
        # Upper bound on rate over the next interval (use current state + max S contribution)
        # Current background rate at t:
        # Find S at t (index into forecast)
        idx = min(int(t - t_start), len(S_forecast_vals) - 1)
        if idx < 0: idx = 0
        S_t = S_forecast_vals[idx]
        mu_t = mu0 * (S_t / S_bar) ** gamma_g
        # Excitation at t from all past events
        if events_t:
            dt = t - np.array(events_t)
            exc = alpha * np.sum(np.array(g_vals) * np.exp(-beta * dt))
        else:
            exc = 0.0
        lam_upper = mu_t + exc + 1e-9  # tiny safety
        # Sample next candidate time
        dt_next = rng.exponential(1.0 / lam_upper)
        t_cand = t + dt_next
        if t_cand >= t_end:
            break
        # Recompute rate at t_cand for thinning
        idx = min(int(t_cand - t_start), len(S_forecast_vals) - 1)
        if idx < 0: idx = 0
        S_tc = S_forecast_vals[idx]
        mu_tc = mu0 * (S_tc / S_bar) ** gamma_g
        if events_t:
            dt = t_cand - np.array(events_t)
            exc_tc = alpha * np.sum(np.array(g_vals) * np.exp(-beta * dt))
        else:
            exc_tc = 0.0
        lam_tc = mu_tc + exc_tc
        # Accept with probability lam_tc / lam_upper
        if rng.random() < lam_tc / lam_upper:
            # Draw a mark from empirical distribution
            m_new = float(rng.choice(sc_recent_marks))
            events_t.append(t_cand)
            events_m.append(m_new)
            g_vals.append(np.exp(kappa * (m_new - m0)))
        t = t_cand
    # Return only the FORWARD events (after t_start, which is past T_obs)
    return [(events_t[i], events_m[i]) for i in range(len(obs_t), len(events_t))]

print("  running posterior predictive simulations…")
forward_counts = np.zeros(B_sims, dtype=int)
all_forward_events = []
for k, idx in enumerate(sim_indices):
    fwd = simulate_hawkes_forward(
        mu0=mu0_sc25_post[idx], gamma_g=gamma_post[idx],
        alpha=alpha_sc25_post[idx], beta=beta_sc25_post[idx],
        kappa=kappa_sc25_post[idx],
        obs_t=sc25_obs_t, obs_m=sc25_obs_m,
        t_start=t_start_forecast, t_end=t_fc_end,
        S_forecast_vals=S_forecast, forecast_day_arr=forecast_days,
        rng=rng_sim,
    )
    forward_counts[k] = len(fwd)
    if k < 200:  # save first 200 simulation traces for plotting
        all_forward_events.append(fwd)
    if (k+1) % 500 == 0:
        print(f"    sim {k+1}/{B_sims}: median forecast G4+ count = {np.median(forward_counts[:k+1]):.1f}")

print(f"\n  G4+ count forecast 2025-06 → 2030-12 (post-SC25 observed):")
print(f"    median = {np.median(forward_counts):.1f}")
print(f"    50% HDI = [{np.quantile(forward_counts, 0.25):.1f}, "
      f"{np.quantile(forward_counts, 0.75):.1f}]")
print(f"    95% HDI = [{np.quantile(forward_counts, 0.025):.1f}, "
      f"{np.quantile(forward_counts, 0.975):.1f}]")
print(f"    P(more than 20 G4+ events) = {np.mean(forward_counts > 20):.3f}")
print(f"    P(at least one G5+ in window) = (mark > 9 in empirical) — see plots")

# G5 estimate: G5 corresponds to peak mark ≥ 9 (from solar definition; v6+ used)
# For each sim, count G5 events. Re-run a quick pass to extract.
sc25_g5_count = np.zeros(B_sims, dtype=int)
for k, sim_evs in enumerate(all_forward_events):
    sc25_g5_count[k] = sum(1 for _, m in sim_evs if m >= 9)
# Extend to all sims using the events stored — but we only kept 200. Use the 200 as estimate.
n_kept = min(200, len(all_forward_events))
p_any_g5 = float(np.mean(sc25_g5_count[:n_kept] >= 1))
print(f"    P(at least one G5 by 2030) ≈ {p_any_g5:.3f} (from {n_kept} retained sims)")

# ----------------------------------------------------------------------
# 8. Plots
# ----------------------------------------------------------------------
print("\n[8] generating plots…")

# F1: Per-cycle mu0 posterior with population reference
fig, ax = plt.subplots(figsize=(10, 5))
cycle_labels = [f"SC{cycle_nums[i]}" for i in used_cycles]
mu0_med = np.array([np.median(mu0_post[..., k]) for k in range(n_used)])
mu0_lo  = np.array([np.quantile(mu0_post[..., k], 0.025) for k in range(n_used)])
mu0_hi  = np.array([np.quantile(mu0_post[..., k], 0.975) for k in range(n_used)])
x = np.arange(n_used)
ax.errorbar(x, mu0_med, yerr=[mu0_med-mu0_lo, mu0_hi-mu0_med], fmt="o",
            color="C0", capsize=4, label="per-cycle posterior")
pop_med = float(np.median(np.exp(mu_mu_post.ravel())))
ax.axhline(pop_med, color="C3", ls="--", label=f"population median = {pop_med:.4f}")
ax.set_xticks(x); ax.set_xticklabels(cycle_labels, rotation=45)
ax.set_ylabel("μ₀ (G4+ events/day at S̄)")
ax.set_title("v15 per-cycle background rate μ₀ posterior (95% HDI)")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "50_v15_per_cycle_mu0.png"), dpi=150)
plt.close()

# F2: Per-cycle beta (kernel half-life)
fig, ax = plt.subplots(figsize=(10, 5))
beta_med = np.array([np.median(beta_post[..., k]) for k in range(n_used)])
hl_med = 1.0 / beta_med
hl_lo  = 1.0 / np.array([np.quantile(beta_post[..., k], 0.975) for k in range(n_used)])
hl_hi  = 1.0 / np.array([np.quantile(beta_post[..., k], 0.025) for k in range(n_used)])
ax.errorbar(x, hl_med, yerr=[hl_med-hl_lo, hl_hi-hl_med], fmt="o",
            color="C2", capsize=4, label="per-cycle 1/β")
ax.axhline(1/v12_beta, color="C3", ls="--", label=f"v12 pooled = {1/v12_beta:.2f} d")
ax.set_xticks(x); ax.set_xticklabels(cycle_labels, rotation=45)
ax.set_ylabel("kernel half-life 1/β (days)")
ax.set_title("v15 per-cycle excitation half-life posterior (95% HDI)")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "51_v15_per_cycle_halflife.png"), dpi=150)
plt.close()

# F3: SC25 forecast distribution
fig, ax = plt.subplots(figsize=(9, 5))
ax.hist(forward_counts, bins=np.arange(forward_counts.max()+2)-0.5,
        color="C3", alpha=0.7, edgecolor="white")
med = np.median(forward_counts)
ax.axvline(med, color="k", lw=2, label=f"median = {med:.0f}")
ax.axvline(np.quantile(forward_counts, 0.025), color="k", ls=":", label="95% HDI")
ax.axvline(np.quantile(forward_counts, 0.975), color="k", ls=":")
ax.set_xlabel("number of G4+ events 2025-06 → 2030-12 (forward only)")
ax.set_ylabel("posterior probability")
ax.set_title(f"v15 SC25 forward forecast: G4+ events through 2030\n"
             f"({B_sims} posterior predictive draws, conditioned on {len(sc25_obs_events)} observed SC25 events)")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "52_v15_sc25_forecast.png"), dpi=150)
plt.close()

# F4: Hyperposterior summaries
fig, axes = plt.subplots(2, 4, figsize=(14, 7))
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
    ax_i.set_title(f"{label}\nmed={np.median(vals):.3f}, 95% HDI=[{np.quantile(vals,0.025):.3f}, {np.quantile(vals,0.975):.3f}]",
                   fontsize=9)
    ax_i.grid(alpha=0.3)
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
    "n_params_total": int(4*n_used + 1 + 8),  # per-cycle + gamma + 8 hyperparams
    "sampling": {
        "tune": 1500, "draws": 1500, "chains": 4,
        "max_rhat": max_rhat, "min_ess_bulk": min_ess, "divergent": n_div,
        "wallclock_min": (time.time()-t_sample)/60,
    },
    "hyperposteriors": {},
    "per_cycle": {},
    "sc25_forecast": {
        "window_start": str((T_END + pd.Timedelta(days=1)).date()),
        "window_end": str(T_FORECAST_END.date()),
        "S_forecast_assumed_sfu": recent_S,
        "n_sc25_observed_to_date": len(sc25_obs_events),
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
# Hyper-posteriors
for name in ["mu_mu","sigma_mu","mu_alpha","sigma_alpha","mu_beta","sigma_beta",
             "mu_kappa","sigma_kappa","gamma"]:
    vals = post[name].values.ravel()
    summary["hyperposteriors"][name] = {
        "median": float(np.median(vals)),
        "hdi_95": [float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975))],
    }
# Per-cycle
for k, ci in enumerate(used_cycles):
    summary["per_cycle"][f"SC{cycle_nums[ci]}"] = {
        "n_events": int(np.sum(c_of_event == ci)),
        "mu0_median": float(np.median(mu0_post[..., k])),
        "mu0_hdi95":  [float(np.quantile(mu0_post[..., k], 0.025)),
                       float(np.quantile(mu0_post[..., k], 0.975))],
        "alpha_median": float(np.median(alpha_post[..., k])),
        "beta_median":  float(np.median(beta_post[..., k])),
        "halflife_days_median": float(1.0 / np.median(beta_post[..., k])),
        "kappa_median": float(np.median(kappa_post[..., k])),
    }
# SC24 anomaly result
if sc24_pos is not None:
    summary["sc24_anomaly_test"] = {
        "z_score_median": float(sc24_z_med),
        "z_score_hdi95": [float(sc24_z_p025), float(sc24_z_p975)],
        "P_mu0_below_population": float(sc24_lt_pop),
        "is_outlier_abs_z_over_1p5": bool(abs(sc24_z_med) > 1.5),
    }

with open(os.path.join(DATA, "v15_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("  saved data/v15_summary.json")
print("\n[done v15]")
