#!/usr/bin/env python3
"""
v5: Non-stationary Hawkes-process fit for G4+ geomagnetic storms, 1932-2025.

Model: univariate exponential-kernel Hawkes with SSN-modulated background:

    λ(t) = μ(t) + Σ_{t_i < t}  α * exp(-β * (t - t_i))

    μ(t) = μ0 * ( S(t) / S_bar )^γ

where S(t) is the 13-month-smoothed monthly mean sunspot number from SILSO
and S_bar is its long-term mean (1749-present subset over 1932-2025).

Parameters:
    μ0 : reference background rate at S(t) = S_bar  (events / day)
    γ  : SSN-modulation exponent (γ = 0 ⇒ stationary background; γ > 0 ⇒ more
         events near solar max as expected)
    α  : excitation strength
    β  : decay rate (1/day)
    η  = α/β : branching ratio

Stationarity (subcriticality) requires η < 1; non-stationarity here refers
ONLY to the background μ(t), since the kernel is still exponential.

References (in addition to v4):
    Ogata (1988) JASA — non-stationary Hawkes ETAS for earthquakes
    Clauset & Woodard (2013) — heavy-tailed inter-arrivals
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize, stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FIG  = os.path.join(ROOT, "figures")

# ----------------------------------------------------------------------
# 1. Load G4+ event times and SSN(t)
# ----------------------------------------------------------------------
daily = pd.read_csv(os.path.join(DATA, "derived_daily_with_phase.csv"),
                    parse_dates=["date"])
daily = daily.sort_values("date").reset_index(drop=True)
T0 = daily.date.min()
T_end = (daily.date.max() - T0).days

g4 = daily[daily.Kp_max >= 8].sort_values("date").reset_index(drop=True)
rng = np.random.default_rng(20260523)
t_int = (g4.date - T0).dt.days.values.astype(float)
t = np.sort(t_int + rng.uniform(0.0, 1.0, size=len(t_int)))
N = len(t)
print(f"[load] N = {N} G4+ events, T = {T_end:,} days = {T_end/365.25:.1f} years")

# SILSO monthly SN
ssn = pd.read_csv(os.path.join(DATA, "SN_m_tot_V2.0.txt"),
                  sep=r"\s+", header=None,
                  names=["year", "month", "yfrac", "sn", "sd", "n", "prov"],
                  engine="python")
ssn = ssn[["yfrac", "sn"]].copy()
ssn.loc[ssn.sn < 0, "sn"] = np.nan
# 13-month centered smoother (standard SILSO convention)
ssn["sn_smooth"] = ssn.sn.rolling(window=13, center=True, min_periods=7).mean()
ssn = ssn.dropna(subset=["sn_smooth"]).reset_index(drop=True)
print(f"[ssn] {len(ssn)} smoothed monthly SSN values, "
      f"range {ssn.yfrac.min():.2f} – {ssn.yfrac.max():.2f}")

# Convert SSN yfrac to days since T0 (1932-01-01)
T0_year = T0.year + (T0.month - 1)/12.0
ssn_t = (ssn.yfrac.values - T0_year) * 365.25
ssn_s = ssn.sn_smooth.values

# Restrict to the 1932-2025 window and compute mean
mask = (ssn_t >= 0) & (ssn_t <= T_end)
S_bar = ssn_s[mask].mean()
print(f"[ssn] mean smoothed SSN over 1932-2025 = {S_bar:.2f}")

def S(t_days):
    """Smoothed SSN at time t (days since 1932-01-01), linear interpolation."""
    return np.interp(t_days, ssn_t, ssn_s, left=ssn_s[0], right=ssn_s[-1])

# Precompute S at event times and on a fine grid for integration of compensator
S_events = S(t)
# Daily grid for integrating μ(t) over [0, T_end]
grid = np.arange(0, T_end + 1.0, 1.0)
S_grid = S(grid)

# ----------------------------------------------------------------------
# 2. Log-likelihood for non-stationary Hawkes (exponential kernel)
#
#    μ(t) = μ0 * (S(t)/S_bar)^γ
#    ln L = Σ_i ln(μ(t_i) + α R_i)  -  ∫_0^T μ(t) dt
#                                    -  (α/β) Σ_i (1 - exp(-β(T - t_i)))
#
#    Excitation compensator unchanged from v4.
# ----------------------------------------------------------------------
def loglike(params, t, T, S_events, grid, S_grid, S_bar):
    mu0, gamma, alpha, beta = params
    if mu0 <= 0 or alpha < 0 or beta <= 0:
        return -1e18
    # μ(t_i)
    mu_events = mu0 * (S_events / S_bar)**gamma
    # Ogata recursion for kernel sum
    N = len(t)
    R = 0.0
    s_log = 0.0
    for i in range(N):
        if i > 0:
            dt = t[i] - t[i-1]
            R = np.exp(-beta*dt) * (1.0 + R)
        rate = mu_events[i] + alpha * R
        if rate <= 0:
            return -1e18
        s_log += np.log(rate)
    # ∫ μ(t) dt via trapezoidal rule on daily grid
    mu_grid = mu0 * (S_grid / S_bar)**gamma
    s_int_mu = np.trapezoid(mu_grid, grid)
    # Excitation compensator
    s_comp_exc = (alpha/beta) * np.sum(1.0 - np.exp(-beta*(T - t)))
    return s_log - s_int_mu - s_comp_exc

def neg_ll(params, *args):
    return -loglike(params, *args)

# ----------------------------------------------------------------------
# 3. Multi-start MLE for (μ0, γ, α, β)
# ----------------------------------------------------------------------
print("\n[fit] multi-start MLE for non-stationary Hawkes (μ0, γ, α, β)…")
best = None
starts = [
    (0.00513, 0.5, 0.1634, 0.5749),  # v4 best + γ=0.5
    (0.00513, 1.0, 0.1634, 0.5749),  # γ=1 (linear)
    (0.00513, 1.5, 0.1634, 0.5749),
    (0.004,   0.7, 0.10,   0.30),
    (0.006,   0.3, 0.20,   0.80),
    (0.005,   2.0, 0.15,   0.50),
    (0.003,   1.2, 0.18,   0.60),
    (0.007,   0.0, 0.15,   0.55),    # γ=0 ⇒ should recover v4
]
for trial, x0 in enumerate(starts):
    res = optimize.minimize(neg_ll, x0,
                            args=(t, T_end, S_events, grid, S_grid, S_bar),
                            method="Nelder-Mead",
                            options={"xatol":1e-7, "fatol":1e-7, "maxiter":40000})
    if res.fun < 1e17:
        mu0r, gr, ar, br = res.x
        print(f"  trial {trial+1}: μ0={mu0r:.5f}  γ={gr:.3f}  α={ar:.4f}  "
              f"β={br:.4f}  η={ar/br:.3f}  -LL={res.fun:.2f}")
        if best is None or res.fun < best.fun:
            best = res

mu0_hat, gamma_hat, alpha_hat, beta_hat = best.x
eta_hat = alpha_hat / beta_hat
ll_v5 = -best.fun
aic_v5 = 2*4 - 2*ll_v5
print(f"\n[MLE v5] best fit:")
print(f"  μ0 = {mu0_hat:.5f} events/day  =  {mu0_hat*365.25:.3f} events/year (background at S=S_bar)")
print(f"  γ  = {gamma_hat:.3f}  (SSN-modulation exponent)")
print(f"  α  = {alpha_hat:.4f}  (excitation amplitude)")
print(f"  β  = {beta_hat:.4f} 1/day,  decay timescale 1/β = {1/beta_hat:.2f} days")
print(f"  η  = α/β = {eta_hat:.3f}  (branching ratio)")
print(f"  log-likelihood = {ll_v5:.2f},  AIC = {aic_v5:.2f}")

# ----------------------------------------------------------------------
# 4. Compare to v4 stationary Hawkes and homogeneous Poisson
# ----------------------------------------------------------------------
print("\n[compare] refitting v4 stationary Hawkes (γ=0) for direct comparison…")
def loglike_stat(params, t, T):
    mu, alpha, beta = params
    if mu <= 0 or alpha < 0 or beta <= 0: return -1e18
    N = len(t); R = 0.0; s_log = 0.0
    for i in range(N):
        if i > 0:
            R = np.exp(-beta*(t[i]-t[i-1]))*(1.0+R)
        rate = mu + alpha*R
        if rate <= 0: return -1e18
        s_log += np.log(rate)
    s_comp = mu*T + (alpha/beta)*np.sum(1.0 - np.exp(-beta*(T - t)))
    return s_log - s_comp
res_v4 = optimize.minimize(lambda p: -loglike_stat(p, t, T_end),
                           [0.00513, 0.1634, 0.5749],
                           method="Nelder-Mead",
                           options={"xatol":1e-7,"fatol":1e-7,"maxiter":20000})
ll_v4 = -res_v4.fun
aic_v4 = 2*3 - 2*ll_v4

lam_pois = N / T_end
ll_pois = N*np.log(lam_pois) - lam_pois*T_end
aic_pois = 2*1 - 2*ll_pois

print(f"  v5 (non-stat Hawkes, 4 params):  LL={ll_v5:.2f}  AIC={aic_v5:.2f}")
print(f"  v4 (stationary Hawkes, 3 params): LL={ll_v4:.2f}  AIC={aic_v4:.2f}")
print(f"  v0 (Poisson, 1 param):            LL={ll_pois:.2f}  AIC={aic_pois:.2f}")
print(f"  ΔAIC (v5 − v4) = {aic_v5 - aic_v4:.2f}   (negative = v5 wins)")
print(f"  ΔAIC (v5 − Poisson) = {aic_v5 - aic_pois:.2f}")
# Likelihood-ratio test v5 vs v4 (γ = 0 is a boundary of the simpler model)
lr_v5_v4 = 2*(ll_v5 - ll_v4)
p_v5_v4 = 1 - stats.chi2.cdf(lr_v5_v4, df=1)
print(f"  LR (v5 vs v4) = {lr_v5_v4:.2f},  χ²(1) p = {p_v5_v4:.3e}")

# ----------------------------------------------------------------------
# 5. Goodness-of-fit: time-rescaling theorem for non-stationary case.
#    Λ(t_i) = ∫_0^{t_i} μ(s) ds  +  (α/β) Σ_{j<i}(1 − exp(−β(t_i − t_j)))
# ----------------------------------------------------------------------
print("\n[GOF] time-rescaling residuals…")
mu_grid_hat = mu0_hat * (S_grid / S_bar)**gamma_hat
# Cumulative integral of μ on the daily grid
cum_mu = np.concatenate([[0.0], np.cumsum(0.5*(mu_grid_hat[:-1]+mu_grid_hat[1:])*np.diff(grid))])
def int_mu(t_query):
    return np.interp(t_query, grid, cum_mu)

# Excitation compensator at event times (recursive, same trick as v4)
def Lam_events_v5(t, mu0, gamma, alpha, beta, S_events_, int_mu_func):
    Lam = np.zeros(len(t))
    S = 0.0
    for i, ti in enumerate(t):
        if i == 0:
            Lam[i] = int_mu_func(ti)
        else:
            dt = ti - t[i-1]
            S = np.exp(-beta*dt)*(S + 1.0)
            Lam[i] = int_mu_func(ti) + (alpha/beta)*(i - S)
    return Lam

Lam_v5 = Lam_events_v5(t, mu0_hat, gamma_hat, alpha_hat, beta_hat, S_events, int_mu)
tau_v5 = np.diff(Lam_v5)
ks_v5 = stats.kstest(tau_v5, "expon", args=(0, 1.0))
print(f"  v5 rescaled τ: mean={tau_v5.mean():.3f}, var={tau_v5.var():.3f}")
print(f"  v5 KS vs Exp(1): D={ks_v5.statistic:.3f}, p={ks_v5.pvalue:.3e}")

# v4 baseline KS for reference
mu_v4, a_v4, b_v4 = res_v4.x
def Lam_events_v4(t, mu, alpha, beta):
    Lam = np.zeros(len(t)); S = 0.0
    for i, ti in enumerate(t):
        if i == 0: Lam[i] = mu*ti
        else:
            dt = ti - t[i-1]
            S = np.exp(-beta*dt)*(S + 1.0)
            Lam[i] = mu*ti + (alpha/beta)*(i - S)
    return Lam
tau_v4 = np.diff(Lam_events_v4(t, mu_v4, a_v4, b_v4))
ks_v4 = stats.kstest(tau_v4, "expon", args=(0, 1.0))
print(f"  v4 KS vs Exp(1) (ref): D={ks_v4.statistic:.3f}, p={ks_v4.pvalue:.3e}")

# Ljung-Box on τ to test independence (excess clustering remaining?)
from scipy.stats import pearsonr
lag1_v5 = pearsonr(tau_v5[:-1], tau_v5[1:])
lag1_v4 = pearsonr(tau_v4[:-1], tau_v4[1:])
print(f"  Lag-1 autocorr of τ: v5 r={lag1_v5.statistic:+.3f} (p={lag1_v5.pvalue:.2e}),  "
      f"v4 r={lag1_v4.statistic:+.3f} (p={lag1_v4.pvalue:.2e})")

# ----------------------------------------------------------------------
# 6. Stochastic declustering under v5
# ----------------------------------------------------------------------
print("\n[declust] stochastic declustering (v5)…")
S_excite = 0.0
bg_prob_v5 = np.zeros(N)
mu_events_hat = mu0_hat * (S_events / S_bar)**gamma_hat
for i, ti in enumerate(t):
    if i == 0:
        rate_excite = 0.0
    else:
        dt = ti - t[i-1]
        S_excite = np.exp(-beta_hat*dt) * (S_excite + 1.0)
        rate_excite = alpha_hat * S_excite
    rate_total = mu_events_hat[i] + rate_excite
    bg_prob_v5[i] = mu_events_hat[i] / rate_total
expected_bg_v5 = bg_prob_v5.sum()
print(f"  Expected background events: {expected_bg_v5:.1f}  ({expected_bg_v5/N*100:.1f}%)")
print(f"  Expected offspring events:  {N-expected_bg_v5:.1f}  ({(N-expected_bg_v5)/N*100:.1f}%)")

# ----------------------------------------------------------------------
# 7. Forward simulation: 5,000 decades.  Need a forecast SSN profile for
#    "the next decade". We use two scenarios:
#      (A) "average decade": resample a historical SSN segment uniformly
#      (B) "current decade": use the actual SC25 → SC26 trajectory by
#          extending the most recent 10 years of SSN with its mean+amp
# ----------------------------------------------------------------------
print("\n[sim] Monte Carlo decadal simulations under v5…")

DECADE = 3652.5

def simulate_nshawkes(mu_fn, alpha, beta, T, rng, mu_max):
    """Ogata thinning with non-stationary μ(t). mu_max must upper-bound μ(t)
    for the simulation horizon."""
    events = []
    s = 0.0
    while s < T:
        if events:
            lam_star_excite = alpha*np.sum(np.exp(-beta*(s - np.array(events))))
        else:
            lam_star_excite = 0.0
        lam_upper = mu_max + lam_star_excite
        if lam_upper <= 0:
            break
        u = rng.random()
        w = -np.log(u) / lam_upper
        s += w
        if s >= T:
            break
        D = rng.random()
        if events:
            lam_excite = alpha*np.sum(np.exp(-beta*(s - np.array(events))))
        else:
            lam_excite = 0.0
        lam_s = mu_fn(s) + lam_excite
        if D * lam_upper <= lam_s:
            events.append(s)
    return np.array(events)

# Scenario A: random historical 10-yr SSN windows
all_starts = np.arange(0, T_end - DECADE, 30.0)
N_TRIALS = 5_000

def make_mu_fn_A(rng):
    start = rng.choice(all_starts)
    # μ(s) = μ0 * (S(start + s)/S_bar)^γ
    def fn(s):
        return mu0_hat * (S(start + s) / S_bar)**gamma_hat
    # upper bound: max SSN in segment
    seg = np.arange(start, start+DECADE+1, 1.0)
    mu_max_seg = mu0_hat * (np.nanmax(S(seg)) / S_bar)**gamma_hat
    return fn, mu_max_seg

# Scenario B: "current epoch" — use 2016-2025 (SC25 ascending → peak)
sc25_start = (pd.Timestamp("2016-01-01") - T0).days
def mu_fn_current(s):
    return mu0_hat * (S(sc25_start + s) / S_bar)**gamma_hat
mu_max_current = mu0_hat * (np.nanmax(S(np.arange(sc25_start, sc25_start+DECADE+1, 1.0))) / S_bar)**gamma_hat

counts_A = np.zeros(N_TRIALS, dtype=int)
counts_B = np.zeros(N_TRIALS, dtype=int)
counts_v4 = np.zeros(N_TRIALS, dtype=int)
counts_pois = rng.poisson(lam_pois*DECADE, size=N_TRIALS)
burst7_A = np.zeros(N_TRIALS, dtype=int)
burst7_B = np.zeros(N_TRIALS, dtype=int)

# v4 stationary Hawkes sim (reference)
def sim_stat_hawkes(mu, alpha, beta, T, rng):
    ev = []; s = 0.0
    while s < T:
        lam_star = mu + (alpha*np.sum(np.exp(-beta*(s - np.array(ev)))) if ev else 0.0)
        s += -np.log(rng.random())/lam_star
        if s >= T: break
        lam_s = mu + (alpha*np.sum(np.exp(-beta*(s - np.array(ev)))) if ev else 0.0)
        if rng.random()*lam_star <= lam_s:
            ev.append(s)
    return np.array(ev)

for i in range(N_TRIALS):
    # A: average historical decade
    fnA, muA_max = make_mu_fn_A(rng)
    simA = simulate_nshawkes(fnA, alpha_hat, beta_hat, DECADE, rng, muA_max)
    counts_A[i] = len(simA)
    if len(simA) >= 1:
        burst7_A[i] = max(((simA >= simA[j]) & (simA <= simA[j]+7)).sum() for j in range(len(simA)))
    # B: SC25-like
    simB = simulate_nshawkes(mu_fn_current, alpha_hat, beta_hat, DECADE, rng, mu_max_current)
    counts_B[i] = len(simB)
    if len(simB) >= 1:
        burst7_B[i] = max(((simB >= simB[j]) & (simB <= simB[j]+7)).sum() for j in range(len(simB)))
    # v4 ref
    counts_v4[i] = len(sim_stat_hawkes(mu_v4, a_v4, b_v4, DECADE, rng))

print(f"  Scenario A (random historical decade):")
print(f"    mean={counts_A.mean():.1f}  sd={counts_A.std():.1f}  "
      f"Var/Mean={counts_A.var()/counts_A.mean():.2f}  "
      f"95% CI=[{np.quantile(counts_A,0.025):.0f}, {np.quantile(counts_A,0.975):.0f}]")
print(f"    P(>=4 G4+ in any 7-day window) = {(burst7_A>=4).mean():.3f}")
print(f"  Scenario B (SC25-like ascending+peak decade):")
print(f"    mean={counts_B.mean():.1f}  sd={counts_B.std():.1f}  "
      f"Var/Mean={counts_B.var()/counts_B.mean():.2f}  "
      f"95% CI=[{np.quantile(counts_B,0.025):.0f}, {np.quantile(counts_B,0.975):.0f}]")
print(f"    P(>=4 G4+ in any 7-day window) = {(burst7_B>=4).mean():.3f}")
print(f"  v4 stationary Hawkes ref:")
print(f"    mean={counts_v4.mean():.1f}  sd={counts_v4.std():.1f}  "
      f"Var/Mean={counts_v4.var()/counts_v4.mean():.2f}")
print(f"  Poisson ref: mean={counts_pois.mean():.1f}  sd={counts_pois.std():.1f}")

# ----------------------------------------------------------------------
# 8. Figures
# ----------------------------------------------------------------------
# (12) μ(t) overlaid on events + SSN
fig, ax = plt.subplots(figsize=(12, 4.5))
years_grid = 1932 + grid/365.25
mu_year_grid = mu_grid_hat * 365.25
ax.fill_between(years_grid, 0, mu_year_grid, color="#5b8dba", alpha=0.45,
                label=f"μ(t) = μ0·(SSN/⟨SSN⟩)^γ,  μ0={mu0_hat*365.25:.2f}/yr, γ={gamma_hat:.2f}")
event_years = 1932 + t/365.25
ax.vlines(event_years, 0, mu_year_grid.max()*0.08, color="#b1361e", lw=0.6, alpha=0.7,
          label="G4+ events")
ax.set_xlabel("Year")
ax.set_ylabel("Background rate μ(t)  (events/year)")
ax.set_title(f"v5: SSN-modulated background, γ={gamma_hat:.2f}, η={eta_hat:.2f}")
ax.legend(loc="upper right")
ax.grid(alpha=0.3)
ax2 = ax.twinx()
ax2.plot(1932 + ssn_t/365.25, ssn_s, color="#888", lw=0.8, alpha=0.8)
ax2.set_ylabel("Smoothed SSN", color="#666")
plt.tight_layout()
plt.savefig(os.path.join(FIG, "12_v5_mu_of_t.png"), dpi=140)
plt.close()
print("\n[fig] 12_v5_mu_of_t.png")

# (13) QQ residuals v5 vs v4
fig, ax = plt.subplots(figsize=(6, 6))
sorted_v5 = np.sort(tau_v5); theor = stats.expon.ppf((np.arange(1,len(tau_v5)+1)-0.5)/len(tau_v5))
ax.plot(theor, sorted_v5, "o", ms=3, color="#2a7", alpha=0.7,
        label=f"v5 non-stat   KS p={ks_v5.pvalue:.2e}")
sorted_v4 = np.sort(tau_v4); theor4 = stats.expon.ppf((np.arange(1,len(tau_v4)+1)-0.5)/len(tau_v4))
ax.plot(theor4, sorted_v4, "x", ms=4, color="#5b8dba", alpha=0.6,
        label=f"v4 stationary   KS p={ks_v4.pvalue:.2e}")
xx = np.linspace(0, max(sorted_v5.max(), sorted_v4.max()), 2)
ax.plot(xx, xx, "--", color="#b1361e", lw=1.5, label="y = x")
ax.set_xlabel("Theoretical Exp(1) quantile")
ax.set_ylabel("Rescaled inter-arrival τ")
ax.set_title("Time-rescaling residuals: v5 vs v4")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "13_v5_residuals.png"), dpi=140)
plt.close()
print("[fig] 13_v5_residuals.png")

# (14) Decadal count distributions
fig, ax = plt.subplots(figsize=(10, 5))
bins = np.arange(0, 80, 2)
ax.hist(counts_pois, bins=bins, color="#bbb", alpha=0.55, edgecolor="white",
        label=f"Poisson μ={counts_pois.mean():.1f}")
ax.hist(counts_v4, bins=bins, color="#5b8dba", alpha=0.55, edgecolor="white",
        label=f"v4 stationary Hawkes μ={counts_v4.mean():.1f}")
ax.hist(counts_A, bins=bins, color="#2a7", alpha=0.55, edgecolor="white",
        label=f"v5 random historical decade μ={counts_A.mean():.1f}")
ax.hist(counts_B, bins=bins, color="#b1361e", alpha=0.55, edgecolor="white",
        label=f"v5 SC25-like decade μ={counts_B.mean():.1f}")
ax.set_xlabel("G4+ events per decade")
ax.set_ylabel("Monte Carlo trials")
ax.set_title("Decadal G4+ counts: Poisson vs v4 stationary vs v5 non-stationary Hawkes")
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "14_v5_decadal_counts.png"), dpi=140)
plt.close()
print("[fig] 14_v5_decadal_counts.png")

# ----------------------------------------------------------------------
# 9. Save summaries
# ----------------------------------------------------------------------
with open(os.path.join(DATA, "hawkes_v5_summary.txt"), "w") as f:
    f.write("v5: Non-stationary Hawkes (SSN-modulated background) — G4+ storms 1932-2025\n")
    f.write("="*70 + "\n\n")
    f.write(f"N events: {N}\nObservation window: {T_end/365.25:.2f} years\n")
    f.write(f"Mean smoothed SSN over window: {S_bar:.2f}\n\n")
    f.write("MLE parameters:\n")
    f.write(f"  μ0 = {mu0_hat:.5f} events/day  ({mu0_hat*365.25:.3f}/yr at S=S_bar)\n")
    f.write(f"  γ  = {gamma_hat:.4f}\n")
    f.write(f"  α  = {alpha_hat:.4f}\n")
    f.write(f"  β  = {beta_hat:.4f} 1/day   (1/β = {1/beta_hat:.2f} d)\n")
    f.write(f"  η  = {eta_hat:.4f}\n\n")
    f.write(f"Log-likelihood: v5={ll_v5:.2f}  v4={ll_v4:.2f}  Poisson={ll_pois:.2f}\n")
    f.write(f"AIC: v5={aic_v5:.2f}  v4={aic_v4:.2f}  Poisson={aic_pois:.2f}\n")
    f.write(f"ΔAIC (v5−v4) = {aic_v5-aic_v4:.2f}\n")
    f.write(f"LR (v5 vs v4) = {lr_v5_v4:.2f},  χ²(1) p = {p_v5_v4:.3e}\n\n")
    f.write("Time-rescaling residuals:\n")
    f.write(f"  v5: mean={tau_v5.mean():.3f} var={tau_v5.var():.3f} KS p={ks_v5.pvalue:.3e}\n")
    f.write(f"  v4: mean={tau_v4.mean():.3f} var={tau_v4.var():.3f} KS p={ks_v4.pvalue:.3e}\n")
    f.write(f"  Lag-1 autocorr v5: r={lag1_v5.statistic:+.3f} p={lag1_v5.pvalue:.2e}\n")
    f.write(f"  Lag-1 autocorr v4: r={lag1_v4.statistic:+.3f} p={lag1_v4.pvalue:.2e}\n\n")
    f.write("Stochastic declustering (v5):\n")
    f.write(f"  Background: {expected_bg_v5:.1f} ({expected_bg_v5/N*100:.1f}%)\n")
    f.write(f"  Offspring:  {N-expected_bg_v5:.1f} ({(N-expected_bg_v5)/N*100:.1f}%)\n\n")
    f.write(f"Monte Carlo ({N_TRIALS} decades, decade = {DECADE:.1f} days):\n")
    f.write(f"  v5 random historical: mean={counts_A.mean():.2f} sd={counts_A.std():.2f} "
            f"Var/Mean={counts_A.var()/counts_A.mean():.2f}\n")
    f.write(f"  v5 SC25-like decade:  mean={counts_B.mean():.2f} sd={counts_B.std():.2f} "
            f"Var/Mean={counts_B.var()/counts_B.mean():.2f}\n")
    f.write(f"  v4 stationary Hawkes: mean={counts_v4.mean():.2f} sd={counts_v4.std():.2f} "
            f"Var/Mean={counts_v4.var()/counts_v4.mean():.2f}\n")
    f.write(f"  Poisson reference:    mean={counts_pois.mean():.2f} sd={counts_pois.std():.2f}\n")
    f.write(f"  P(>=4 G4+ in any 7-day window per decade):\n")
    f.write(f"     v5 historical = {(burst7_A>=4).mean():.3f}\n")
    f.write(f"     v5 SC25-like  = {(burst7_B>=4).mean():.3f}\n")

pd.DataFrame({
    "event_index": np.arange(N),
    "t_days": t,
    "year": 1932 + t/365.25,
    "ssn_at_event": S_events,
    "mu_at_event_per_year": mu_events_hat*365.25,
    "p_background_v5": bg_prob_v5,
}).to_csv(os.path.join(DATA, "derived_hawkes_v5_declustering.csv"), index=False)

print("\n[done] v5 non-stationary Hawkes analysis complete")
