#!/usr/bin/env python3
"""
Formal Hawkes-process fit for G4+ geomagnetic storms, 1932-2025.

Model: univariate exponential-kernel Hawkes
    λ(t) = μ + Σ_{t_i < t}  α * exp(-β * (t - t_i))

Parameters:
    μ : background rate (events / day)
    α : excitation strength
    β : decay rate (1/days)
    η = α / β : branching ratio (mean number of direct offspring per event)

Stationarity requires η < 1.

References:
    Hawkes (1971) — original spectra paper
    Ogata (1981) — fast recursive log-likelihood for exp kernel
    Laub, Lee, Pollett, Taimre (2024) — review at arXiv:2405.10527
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
# 1. Load G4+ event times (days since 1932-01-01)
# ----------------------------------------------------------------------
daily = pd.read_csv(os.path.join(DATA, "derived_daily_with_phase.csv"),
                    parse_dates=["date"])
daily = daily.sort_values("date").reset_index(drop=True)
T0 = daily.date.min()
T_end = (daily.date.max() - T0).days  # observation window length, days

g4 = daily[daily.Kp_max >= 8].sort_values("date").reset_index(drop=True)
# Sub-day resolution: jitter ties within a day uniformly to avoid degeneracies
# (Kp is reported in 3-hour bins; daily aggregation creates ties at integer days)
rng = np.random.default_rng(20260523)
t_int = (g4.date - T0).dt.days.values.astype(float)
t = np.sort(t_int + rng.uniform(0.0, 1.0, size=len(t_int)))
N = len(t)
print(f"[load] N = {N} G4+ events, T = {T_end:,} days = {T_end/365.25:.1f} years")
print(f"       baseline rate (if Poisson) = {N/T_end*365.25:.2f} events/year")

# ----------------------------------------------------------------------
# 2. Ogata's recursive log-likelihood for exponential Hawkes
#
#  ln L = -μ T  +  Σ_i ln(μ + α R_i)  -  Σ_i (α/β)(1 - exp(-β(T - t_i)))
#
#  where R_1 = 0  and  R_i = exp(-β (t_i - t_{i-1})) * (1 + R_{i-1})
# ----------------------------------------------------------------------
def loglike(params, t, T):
    mu, alpha, beta = params
    if mu <= 0 or alpha < 0 or beta <= 0:
        return -1e18
    N = len(t)
    R = 0.0
    s_log = 0.0
    for i in range(N):
        if i > 0:
            dt = t[i] - t[i-1]
            if dt < 0: return -1e18
            R = np.exp(-beta*dt) * (1.0 + R)
        rate = mu + alpha * R
        if rate <= 0:
            return -1e18
        s_log += np.log(rate)
    # compensator
    s_comp = mu*T + (alpha/beta) * np.sum(1.0 - np.exp(-beta*(T - t)))
    return s_log - s_comp

def neg_ll(params, t, T):
    return -loglike(params, t, T)

# ----------------------------------------------------------------------
# 3. Multi-start MLE
# ----------------------------------------------------------------------
print("\n[fit] multi-start MLE for Hawkes (μ, α, β)…")
best = None
for trial, (mu0, eta0, beta0) in enumerate([
    (0.002, 0.30, 0.20),
    (0.003, 0.50, 0.50),
    (0.001, 0.70, 1.00),
    (0.004, 0.20, 0.10),
    (0.0015, 0.45, 0.30),
    (0.005, 0.60, 0.80),
]):
    x0 = [mu0, eta0*beta0, beta0]   # α = η * β
    res = optimize.minimize(neg_ll, x0, args=(t, T_end),
                            method="Nelder-Mead",
                            options={"xatol":1e-7, "fatol":1e-7, "maxiter":20000})
    if res.success or res.fun < 1e17:
        print(f"  trial {trial+1}: x0={x0}  →  μ={res.x[0]:.5f}  α={res.x[1]:.4f}  "
              f"β={res.x[2]:.4f}  η={res.x[1]/res.x[2]:.3f}  -LL={res.fun:.2f}")
        if best is None or res.fun < best.fun:
            best = res

mu_hat, alpha_hat, beta_hat = best.x
eta_hat = alpha_hat / beta_hat
ll_hawkes = -best.fun
aic_hawkes = 2*3 - 2*ll_hawkes
print(f"\n[MLE] best fit:")
print(f"  μ  = {mu_hat:.5f} events/day  =  {mu_hat*365.25:.3f} events/year (background)")
print(f"  α  = {alpha_hat:.4f}  (excitation amplitude)")
print(f"  β  = {beta_hat:.4f} 1/day,  i.e. decay timescale 1/β = {1/beta_hat:.2f} days")
print(f"  η  = α/β = {eta_hat:.3f}  (branching ratio: avg # direct offspring per event)")
print(f"  log-likelihood = {ll_hawkes:.2f},  AIC = {aic_hawkes:.2f}")

# Compare to Poisson baseline
lam_pois = N / T_end
ll_pois = N*np.log(lam_pois) - lam_pois*T_end
aic_pois = 2*1 - 2*ll_pois
print(f"\n  Homogeneous Poisson reference:")
print(f"    λ = {lam_pois:.5f} events/day,  log-likelihood = {ll_pois:.2f}, AIC = {aic_pois:.2f}")
print(f"    ΔAIC (Hawkes − Poisson) = {aic_hawkes - aic_pois:.1f}  (negative = Hawkes wins)")

# Likelihood ratio test (Hawkes nested under Poisson when α=0)
lr = 2*(ll_hawkes - ll_pois)
print(f"\n  Likelihood ratio statistic = {lr:.2f}  (df=2)")
print(f"  χ²(2) p-value = {1 - stats.chi2.cdf(lr, df=2):.3e}")

# ----------------------------------------------------------------------
# 4. Goodness-of-fit via time-rescaling theorem
#
#    Λ(t_i) = μ t_i + (α/β) Σ_{j<i} (1 − exp(−β(t_i − t_j)))
#    The transformed inter-arrivals τ_k = Λ(t_k) − Λ(t_{k-1})
#    should be i.i.d. Exp(1) if the model is correct.
# ----------------------------------------------------------------------
print("\n[GOF] time-rescaling residuals…")
def compensator_at_events(t, mu, alpha, beta):
    """Compute Λ at each event time."""
    Lam = np.zeros(len(t))
    # Σ_{j<i} (1 − exp(−β(t_i − t_j))) computed recursively
    # Let A_i = Σ_{j<i} exp(β t_j). Then Σ_{j<i} exp(−β(t_i − t_j)) = exp(−β t_i) * A_i
    # And the count of prior events is i.
    A = 0.0
    for i, ti in enumerate(t):
        Lam[i] = mu*ti + (alpha/beta) * (i - np.exp(-beta*ti)*A)
        A += np.exp(beta*ti)
    return Lam

# Numerical care for very large t with exp(beta*t): rescale via running formula
def compensator_at_events_stable(t, mu, alpha, beta):
    Lam = np.zeros(len(t))
    # Recursion: let S_i = Σ_{j<i} exp(-β (t_{i} − t_j)).  S_1 = 0.
    # Then S_{i+1} = exp(-β (t_{i+1} − t_i)) * (S_i + 1).
    # Σ_{j<i}(1 − exp(−β(t_i − t_j))) = i − S_i  (for i ≥ 1; treat first event as 0)
    S = 0.0
    for i, ti in enumerate(t):
        if i == 0:
            Lam[i] = mu*ti
        else:
            dt = ti - t[i-1]
            S = np.exp(-beta*dt)*(S + 1.0)
            Lam[i] = mu*ti + (alpha/beta)*(i - S)
    return Lam

Lam_events = compensator_at_events_stable(t, mu_hat, alpha_hat, beta_hat)
tau = np.diff(Lam_events)
print(f"  Number of rescaled inter-arrivals: {len(tau)}")
print(f"  Mean (should be 1):  {tau.mean():.3f}")
print(f"  Variance (should be 1): {tau.var():.3f}")

ks_resid = stats.kstest(tau, "expon", args=(0, 1.0))
print(f"  KS vs Exp(1):  D = {ks_resid.statistic:.3f},  p = {ks_resid.pvalue:.3e}")

# Same test for Poisson model
Lam_pois = lam_pois * t
tau_pois = np.diff(Lam_pois)
ks_pois = stats.kstest(tau_pois, "expon", args=(0, 1.0))
print(f"  Poisson reference KS: D = {ks_pois.statistic:.3f}, p = {ks_pois.pvalue:.3e}")

# ----------------------------------------------------------------------
# 5. Stochastic declustering — for each event, probability that it is
#    background (immigrant) vs. an offspring of an earlier event.
# ----------------------------------------------------------------------
print("\n[declust] stochastic declustering…")
# Probability event i is background = μ / λ*(t_i)
# Probability event i was triggered by event j (j<i) = α exp(-β(t_i - t_j)) / λ*(t_i)
S = 0.0
bg_prob = np.zeros(N)
parent = np.full(N, -1, dtype=int)  # most likely parent
for i, ti in enumerate(t):
    if i == 0:
        rate_excite = 0.0
    else:
        dt = ti - t[i-1]
        S = np.exp(-beta_hat*dt) * (S + 1.0)
        rate_excite = alpha_hat * S
    rate_total = mu_hat + rate_excite
    bg_prob[i] = mu_hat / rate_total
    if rate_excite > 0 and i > 0:
        # parent assignment: highest-contribution prior event
        contribs = alpha_hat * np.exp(-beta_hat*(ti - t[:i]))
        parent[i] = np.argmax(contribs)

expected_bg = bg_prob.sum()
expected_off = N - expected_bg
print(f"  Expected background (immigrant) events: {expected_bg:.1f}  ({expected_bg/N*100:.1f}%)")
print(f"  Expected offspring (excited) events:    {expected_off:.1f}  ({expected_off/N*100:.1f}%)")
print(f"  (Branching ratio η = {eta_hat:.3f} → asymptotic offspring fraction = {eta_hat:.3f})")

# ----------------------------------------------------------------------
# 6. Simulate Hawkes forward: 10,000 decade-long realizations
# ----------------------------------------------------------------------
print("\n[sim] Monte Carlo Hawkes simulations (Ogata thinning)…")
def simulate_hawkes(mu, alpha, beta, T, rng):
    """Ogata's thinning algorithm."""
    events = []
    s = 0.0
    while s < T:
        # current λ*(s)
        if events:
            lam_star = mu + alpha*np.sum(np.exp(-beta*(s - np.array(events))))
        else:
            lam_star = mu
        # candidate
        u = rng.random()
        w = -np.log(u) / lam_star
        s += w
        if s >= T:
            break
        D = rng.random()
        if events:
            lam_s = mu + alpha*np.sum(np.exp(-beta*(s - np.array(events))))
        else:
            lam_s = mu
        if D * lam_star <= lam_s:
            events.append(s)
    return np.array(events)

N_TRIALS = 5_000
DECADE = 3652.5  # days
counts_hawkes = np.zeros(N_TRIALS, dtype=int)
max_burst_hawkes = np.zeros(N_TRIALS, dtype=int)  # max events in any 7-day window
for i in range(N_TRIALS):
    sim = simulate_hawkes(mu_hat, alpha_hat, beta_hat, DECADE, rng)
    counts_hawkes[i] = len(sim)
    if len(sim) >= 2:
        # max events in any 7-day sliding window
        best = 1
        for j in range(len(sim)):
            cnt = ((sim >= sim[j]) & (sim <= sim[j]+7)).sum()
            if cnt > best: best = cnt
        max_burst_hawkes[i] = best
    else:
        max_burst_hawkes[i] = len(sim)

counts_pois = rng.poisson(lam_pois*DECADE, size=N_TRIALS)
print(f"  Decadal event counts:")
print(f"    Hawkes:  mean={counts_hawkes.mean():.1f}  sd={counts_hawkes.std():.1f}  "
      f"95% CI = [{np.quantile(counts_hawkes,0.025):.0f}, {np.quantile(counts_hawkes,0.975):.0f}]")
print(f"    Poisson: mean={counts_pois.mean():.1f}  sd={counts_pois.std():.1f}  "
      f"95% CI = [{np.quantile(counts_pois,0.025):.0f}, {np.quantile(counts_pois,0.975):.0f}]")
print(f"  Hawkes variance/mean = {counts_hawkes.var()/counts_hawkes.mean():.2f}  "
      f"(Poisson = 1.0; ratio shows over-dispersion)")

# Max 7-day burst
emp_burst = []
all_g4 = t.copy()
for j in range(len(all_g4)):
    cnt = ((all_g4 >= all_g4[j]) & (all_g4 <= all_g4[j]+7)).sum()
    emp_burst.append(cnt)
emp_max_burst_per_decade = max(emp_burst)  # over the full 94 yrs
print(f"\n  Max observed 7-day cluster (any time 1932-2025): {emp_max_burst_per_decade} G4+ days")
print(f"  Hawkes decadal max 7-day cluster: mean = {max_burst_hawkes.mean():.2f}, "
      f"P(>=3) = {(max_burst_hawkes>=3).mean():.3f},  P(>=4) = {(max_burst_hawkes>=4).mean():.3f}")

# ----------------------------------------------------------------------
# 7. Figures
# ----------------------------------------------------------------------
# (a) QQ plot of rescaled residuals
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax = axes[0]
sorted_tau = np.sort(tau)
theor = stats.expon.ppf((np.arange(1,len(tau)+1)-0.5)/len(tau))
ax.plot(theor, sorted_tau, "o", ms=3, color="#5b8dba", alpha=0.7, label="Hawkes residuals")
sorted_tau_p = np.sort(tau_pois)
theor_p = stats.expon.ppf((np.arange(1,len(tau_pois)+1)-0.5)/len(tau_pois))
ax.plot(theor_p, sorted_tau_p, "x", ms=4, color="#888", alpha=0.6, label="Poisson residuals")
xx = np.linspace(0, max(sorted_tau.max(), sorted_tau_p.max()), 2)
ax.plot(xx, xx, "--", color="#b1361e", lw=1.5, label="y = x  (perfect fit)")
ax.set_xlabel("Theoretical quantile (Exp(1))")
ax.set_ylabel("Rescaled inter-arrival τ")
ax.set_title(f"QQ-plot: time-rescaling residuals\nHawkes KS p={ks_resid.pvalue:.2e},  "
             f"Poisson KS p={ks_pois.pvalue:.2e}")
ax.legend()
ax.grid(alpha=0.3)

# (b) Decadal count comparison
ax = axes[1]
bins = np.arange(0, 50, 2)
ax.hist(counts_pois, bins=bins, color="#888", alpha=0.55, edgecolor="white",
        label=f"Poisson  μ={counts_pois.mean():.1f}, σ={counts_pois.std():.1f}")
ax.hist(counts_hawkes, bins=bins, color="#b1361e", alpha=0.65, edgecolor="white",
        label=f"Hawkes  μ={counts_hawkes.mean():.1f}, σ={counts_hawkes.std():.1f}")
ax.set_xlabel("G4+ events per decade")
ax.set_ylabel("Monte Carlo trials")
ax.set_title("Simulated decadal G4+ event counts: Hawkes vs Poisson")
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIG, "10_hawkes_gof_and_sim.png"), dpi=140)
plt.close()
print("\n[fig] 10_hawkes_gof_and_sim.png")

# (c) Conditional intensity λ*(t) over the 1932-2025 window with events overlaid
fig, ax = plt.subplots(figsize=(12, 4))
# Compute λ*(t) on a coarse grid in days, with running sum
tgrid = np.arange(0, T_end, 7.0)  # weekly
# For each grid point, sum kernel from prior events
lam_star = np.zeros_like(tgrid)
for k, tg in enumerate(tgrid):
    pri = t[t < tg]
    if len(pri):
        lam_star[k] = mu_hat + alpha_hat*np.sum(np.exp(-beta_hat*(tg - pri)))
    else:
        lam_star[k] = mu_hat
# Convert tgrid to years
years_grid = 1932 + tgrid/365.25
ax.fill_between(years_grid, 0, lam_star*365.25, color="#5b8dba", alpha=0.5,
                label="Conditional intensity λ*(t) (events/yr)")
event_years = 1932 + t/365.25
ax.vlines(event_years, 0, mu_hat*365.25*0.5, color="#b1361e", lw=0.5, alpha=0.8,
          label="G4+ events")
ax.axhline(mu_hat*365.25, color="black", ls="--", lw=1,
           label=f"Background μ = {mu_hat*365.25:.2f} /yr")
ax.set_xlabel("Year")
ax.set_ylabel("Intensity (events / year)")
ax.set_title(f"Hawkes conditional intensity, 1932–2025  "
             f"(μ={mu_hat*365.25:.2f}/yr, η={eta_hat:.2f}, 1/β={1/beta_hat:.1f}d)")
ax.legend(loc="upper right")
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "11_hawkes_intensity.png"), dpi=140)
plt.close()
print("[fig] 11_hawkes_intensity.png")

# ----------------------------------------------------------------------
# 8. Save summary
# ----------------------------------------------------------------------
with open(os.path.join(DATA, "hawkes_summary.txt"), "w") as f:
    f.write("Hawkes-process MLE fit for G4+ geomagnetic storms, 1932-2025\n")
    f.write("="*60 + "\n")
    f.write(f"N events: {N}\n")
    f.write(f"Observation window: {T_end/365.25:.2f} years\n\n")
    f.write(f"MLE parameters (exponential kernel):\n")
    f.write(f"  μ = {mu_hat:.5f} events/day  ({mu_hat*365.25:.3f}/yr) — background rate\n")
    f.write(f"  α = {alpha_hat:.4f} — excitation amplitude\n")
    f.write(f"  β = {beta_hat:.4f} 1/day  (decay timescale {1/beta_hat:.2f} d)\n")
    f.write(f"  η = α/β = {eta_hat:.3f} — branching ratio\n\n")
    f.write(f"Log-likelihood: {ll_hawkes:.2f}\n")
    f.write(f"AIC: Hawkes={aic_hawkes:.2f}   Poisson={aic_pois:.2f}   "
            f"Δ={aic_hawkes-aic_pois:.2f}\n")
    f.write(f"LR statistic: {lr:.2f}   χ²(2) p = {1 - stats.chi2.cdf(lr, df=2):.3e}\n\n")
    f.write(f"Time-rescaling residuals:\n")
    f.write(f"  mean={tau.mean():.3f} (expect 1)  var={tau.var():.3f} (expect 1)\n")
    f.write(f"  KS vs Exp(1): D={ks_resid.statistic:.3f}, p={ks_resid.pvalue:.3e}\n")
    f.write(f"  Poisson reference KS: D={ks_pois.statistic:.3f}, p={ks_pois.pvalue:.3e}\n\n")
    f.write(f"Stochastic declustering:\n")
    f.write(f"  Expected background events: {expected_bg:.1f}  ({expected_bg/N*100:.1f}%)\n")
    f.write(f"  Expected offspring events:  {expected_off:.1f}  ({expected_off/N*100:.1f}%)\n\n")
    f.write(f"Monte Carlo (N={N_TRIALS} decades):\n")
    f.write(f"  Hawkes mean count/decade  = {counts_hawkes.mean():.2f}, sd={counts_hawkes.std():.2f}\n")
    f.write(f"  Poisson mean count/decade = {counts_pois.mean():.2f}, sd={counts_pois.std():.2f}\n")
    f.write(f"  Hawkes Var/Mean = {counts_hawkes.var()/counts_hawkes.mean():.2f}  (Poisson=1.0)\n")
    f.write(f"  P(>=4 G4+ in any 7-day window per decade): {(max_burst_hawkes>=4).mean():.3f}\n")

# Save declustering output
pd.DataFrame({
    "event_index": np.arange(N),
    "t_days": t,
    "year": 1932 + t/365.25,
    "p_background": bg_prob,
    "most_likely_parent": parent,
}).to_csv(os.path.join(DATA, "derived_hawkes_declustering.csv"), index=False)

print("\n[done] Hawkes analysis complete")
