#!/usr/bin/env python3
"""
v6: Marked Hawkes process for G4+ geomagnetic storms, 1932-2025.

Extends v5 by giving each event a magnitude *mark* m_i (the Kp_max of that day,
which is 8.0, 8.33, 8.67, or 9.0). A larger event excites future activity more
strongly:

    λ(t, m) = μ(t) + Σ_{t_i < t}  α * exp(κ (m_i - m0)) * exp(-β (t - t_i))

    μ(t)   = μ0 * (S(t) / S_bar)^γ          (v5 background)
    α, β   : v5 kernel
    κ      : NEW — magnitude productivity exponent
    m0     : reference mark (we use 8.0, i.e. the G4 threshold)

So for a G4 day (m = 8.0)  the kernel amplitude is α.
For a G5 day (m = 9.0)     it is α * exp(κ).

κ > 0 ⇒ stronger storms produce more aftershocks (the seismology analogy
        with Båth's law / productivity scaling).
κ = 0 ⇒ recovers v5 (all events excite equally).
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
# 1. Load G4+ event times *with marks*
# ----------------------------------------------------------------------
daily = pd.read_csv(os.path.join(DATA, "derived_daily_with_phase.csv"),
                    parse_dates=["date"])
daily = daily.sort_values("date").reset_index(drop=True)
T0 = daily.date.min()
T_end = (daily.date.max() - T0).days

g4 = daily[daily.Kp_max >= 8].sort_values("date").reset_index(drop=True)
rng = np.random.default_rng(20260523)
t_int = (g4.date - T0).dt.days.values.astype(float)
jitter = rng.uniform(0.0, 1.0, size=len(t_int))
order = np.argsort(t_int + jitter)
t = (t_int + jitter)[order]
m = g4.Kp_max.values[order].astype(float)
N = len(t)
m0 = 8.0   # G4 threshold = reference mark
print(f"[load] N = {N} G4+ events  ({(m>=9).sum()} G5,  {(m<9).sum()} G4)")
print(f"       Kp_max distribution: 8.0={int((m==8.0).sum())}  "
      f"8.33={int((np.isclose(m,8.333)).sum())}  "
      f"8.67={int((np.isclose(m,8.667)).sum())}  9.0={int((m==9.0).sum())}")

# SSN (same as v5)
ssn = pd.read_csv(os.path.join(DATA, "SN_m_tot_V2.0.txt"),
                  sep=r"\s+", header=None,
                  names=["year","month","yfrac","sn","sd","n","prov"],
                  engine="python")
ssn = ssn[["yfrac","sn"]].copy()
ssn.loc[ssn.sn < 0, "sn"] = np.nan
ssn["sn_smooth"] = ssn.sn.rolling(window=13, center=True, min_periods=7).mean()
ssn = ssn.dropna().reset_index(drop=True)
T0_year = T0.year + (T0.month - 1)/12.0
ssn_t = (ssn.yfrac.values - T0_year) * 365.25
ssn_s = ssn.sn_smooth.values
mask = (ssn_t >= 0) & (ssn_t <= T_end)
S_bar = ssn_s[mask].mean()
print(f"[ssn] S_bar = {S_bar:.2f}")

def S(td):
    return np.interp(td, ssn_t, ssn_s, left=ssn_s[0], right=ssn_s[-1])

S_events = S(t)
grid = np.arange(0, T_end + 1.0, 1.0)
S_grid = S(grid)

# ----------------------------------------------------------------------
# 2. Marked Hawkes log-likelihood
#
#  ln L = Σ_i ln(μ(t_i) + α R_i)  -  ∫_0^T μ(t) dt
#                                   -  Σ_i (α/β) g_i (1 - exp(-β(T - t_i)))
#
#  where g_i = exp(κ (m_i - m0))  is the per-event productivity factor,
#  and the Ogata recursion is generalized to include the marks:
#
#     R_i  =  exp(-β (t_i - t_{i-1})) * (R_{i-1} + g_{i-1})
#     R_1  =  0
#
#  so that  α R_i = Σ_{j<i} α g_j exp(-β (t_i - t_j))  exactly.
# ----------------------------------------------------------------------
def loglike_v6(params, t, m, T, S_events, grid, S_grid, S_bar, m0):
    mu0, gamma, alpha, beta, kappa = params
    if mu0 <= 0 or alpha < 0 or beta <= 0:
        return -1e18
    mu_events = mu0 * (S_events / S_bar) ** gamma
    g = np.exp(kappa * (m - m0))
    R = 0.0
    s_log = 0.0
    for i in range(len(t)):
        if i > 0:
            dt = t[i] - t[i-1]
            R = np.exp(-beta * dt) * (R + g[i-1])
        rate = mu_events[i] + alpha * R
        if rate <= 0:
            return -1e18
        s_log += np.log(rate)
    mu_grid = mu0 * (S_grid / S_bar) ** gamma
    s_int_mu = np.trapezoid(mu_grid, grid)
    s_comp_exc = (alpha / beta) * np.sum(g * (1.0 - np.exp(-beta * (T - t))))
    return s_log - s_int_mu - s_comp_exc

def neg_ll_v6(p, *a):
    return -loglike_v6(p, *a)

# ----------------------------------------------------------------------
# 3. Multi-start MLE
# ----------------------------------------------------------------------
print("\n[fit] multi-start MLE for marked Hawkes (μ0, γ, α, β, κ)…")
starts = [
    (0.00548, 0.846, 0.1707, 0.6413, 0.0),   # v5 best + κ=0
    (0.00548, 0.846, 0.1707, 0.6413, 0.5),
    (0.00548, 0.846, 0.1707, 0.6413, 1.0),
    (0.005,   0.7,   0.15,   0.50,   0.3),
    (0.006,   1.0,   0.18,   0.70,   0.7),
    (0.004,   1.2,   0.12,   0.40,   1.5),
    (0.007,   0.5,   0.20,   0.80,   -0.5),  # negative κ as sanity check
    (0.005,   0.9,   0.16,   0.60,   2.0),
]
best = None
for trial, x0 in enumerate(starts):
    res = optimize.minimize(neg_ll_v6, x0,
                            args=(t, m, T_end, S_events, grid, S_grid, S_bar, m0),
                            method="Nelder-Mead",
                            options={"xatol":1e-7, "fatol":1e-7, "maxiter":60000})
    if res.fun < 1e17:
        mu0r, gr, ar, br, kr = res.x
        print(f"  trial {trial+1}: μ0={mu0r:.5f} γ={gr:.3f} α={ar:.4f} "
              f"β={br:.4f} κ={kr:+.3f}  -LL={res.fun:.2f}")
        if best is None or res.fun < best.fun:
            best = res

mu0_h, gamma_h, alpha_h, beta_h, kappa_h = best.x
eta_g4 = alpha_h / beta_h
eta_g5 = (alpha_h * np.exp(kappa_h)) / beta_h
ll_v6 = -best.fun
aic_v6 = 2*5 - 2*ll_v6
print(f"\n[MLE v6] best fit:")
print(f"  μ0 = {mu0_h:.5f} /day  ({mu0_h*365.25:.3f}/yr at S = S_bar)")
print(f"  γ  = {gamma_h:.3f}")
print(f"  α  = {alpha_h:.4f}")
print(f"  β  = {beta_h:.4f} /day  (1/β = {1/beta_h:.2f} d)")
print(f"  κ  = {kappa_h:+.4f}  → G5 productivity factor exp(κ) = {np.exp(kappa_h):.3f}×")
print(f"  branching ratio η for G4 parent = α/β       = {eta_g4:.3f}")
print(f"  branching ratio η for G5 parent = α e^κ / β = {eta_g5:.3f}")
print(f"  log-likelihood = {ll_v6:.2f},  AIC = {aic_v6:.2f}")

# v5 reference (refit so AIC comparison is on the same data jitter)
print("\n[compare] refitting v5 (κ ≡ 0)…")
def loglike_v5(p, *a): return loglike_v6(list(p)+[0.0], *a)
res_v5 = optimize.minimize(lambda p: -loglike_v5(p, t, m, T_end, S_events, grid, S_grid, S_bar, m0),
                           [0.00548, 0.846, 0.1707, 0.6413],
                           method="Nelder-Mead",
                           options={"xatol":1e-7,"fatol":1e-7,"maxiter":40000})
ll_v5 = -res_v5.fun
aic_v5 = 2*4 - 2*ll_v5
print(f"  v6: LL={ll_v6:.2f}  AIC={aic_v6:.2f}")
print(f"  v5: LL={ll_v5:.2f}  AIC={aic_v5:.2f}")
print(f"  ΔAIC (v6 − v5) = {aic_v6 - aic_v5:.2f}   (negative ⇒ v6 wins)")
lr = 2 * (ll_v6 - ll_v5)
p_lr = 1 - stats.chi2.cdf(lr, df=1)
print(f"  LR statistic = {lr:.2f},  χ²(1) p = {p_lr:.3e}")

# ----------------------------------------------------------------------
# 4. Time-rescaling GOF for v6
# ----------------------------------------------------------------------
print("\n[GOF] time-rescaling residuals…")
mu_grid_h = mu0_h * (S_grid / S_bar) ** gamma_h
cum_mu = np.concatenate([[0.0], np.cumsum(0.5*(mu_grid_h[:-1]+mu_grid_h[1:])*np.diff(grid))])
def int_mu(tq): return np.interp(tq, grid, cum_mu)

g_arr = np.exp(kappa_h * (m - m0))
Lam = np.zeros(N)
# Compensator at event times:  Λ(t_i) = ∫μ + (α/β) Σ_{j<i} g_j (1 - exp(-β(t_i - t_j)))
# Recursion: track A_i = Σ_{j<i} g_j exp(-β(t_i - t_j)) and G_i = Σ_{j<i} g_j
A = 0.0
G = 0.0
for i, ti in enumerate(t):
    if i == 0:
        Lam[i] = int_mu(ti)
    else:
        dt = ti - t[i-1]
        A = np.exp(-beta_h * dt) * (A + g_arr[i-1])
        G = G + g_arr[i-1]
        Lam[i] = int_mu(ti) + (alpha_h / beta_h) * (G - A)
tau_v6 = np.diff(Lam)
ks_v6 = stats.kstest(tau_v6, "expon", args=(0, 1.0))
print(f"  v6 τ: mean={tau_v6.mean():.3f}, var={tau_v6.var():.3f}")
print(f"  v6 KS vs Exp(1): D={ks_v6.statistic:.3f}, p={ks_v6.pvalue:.3e}")
from scipy.stats import pearsonr
lag1 = pearsonr(tau_v6[:-1], tau_v6[1:])
print(f"  v6 lag-1 autocorr: r={lag1.statistic:+.3f}, p={lag1.pvalue:.2e}")

# ----------------------------------------------------------------------
# 5. Empirical sanity check: do G5s actually have more close-followers?
#    Look at the count of G4+ events in the 7 days *following* each event,
#    split by parent magnitude.
# ----------------------------------------------------------------------
print("\n[empirical] aftershock window (next 7 days), by parent mark:")
for m_val in [8.0, 8.333, 8.667, 9.0]:
    idx = np.where(np.isclose(m, m_val))[0]
    follow = []
    for i in idx:
        ti = t[i]
        cnt = ((t > ti) & (t <= ti + 7)).sum()
        follow.append(cnt)
    print(f"  parent Kp={m_val:.2f}  (n={len(idx):3d}):  "
          f"mean followers in next 7 d = {np.mean(follow):.3f}  "
          f"sd = {np.std(follow):.3f}  max = {max(follow) if follow else 0}")

# ----------------------------------------------------------------------
# 6. Simulate marked Hawkes forward — what does a G5 trigger?
#    For each scenario, plant a single G5 at t=0 and simulate offspring
#    via the immigration-and-birth process.
# ----------------------------------------------------------------------
print("\n[sim] direct-offspring count by parent magnitude:")
N_SIM = 50_000
def sim_offspring(parent_mark, kappa, alpha, beta, rng, horizon=120.0):
    """Generations: G_0 starts with one event at t=0 of given mark.
    For each event, number of direct offspring is Poisson(α e^{κ(m_p-m0)} / β),
    and their times are exponential(β) after the parent.
    Children inherit the *empirical* mark distribution (not modeled).
    Returns the total descendant count over [0, horizon] (excluding root)."""
    # we only track direct offspring for productivity sanity check
    mu_off = (alpha * np.exp(kappa * (parent_mark - m0))) / beta
    nchild = rng.poisson(mu_off)
    # ages of children
    if nchild == 0:
        return 0
    ages = rng.exponential(1/beta, size=nchild)
    return int((ages <= horizon).sum())

for m_val in [8.0, 8.333, 8.667, 9.0]:
    offs = [sim_offspring(m_val, kappa_h, alpha_h, beta_h, rng) for _ in range(N_SIM)]
    print(f"  Kp parent = {m_val:.2f}:  E[direct offspring] = {np.mean(offs):.3f}  "
          f"P(≥1 follower) = {np.mean(np.array(offs)>=1):.3f}  "
          f"P(≥3 followers) = {np.mean(np.array(offs)>=3):.3f}")

# Full forward decadal sim (with mark resampling)
print("\n[sim] full decadal Monte Carlo (v6, with empirical mark distribution)…")
mark_pool = m.copy()  # bootstrap pool for new event marks
DECADE = 3652.5
N_TRIALS = 3000

def simulate_v6_full(mu_fn, alpha, beta, kappa, T, rng, mu_max, mark_pool):
    events_t, events_m = [], []
    s = 0.0
    while s < T:
        if events_t:
            ev_t = np.array(events_t); ev_m = np.array(events_m)
            g_e = np.exp(kappa * (ev_m - m0))
            excite = alpha * np.sum(g_e * np.exp(-beta * (s - ev_t)))
        else:
            excite = 0.0
        lam_up = mu_max + excite
        if lam_up <= 0: break
        s += -np.log(rng.random()) / lam_up
        if s >= T: break
        if events_t:
            ev_t = np.array(events_t); ev_m = np.array(events_m)
            g_e = np.exp(kappa * (ev_m - m0))
            excite_s = alpha * np.sum(g_e * np.exp(-beta * (s - ev_t)))
        else:
            excite_s = 0.0
        lam_s = mu_fn(s) + excite_s
        if rng.random() * lam_up <= lam_s:
            events_t.append(s)
            events_m.append(rng.choice(mark_pool))
    return np.array(events_t), np.array(events_m)

# Two scenarios same as v5: historical random decade vs SC25-like
all_starts = np.arange(0, T_end - DECADE, 30.0)
sc25_start = (pd.Timestamp("2016-01-01") - T0).days

counts_A = np.zeros(N_TRIALS, dtype=int)
counts_B = np.zeros(N_TRIALS, dtype=int)
g5_A = np.zeros(N_TRIALS, dtype=int)
g5_B = np.zeros(N_TRIALS, dtype=int)
burst_A = np.zeros(N_TRIALS, dtype=int)
burst_B = np.zeros(N_TRIALS, dtype=int)

for i in range(N_TRIALS):
    # A: random historical decade
    start = rng.choice(all_starts)
    seg = np.arange(start, start+DECADE+1, 1.0)
    mu_max_A = mu0_h * (np.nanmax(S(seg)) / S_bar) ** gamma_h
    def mu_fn_A(s, start=start):
        return mu0_h * (S(start + s) / S_bar) ** gamma_h
    eA_t, eA_m = simulate_v6_full(mu_fn_A, alpha_h, beta_h, kappa_h, DECADE, rng, mu_max_A, mark_pool)
    counts_A[i] = len(eA_t)
    g5_A[i] = int((eA_m >= 9.0).sum())
    if len(eA_t):
        burst_A[i] = max(((eA_t >= eA_t[j]) & (eA_t <= eA_t[j]+7)).sum() for j in range(len(eA_t)))
    # B: SC25-like
    seg = np.arange(sc25_start, sc25_start+DECADE+1, 1.0)
    mu_max_B = mu0_h * (np.nanmax(S(seg)) / S_bar) ** gamma_h
    def mu_fn_B(s): return mu0_h * (S(sc25_start + s) / S_bar) ** gamma_h
    eB_t, eB_m = simulate_v6_full(mu_fn_B, alpha_h, beta_h, kappa_h, DECADE, rng, mu_max_B, mark_pool)
    counts_B[i] = len(eB_t)
    g5_B[i] = int((eB_m >= 9.0).sum())
    if len(eB_t):
        burst_B[i] = max(((eB_t >= eB_t[j]) & (eB_t <= eB_t[j]+7)).sum() for j in range(len(eB_t)))

print(f"  Scenario A (random historical decade):")
print(f"    G4+ count mean={counts_A.mean():.1f}  G5 mean={g5_A.mean():.2f}  "
      f"P(>=1 G5/decade)={(g5_A>=1).mean():.3f}  P(>=2 G5/decade)={(g5_A>=2).mean():.3f}")
print(f"    P(>=4 in 7d) = {(burst_A>=4).mean():.3f}")
print(f"  Scenario B (SC25-like decade):")
print(f"    G4+ count mean={counts_B.mean():.1f}  G5 mean={g5_B.mean():.2f}  "
      f"P(>=1 G5/decade)={(g5_B>=1).mean():.3f}  P(>=2 G5/decade)={(g5_B>=2).mean():.3f}")
print(f"    P(>=4 in 7d) = {(burst_B>=4).mean():.3f}")

# ----------------------------------------------------------------------
# 7. Figures
# ----------------------------------------------------------------------
# (15) productivity scaling with mark
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax = axes[0]
m_grid = np.linspace(8.0, 9.0, 100)
prod = np.exp(kappa_h * (m_grid - m0))
ax.plot(m_grid, prod, color="#b1361e", lw=2)
ax.axhline(1.0, color="gray", lw=0.8, ls="--")
for m_val in [8.0, 8.333, 8.667, 9.0]:
    ax.axvline(m_val, color="#aaa", lw=0.5, ls=":")
ax.fill_between(m_grid, prod, color="#b1361e", alpha=0.1)
ax.set_xlabel("Parent Kp magnitude")
ax.set_ylabel("Excitation multiplier exp(κ·(m−8))")
ax.set_title(f"Marked Hawkes productivity: G5 ≈ {np.exp(kappa_h):.2f}× the kick of a G4")
ax.grid(alpha=0.3)

# (b) Empirical mean followers in next 7d by mark — bar plot
emp_mean = []
labels = []
for m_val in [8.0, 8.333, 8.667, 9.0]:
    idx = np.where(np.isclose(m, m_val))[0]
    follow = []
    for i in idx:
        ti = t[i]
        cnt = ((t > ti) & (t <= ti + 7)).sum()
        follow.append(cnt)
    emp_mean.append(np.mean(follow) if follow else 0.0)
    labels.append(f"Kp={m_val:.2f}\nn={len(idx)}")
ax = axes[1]
xs = np.arange(len(labels))
bars = ax.bar(xs, emp_mean, color=["#5b8dba","#5b8dba","#b1361e","#b1361e"], edgecolor="white")
for i, v in enumerate(emp_mean):
    ax.text(i, v+0.01, f"{v:.2f}", ha="center", fontsize=10)
ax.set_xticks(xs); ax.set_xticklabels(labels)
ax.set_ylabel("Mean G4+ events in next 7 days")
ax.set_title("Empirical: G5 days are followed by more close events than G4 days")
ax.grid(alpha=0.3, axis="y")
plt.tight_layout()
plt.savefig(os.path.join(FIG, "15_v6_mark_productivity.png"), dpi=140)
plt.close()
print("\n[fig] 15_v6_mark_productivity.png")

# (16) G5 count distribution per decade
fig, ax = plt.subplots(figsize=(9, 5))
bins = np.arange(0, 12) - 0.5
ax.hist(g5_A, bins=bins, color="#2a7", alpha=0.55, edgecolor="white",
        label=f"v6 random historical decade  mean={g5_A.mean():.2f}")
ax.hist(g5_B, bins=bins, color="#b1361e", alpha=0.55, edgecolor="white",
        label=f"v6 SC25-like decade  mean={g5_B.mean():.2f}")
ax.set_xlabel("G5 (Kp=9) days per decade")
ax.set_ylabel("Monte Carlo trials")
ax.set_title("How many G5 days per decade? v6 forecast")
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "16_v6_g5_decadal.png"), dpi=140)
plt.close()
print("[fig] 16_v6_g5_decadal.png")

# ----------------------------------------------------------------------
# 8. Save summary
# ----------------------------------------------------------------------
with open(os.path.join(DATA, "hawkes_v6_summary.txt"), "w") as f:
    f.write("v6: Marked Hawkes (mark = Kp_max) — G4+ storms 1932-2025\n")
    f.write("="*60 + "\n\n")
    f.write(f"N = {N}  ({(m>=9).sum()} G5, {(m<9).sum()} G4)\n\n")
    f.write("MLE:\n")
    f.write(f"  μ0   = {mu0_h:.5f} /day  ({mu0_h*365.25:.3f}/yr)\n")
    f.write(f"  γ    = {gamma_h:.4f}\n")
    f.write(f"  α    = {alpha_h:.4f}\n")
    f.write(f"  β    = {beta_h:.4f}   (1/β = {1/beta_h:.2f} d)\n")
    f.write(f"  κ    = {kappa_h:+.4f}   ⇒ G5 productivity = {np.exp(kappa_h):.3f}× a G4\n\n")
    f.write(f"  η(G4) = α/β       = {eta_g4:.3f}\n")
    f.write(f"  η(G5) = α e^κ / β = {eta_g5:.3f}\n\n")
    f.write(f"AIC: v6={aic_v6:.2f}  v5={aic_v5:.2f}  Δ={aic_v6-aic_v5:.2f}\n")
    f.write(f"LR (v6 vs v5) = {lr:.2f}   χ²(1) p = {p_lr:.3e}\n\n")
    f.write(f"Time-rescaling τ: mean={tau_v6.mean():.3f} var={tau_v6.var():.3f}\n")
    f.write(f"  KS vs Exp(1): D={ks_v6.statistic:.3f}  p={ks_v6.pvalue:.3e}\n")
    f.write(f"  lag-1 autocorr: r={lag1.statistic:+.3f}  p={lag1.pvalue:.2e}\n\n")
    f.write("Monte Carlo (3,000 decades each):\n")
    f.write(f"  Scenario A (random hist):  G4+/dec mean={counts_A.mean():.2f}  "
            f"G5/dec mean={g5_A.mean():.2f}  P(>=1 G5)={(g5_A>=1).mean():.3f}  "
            f"P(>=4 in 7d)={(burst_A>=4).mean():.3f}\n")
    f.write(f"  Scenario B (SC25-like):    G4+/dec mean={counts_B.mean():.2f}  "
            f"G5/dec mean={g5_B.mean():.2f}  P(>=1 G5)={(g5_B>=1).mean():.3f}  "
            f"P(>=4 in 7d)={(burst_B>=4).mean():.3f}\n")

print("\n[done] v6 marked-Hawkes analysis complete")
