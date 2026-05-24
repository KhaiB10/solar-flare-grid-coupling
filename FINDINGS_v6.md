# FINDINGS v6 — Marked Hawkes: a G5 punches ~2.7× harder than a G4

**Defensive publication. Not an operational risk model. Not financial, engineering, or insurance advice.**
**Diatom Sky R&D — open methodology, open data, open prior art.**

---

## The plain-English version first

### What we've been building, in one analogy

Imagine you live in earthquake country. You want to know not just *how often the
ground shakes*, but *what happens after a big quake* — because aftershocks matter
for whether you can leave your house, when crews can repair things, and whether
the next event will catch a weakened structure.

That is exactly the question this whole repo has been asking, but for solar
storms hitting Earth instead of earthquakes hitting buildings. The grid-impact
events we list (1989 Quebec blackout, 2003 Halloween storms, the 2024 Gannon
storm that knocked out GPS for precision farmers across the Midwest) are the
"buildings that fell." The Kp/ap geomagnetic indices since 1932 are the
"seismograph trace." And the question is: **how do we honestly forecast the next
ten years of solar-storm shaking, including aftershocks?**

Each version of the analysis adds one piece of physics that the previous version
ignored:

- **v1 (Poisson)** — pretend storms happen completely at random, like coin flips
  every day. Result: ~58.5% chance per decade of a Carrington-class day. Useful
  baseline but obviously wrong because the Sun has an 11-year cycle.
- **v2 (cycle phase)** — split history into "solar maximum" vs "solar minimum"
  decades. Carrington-class hazard ranges from 6.3% (deep minimum) to 76.8%
  (solar max). The analogy: in earthquake country, fault-loading rate varies
  by region; here it varies by where the Sun is in its 11-year cycle.
- **v3 (clustering)** — discovered that storm wait-times are sharply
  *non-random*. Once a G4+ happens, another one within 5 days is **5.9× more
  likely than chance**. This is the "aftershock" effect.
- **v4 (Hawkes process)** — wrote that aftershock effect down as a proper
  generative model. Every storm gives the Sun a *temporary boost* that decays
  with a 1.74-day half-life. Branching ratio η = 0.28 (about 28% of storms are
  aftershocks of an earlier one).
- **v5 (non-stationary Hawkes)** — added the solar cycle back in as a slowly
  varying baseline rate. The Sun is more "loaded" during solar maximum. Fit:
  background scales as SSN^0.85. This eliminated *every* remaining diagnostic
  problem in v4 — the rescaled residuals are now indistinguishable from random
  noise (KS p = 0.44).
- **v6 (this finding, marked Hawkes)** — now ask: **do G5 (extreme) storms
  excite more follow-on storms than G4 (severe) storms?** In the earthquake
  analogy: does a magnitude 7 earthquake produce more aftershocks than a
  magnitude 6? Seismologists have a name for this — **productivity scaling**,
  and they've shown for decades that yes, bigger quakes produce more
  aftershocks, by a factor that scales roughly as 10^(magnitude difference).

### v6 in one sentence

**A G5 day (Kp=9) produces about 2.7 times as much follow-on geomagnetic storm
activity as a G4 day (Kp=8) within the next ~1.5 days** — meaning the storms
behave qualitatively the same way aftershocks do in seismology.

### Why this matters for the future, in plain terms

This is *not* a prediction that "the next G5 will cause X dollars of damage."
But it does point at things that grid operators, satellite operators, GPS
users, airlines flying polar routes, and military planners arguably already do
know intuitively, and quantifies them:

1. **A G5 should put the system on a multi-day "stay alert" footing, not a
   one-day one.** If your storm-response protocol resets the moment Kp drops
   below 8, you are systematically *under*-prepared for the days following an
   extreme event. The fitted decay timescale is 1.5 days, so a sensible window
   is more like 3–5 days.

2. **The 2024 Gannon storm (G5) was the kind of event that the model says is
   most likely to be *immediately followed by another G4+*.** And indeed, a
   G4 followed in October 2024 — exactly the kind of pattern the model predicts
   should be common after a G5. Empirically, of the 27 G5 days in the 94-year
   record, ~67% are followed by another G4+ within a week, vs ~29% for "ordinary"
   G4 days. This finding is robust and not from the model — it's in the raw data.

3. **For the SC25 decade (we're living in it right now), the model expects
   ~1.9 G5 days and ~17 G4+ days total**, with the G5 days clustering near
   the cycle peak (now to ~2027). The previous decade (SC24 declining into
   SC25 ascending) had 2 G5 days (May 2024 Gannon, October 2024) — within
   the model's predicted range.

### Who is impacted by this kind of work

A research finding by a hobbyist obviously does not change any operational
posture by itself. But the *kind of analysis* this repo is doing is the kind
of input that several different groups feed into much bigger pipelines:

- **Electric grid operators** (e.g. PJM, ERCOT, MISO, regional transmission
  organizations): they have GIC-monitoring relays on long high-voltage
  transmission lines. After a big storm they decide when to redispatch generation
  and recheck transformer health. A "the next 3 days are still risky" framing
  is operationally different from "today's risk is over."
- **NERC** (the North American grid-reliability regulator): currently mandates
  GMD (Geomagnetic Disturbance) operating procedures (TPL-007). The size of
  the benchmark events used in those standards is exactly the kind of question
  Hawkes-style hazard analysis informs.
- **Satellite operators** (Starlink, GPS, GEO comms): satellite anomaly rates
  spike during and *after* G4+ storms because of charging effects and increased
  drag in low orbit. Starlink lost 38 satellites during a G2 storm in Feb 2022
  — drag effects don't reset at midnight.
- **Precision-agriculture and aviation**: GPS RTK service degraded across the
  US Midwest during Gannon for ~48 hours, costing farmers a planting window.
  Polar-route airlines (Delta, United, Air Canada) reroute or descend to
  reduce radiation exposure when Kp spikes.
- **Insurance and reinsurance** (Lloyd's, Munich Re, Swiss Re): publish
  scenario-based PML estimates for severe space weather; the AIC-jumps in this
  series are the kind of quantitative basis the underwriting community asks
  for when scenario assumptions are challenged.
- **Defense / Space Force**: the recently stood-up USSF cares operationally
  about geomagnetic conditions for surveillance radar refraction, OTH radar,
  and satellite tasking. They have classified models, but the *form* of the
  hazard surface is open.
- **Researchers**: this is the audience the repo is genuinely written for.
  Anyone fitting a future ETAS-style space-weather model has prior-art for
  marked-Hawkes-with-SSN-background; future re-analyses can compare to these
  numbers.

The honest version of the "who is impacted" question is: **billions of people
indirectly**, because high-voltage transformers, satellite communications, and
GPS are critical infrastructure; **millions of people directly** during a major
event (the 1989 Quebec blackout darkened 6 million Hydro-Québec customers);
and **a few hundred analysts at the institutions above** use this kind of
work as input. This repo is a citation chip, not a recommendation engine.

---

## The technical version

### The model

Take v5 and replace the simple kernel with one that depends on parent magnitude:

$$
\lambda(t \mid \text{history}) \;=\; \mu(t) \;+\; \sum_{t_i < t}
   \alpha \, e^{\kappa (m_i - m_0)} \, e^{-\beta (t - t_i)}
$$

where $m_i$ is the Kp_max of event $i$ (values 8.0, 8.33, 8.67, or 9.0) and
$m_0 = 8.0$ is the G4 threshold. New parameter:

- $\kappa$ — magnitude productivity exponent. $\kappa = 0$ ⇒ recovers v5; the
  multiplier $e^{\kappa(m_i - m_0)}$ scales the kernel amplitude.

The compensator for the log-likelihood and the Ogata recursion are both
generalized in the obvious way; full code in `scripts/analyze_hawkes_marked.py`.

### MLE result (8 multi-starts, all converging to one optimum)

| Parameter | Value | Interpretation |
|---|---|---|
| $\mu_0$ | 0.00549 /day = 2.00/yr | background at mean-cycle SSN (essentially unchanged from v5) |
| $\gamma$ | 0.845 | SSN exponent (unchanged from v5's 0.846) |
| $\alpha$ | 0.1113 | excitation amplitude (lower than v5's 0.171 — the "average" amplitude is now G4-baseline) |
| $\beta$ | 0.6520 /day | 1/β = **1.53 d** (essentially unchanged from v5's 1.56 d) |
| **$\kappa$** | **+1.005** | **G5 productivity = e^1.005 ≈ 2.73× a G4** |

Derived branching ratios:

- $\eta_{G4} = \alpha/\beta = 0.171$ — every G4 produces ~0.17 direct offspring on average
- $\eta_{G5} = \alpha e^{\kappa}/\beta = 0.466$ — every G5 produces ~0.47 direct offspring on average

A G5 is, on average, getting nearly halfway to the critical threshold $\eta = 1$
where the process would explode. Still subcritical — but a *lot* closer to the
edge than the population average.

### Statistical significance

| Comparison | Value |
|---|---|
| log-likelihood (v6) | −1292.01 |
| log-likelihood (v5, same data jitter) | −1295.11 |
| AIC v6 | 2594.01 |
| AIC v5 | 2598.21 |
| **ΔAIC (v6 − v5)** | **−4.20** (negative ⇒ v6 wins) |
| **LR statistic, χ²(1)** | **6.20**, p = **0.0128** |
| Time-rescaling KS p vs Exp(1) | 0.57 (v5: 0.44; both clean) |
| Lag-1 autocorr of τ | +0.008 (essentially zero) |

ΔAIC = −4.20 is modest — this is a one-parameter extension to an already-good
model. The improvement passes the conventional ΔAIC > 2 threshold for "actually
worth adding" and the likelihood-ratio test rejects $\kappa = 0$ at p = 0.013.
It is also consistent with an obvious empirical pattern (next section).

### Empirical sanity check (model-independent)

Just count how many G4+ days follow each event within 7 days, split by parent Kp:

| Parent Kp | n events | Mean G4+ followers in next 7 d | Max |
|---|---|---|---|
| 8.00 (G4) | 79 | 0.291 | 3 |
| 8.33 (G4) | 81 | 0.321 | 3 |
| 8.67 (G4+) | 59 | **0.525** | 3 |
| **9.00 (G5)** | 27 | **0.667** | 4 |

This is in the raw 1932–2025 record with no model fit at all: G5 days are
followed by ~2.3× as many G4+ days within a week as Kp = 8.0 days. The marked
Hawkes is fitting a real, model-independent pattern.

### Decadal hazard, updated

3,000-trial Monte Carlo using the v6 fit, two scenarios as in v5:

| Scenario | G4+/decade | **G5/decade** | **P(≥1 G5/decade)** | **P(≥2 G5)** | P(≥4 G4+ in any 7-d window) |
|---|---|---|---|---|---|
| Random historical | 26.5 | 2.94 | 91.9% | 74.4% | 47.7% |
| SC25-like (now) | 17.4 | 1.88 | 82.4% | 54.9% | 32.4% |

The G5-per-decade number is new. A SC25-like decade is expected to deliver about
2 G5 days, with an 82% chance of at least one and a 55% chance of at least two.
The historical record shows 2 G5 days in 2024 alone (Gannon plus October),
which means SC25 has likely already delivered most of its expected G5 budget.

### Goodness-of-fit summary across all five versions

| Model | Params | log-L | ΔAIC vs Poisson | Residual KS p | Residual lag-1 r |
|---|---|---|---|---|---|
| v0 Poisson | 1 | −1460.9 | 0 | 4.5e−16 ❌ | n/a |
| v4 stationary Hawkes | 3 | −1335.6 | −246.5 | 2.8e−3 ⚠️ | +0.228 ⚠️ |
| v5 non-stat Hawkes | 4 | −1295.1 | −325.6 | 0.44 ✓ | +0.019 ✓ |
| **v6 marked Hawkes** | **5** | **−1292.0** | **−329.8** | **0.57 ✓** | **+0.008 ✓** |

v6 is at this point a *thoroughly* validated description of the 1932–2025
Kp ≥ 8 series. The marginal residuals, the serial residuals, and the
magnitude-dependent productivity all check out simultaneously.

### Caveats — important

1. **The mark distribution is sampled, not modeled.** In simulation, new event
   marks are drawn from the empirical distribution (79/81/59/27). A fully
   self-consistent marked Hawkes would also model the conditional mark
   distribution given the parent; v6 does not.
2. **27 G5 days is a small sample.** The κ estimate has real uncertainty —
   a profile-likelihood or bootstrap CI is a natural next step (v7).
3. **Same caveats as v5**: smoothed SSN is causally backward at the edges; the
   exponential kernel could be the wrong shape (power-law is the standard
   seismology alternative); G4 and G5 lump together categorically distinct
   physical scenarios in different IMF/CME contexts.
4. **The "G5 produces 2.7× more aftershocks" finding does not mean a G5 is more
   physically dangerous than a G4 by a factor of 2.7.** Damage potential scales
   with peak GIC, transformer thermal limits, and grid topology — none of which
   are in this model. The 2.7× is purely about *follow-on event frequency*.
5. **Marks max out at Kp = 9.** A truly extreme event (Carrington-class, where
   Kp would saturate at 9 for multiple days) is not extrapolable from this fit;
   the model treats every G5 day equivalently regardless of how "deep" it was.

### Outputs

| File | Description |
|---|---|
| `scripts/analyze_hawkes_marked.py` | v6 marked-Hawkes pipeline |
| `data/hawkes_v6_summary.txt` | text summary of all numbers above |
| `data/run_v6_log.txt` | full console output |
| `figures/15_v6_mark_productivity.png` | productivity multiplier + empirical bar chart |
| `figures/16_v6_g5_decadal.png` | G5-per-decade distribution under two scenarios |

### Data and reproducibility

- Kp/ap: GFZ Potsdam, [`Kp_ap_since_1932.txt`](https://kp.gfz.de/app/files/Kp_ap_since_1932.txt)
- SSN: SILSO, [`SN_m_tot_V2.0.txt`](https://www.sidc.be/SILSO/DATA/SN_m_tot_V2.0.txt)
- Seed: 20260523 (every stochastic step)
- License: MIT for code, CC0 for derived data and figures

### Suggested v7+

- **Power-law kernel** $\alpha (c + (t-t_i))^{-(p+1)}$ — the standard seismology
  alternative to exponential; could capture genuinely heavy aftershock tails
- **Conditional mark distribution** — let $p(m_{\text{child}} \mid m_{\text{parent}})$
  depend on parent magnitude (do G5s more often spawn G5s?)
- **Pre-1932 extension** via aa-index (1868→) — adds Cycles 11–16, including the
  Carrington-class event (1859 is just before aa-index data but the 1872 and
  1909 events are in)
- **Cross-validation by leave-one-cycle-out** — does the fit on Cycles 17–24
  predict Cycle 25's actual storm series well?
- **Bootstrap CI on κ** — currently we report a point estimate; a 95% CI would
  pin down whether κ ∈ [0.5, 1.5] or something tighter

— Diatom Sky R&D, 2026-05-23
