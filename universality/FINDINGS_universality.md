# Hawkes universality: solar storms ⋃ hurricanes ⋃ earthquakes

**Branch:** `universality` of `solar-flare-grid-coupling`
**Date:** 2026-05-27
**Random seed:** 20260523
**License:** MIT (code) · CC0 (data/figures)
**Status:** Defensive publication — Diatom Sky R&D, exploratory. Numbers are honest; framing is humble.

---

## 0. TL;DR (60 seconds)

Three completely different physical systems — solar coronal storms (days), Atlantic + East-Pacific hurricanes (weeks), and tectonic earthquakes (minutes-hours) — when fit with the **same** statistical model (a Hawkes self-exciting point process with marks), produce parameters that **share the same qualitative structure**, despite operating ~7 orders of magnitude apart in their characteristic forcing timescale.

| Domain | n (branching ratio) | t_half (kernel half-life) | κ (mark productivity) | Sub-critical? | Mark-positive? |
|---|---|---|---|---|---|
| Solar G4+ storms (this work, v15) | **0.134** [0.078, 0.241] | **1.16 d** | **+1.06** / Kp unit | ✅ | ✅ |
| Hurricanes Cat-3+ (this work, v3) | **0.006** [0.002, 0.012] | **17.4 d** | **+0.17** / Saffir-Simpson cat | ✅ | ✅ |
| Earthquakes ETAS (literature)    | **0.5 – 0.8**           | **3 min – 7 h**          | implicit (productivity ∝ 10^(α·M)) | ✅ (typically) | ✅ |

**What's actually universal:** the *structural form* — sub-critical clustering (n < 1), monotonically-decreasing kernel, positive mark sensitivity (bigger trigger ⇒ more aftershocks on average). Not the absolute magnitudes.

**What is NOT universal:** the dimensionless ratio t_half / τ_forcing. A naive scale-invariant ansatz (t_half ∝ τ_forcing) fails. OLS slope across 3 domains is **0.24**, not 1.0. This was a hypothesis we tested and falsified.

---

## 1. The popcorn-machine analogy — extended to three machines

In the solar v12 writeup we framed Hawkes as a **popcorn machine**: every kernel pops, and a fraction of pops trigger neighbors. Two parameters control everything visible:

- **n (branching ratio)** — average number of *child* pops triggered per parent. n < 1 ⇒ the machine eventually quiets; n ≥ 1 ⇒ runaway.
- **t_half (kernel half-life)** — how quickly the trigger probability fades. Short half-life ⇒ tight clusters; long half-life ⇒ smeared echoes.

Now we have three machines, and we want to know whether they're the *same machine in different rooms* or *different machines that happen to all pop*.

- **Solar machine.** Pops are G4+ geomagnetic storms. Each pop produces ~0.13 child pops within a kernel half-life of ~1.2 days. Bigger pops (Kp 9) are mildly more productive than smaller pops (Kp 8): κ ≈ +1.06 in the natural log-productivity scale.
- **Hurricane machine.** Pops are Atlantic + East-Pacific Cat-3+ landfalls/intensifications. Each pop produces only ~0.006 child pops, but the kernel half-life is ~17 days — long enough that the cluster covers a meaningful chunk of a hurricane season. Bigger pops (Cat-5) are also more productive: κ ≈ +0.17 per Saffir-Simpson category.
- **Earthquake machine.** Pops are mainshocks above Mmin. ETAS literature (Ogata 1988, Utsu 1993, USGS) reports n in [0.5, 0.8] depending on Mmin, t_half in minutes-to-hours via the Omori c·(2^(1/p)−1) law, and the well-known 10^α·M productivity scaling — exactly mark-positive in the same direction.

All three machines are **sub-critical** (n < 1, no runaway), **mark-positive** (bigger trigger ⇒ more aftershocks), and **monotonically-decaying** kernels. They differ enormously in:
- the *level* of n (solar/hurricane = quiet machines, earthquakes = noisy machines),
- the *timescale* of t_half (minutes vs days vs weeks),
- the *physical mechanism* (Earth's magnetotail vs ocean-atmosphere coupling vs crustal stress transfer).

But once you remove the units and ask "is this thing self-exciting in the same *qualitative* way as the others?", the answer is yes for all three.

---

## 2. What we did

### 2.1 Inputs

- **Solar v15:** PyMC NUTS posterior (4 chains × 1500 draws), 14 solar cycles SC10–SC25, 435 events, hierarchical α/β/κ per cycle with population hyperparameters. R̂ = 1.000, 0 divergences. (Committed to `main` at 7b4c821.)
- **Hurricane v3:** 200 stratified bootstrap fits to NHC HURDAT2 1949–2024, 452 Cat-3+ events, Atlantic + East-Pacific basins, β=0.0388/d, κ=+0.180 in Saffir-Simpson units. (Lives in a private local repo — not committed here; results referenced through `data/hurricane_v3_summary.json`.)
- **Earthquake ETAS:** Literature numbers only ([USGS Appendix S](https://pubs.usgs.gov/of/2013/1165/pdf/ofr2013-1165_appendixS.pdf), [Utsu 1993 slides](http://www-solid.eps.s.u-tokyo.ac.jp/~hassei/2017/slides/1b.pdf)). No refit here — we use the published p, c, n ranges as a third domain anchor.

### 2.2 Pipeline

`universality/scripts/build_universality_table.py`:

1. Pull solar v15 InferenceData. Draw B = 3000 "typical cycle" parameter triples from the population posterior: `α_typ ~ Lognormal(μ_α, σ_α)`, `β_typ ~ Lognormal(μ_β, σ_β)`, `κ_typ ~ Normal(μ_κ, σ_κ)`. This represents what a *new, unobserved* solar cycle would look like under the hierarchical model.
2. Pull hurricane v3 bootstrap params (shape (200, 7): μ_atl, μ_epac, b_soi, b_amo, α, β, κ).
3. For each domain, compute the branching ratio `n = α · E[exp(κ·(m − m₀))]` where the expectation is taken over the empirical mark distribution. For solar, m₀=8 (Kp G4 threshold) and marks are Kp magnitudes in 8.0–9.5. For hurricanes, m₀=3 (Cat-3) and marks are `peak_cat` in {3, 4, 5}.
4. Compute kernel half-life `t_half = ln(2) / β` (exponential kernels). For ETAS, use the Omori power-law half-life `c · (2^(1/p) − 1)` with the Utsu 1993 ranges.
5. Render three-panel figure (branching ratio, half-life, κ) plus a scaling figure (t_half vs τ_forcing on log-log).

### 2.3 The scale-invariance test

If Hawkes parameters were genuinely scale-invariant in time, we'd expect `t_half ∝ τ_forcing` (slope = 1 on log-log axes). We tested this against three honestly-chosen forcing timescales:
- Earthquakes — S-wave traversal of a fault patch, ~10 s
- Solar — Carrington rotation, ~27.3 d
- Hurricanes — MJO oscillation, ~45 d

`universality/scripts/scaling_test.py` reports:

```
slope b = 0.24   (b=1 = perfect scale-invariance)
inter. a = 0.29
R²       = 0.73
ratio range: 0.39 (hur) to 0.04 (solar) to 11500 (quakes)
```

**Verdict:** the strict scaling hypothesis fails. The earthquake t_half is **vastly larger** than the S-wave forcing scale (≈10⁴×), while the solar and hurricane half-lives are *fractions* of their forcing scales. The "5% of forcing scale" heuristic I drew in the first draft of figure 2 was misleading; it has been replaced with the actual OLS fit and a `slope=1` reference for honesty.

This is fine. The universality claim was never "all kernels are dimensionally identical" — it was "all three systems have the same statistical anatomy." That anatomy survives. The dimensionless ratio does not.

---

## 3. Results

### 3.1 Branching ratio (Panel A, `figures/01_universality_three_panel.png`)

- **Solar G4+**: posterior median 0.134, 80% CI [0.078, 0.241]. Roughly **1 in 7** observed storms is, in expectation, a Hawkes-triggered child of an earlier storm.
- **Hurricanes Cat-3+**: 0.006 [0.002, 0.012]. The hurricane catalog is **dominantly background-driven** — almost every major hurricane is a Poisson event modulated by SOI/AMO, not a triggered child of a prior major hurricane. This is consistent with the meteorological intuition that two Cat-3+ storms in the same basin within ~2 weeks is rare even in active years.
- **ETAS quakes**: 0.5 – 0.8 from the literature. Earthquakes are the *noisiest* of the three — a large fraction of catalog events are children, which is why "aftershock" is a household word in seismology but not in heliophysics.

### 3.2 Kernel half-life (Panel B)

- **Solar**: 1.16 d. Tightly constrained.
- **Hurricanes**: 17.4 d.
- **ETAS**: 3 min – 7 h, depending on p, c. Tightly constrained at the *low* end of the figure's x-axis.

### 3.3 Mark productivity (Panel C)

- **Solar κ ≈ +1.06 / Kp unit** ⇒ a Kp-9 storm is `exp(1.06)` ≈ 2.9× as productive as a Kp-8 storm.
- **Hurricane κ ≈ +0.17 / category** ⇒ a Cat-5 is `exp(0.17·2)` ≈ 1.43× as productive as a Cat-3.
- **ETAS**: productivity scales as 10^(α·M) — equivalent to κ = α·ln(10) on the natural-log scale. Standard α ≈ 1 gives κ ≈ 2.3 / magnitude unit, so a magnitude jump of 1 ≈ 10× productivity. Mark-positive in the same direction.

### 3.4 Scaling figure (`figures/02_universality_scaling.png`)

Three points on log-log axes spanning 7 orders of magnitude in τ_forcing and 4 orders in t_half. OLS slope = 0.24, R² = 0.73. The fitted line is shown as a dotted reference; a slope=1 reference (true scale-invariance) is shown as a dashed reference. The data do **not** track slope=1, falsifying the strict scaling ansatz.

---

## 4. Where this could go

- **More domains.** Network-science Hawkes models (Twitter cascades, email storms, financial-market microstructure) have published n, β, and mark exponents. Adding 2–3 more domains would let us test whether the *structural* universality (sub-critical + mark-positive + decaying kernel) holds beyond geophysics, and whether there's a tighter sub-family with a real scaling law.
- **Why is n so different across domains?** A first-principles guess: n probably reflects the *coupling strength* of the driving medium. Solar wind streams clearly drive substorms (n ≈ 0.1 in our data) but most storms are independent CMEs. Hurricane-on-hurricane triggering is essentially zero on the Cat-3+ branch (Madden-Julian-style modulation dominates background, not aftershocks). Crustal stress transfer in earthquakes is geometric and direct — a real physical "kick" that propagates — hence n ≈ 0.5–0.8. A future writeup could relate domain-specific n to a measurable physical coupling parameter.
- **Live update.** Once SC25 enters its decay phase (2026–2028), we'll have new G4+ events to ingest. Re-running the v15 hierarchical fit each year would shrink the solar n posterior and let us watch n drift across the cycle.
- **Marks beyond category.** For hurricanes, replacing peak_cat with integrated kinetic energy (IKE) or ACE-per-storm might tighten the κ posterior considerably. The current Cat-3/4/5 mark is coarse.

## 5. Who this impacts

- **Solar / space-weather forecasters.** The branching ratio n ≈ 0.13 sets a hard ceiling on how much of the next G4+ event is predictable *from the prior G4+ event alone*. Most predictability still has to come from the F10.7-driven background μ(t). Don't sell triggering-based forecasts that imply more.
- **Insurance / catastrophe modelers.** Cat-3+ hurricane n ≈ 0.006 confirms that ENSO/AMO/SOI-driven background dominates seasonal counts at this severity threshold — triggered-cluster reserves are negligible for Cat-3+ at the basin level. This may *not* hold for lower thresholds (Cat-1+, tropical storms) where clustering literature reports more activity.
- **Anyone fitting Hawkes to a new dataset.** Use this table as a sanity-check prior: a fit returning n ≈ 0.95 (near-critical) or n > 1 (super-critical) for a *naturally-sparse* phenomenon (one G4+ per year, one Cat-5 every few years) is almost certainly mis-specified.

---

## 6. Caveats, things this is not, and known weaknesses

This is a **defensive publication**, not a peer-reviewed claim. We are putting numbers and code on the public record to (1) document the work, (2) invite criticism, and (3) avoid future patent overreach by anyone in the space. Specifically:

1. **The "third domain" is literature, not a refit.** We did not refit ETAS on our own earthquake catalog. The n ≈ 0.5–0.8 and t_half ≈ minutes-hours come from Ogata 1988, Utsu 1993, and the USGS UCERF appendix. They are well-attested, but they are not our numbers and our error bars on the green earthquake point in figure 2 reflect the literature range, not a fit.
2. **Hurricane data not in this repo.** The HURDAT2-derived Cat-3+ event catalog and the v3 bootstrap params live in a private local repo. The summary JSON committed here is sufficient to reproduce the universality numbers, but a full re-run of v3 from raw HURDAT2 is not possible from this branch alone. (Reasonable people can disagree about whether this is a feature or a bug.)
3. **Different mark units.** Solar Kp (8–9.5, 0.33-step ordinal) and Saffir-Simpson category (3–5, integer ordinal) and earthquake magnitude (continuous, log-scale already) are not the same kind of mark. κ values are *not* directly comparable across domains; we can only compare their **sign and qualitative magnitude**. We say "mark-positive across all three" — we do *not* say "0.17 < 1.06 < 2.3 implies anything."
4. **Exponential vs power-law kernels.** Solar and hurricane fits use exponential kernels (β·exp(−β·Δt)). ETAS uses Omori (c+t)^(−p). We compared "half-life," which is a common scalar both kernels expose, but the *shape* of the kernel beyond t_half differs — power-law kernels have heavier tails. A future cross-domain paper should compare full kernel CDFs, not just one summary statistic.
5. **Forcing timescales are debatable.** Choosing 10 s for earthquake forcing (S-wave traversal) vs Carrington rotation for solar (27.3 d) vs MJO for hurricanes (45 d) is one of several defensible choices. Different choices would shift the figure-2 points horizontally and possibly recover a different slope. We chose the **fastest** plausible forcing scale in each domain for consistency. Slope estimates are *not* robust to this choice and we don't claim they are.
6. **Solar branching ratio is on the population-typical cycle, not on a specific cycle.** Per-cycle n estimates from v15 vary from ~0.05 (quiet cycles SC23) to ~0.22 (active cycles SC19, SC22). The 0.134 number is the *expected next cycle* under the hierarchical model. Don't use it to forecast a *specific* cycle — use the per-cycle posteriors from `solar_v15_per_cycle.csv` directly.
7. **No formal universality test.** What would a real test of "structural Hawkes universality" look like? Probably a posterior-predictive check: simulate from each domain's fitted model, look at clustering statistics (Ripley's K, log-product spacings, residual point process), and ask whether the *transformed* point process is Poisson(1) in each domain. We did not do this. This is a 1-day cross-domain sketch, not a 1-year unification paper.
8. **One author, no peer review.** Everything in this branch is the work of one human and one Hawkes-friendly LLM, in one continuous coding session, with reproducible random seed 20260523. Treat it like any other arXiv preprint: useful starting point, not gospel.

---

## 7. Reproducibility

```bash
# from repo root, on branch `universality`
cd universality
# Requires solar v15 idata.pkl (10 MB) at one of:
#   universality/data/solar_v15_idata.pkl     (preferred, not committed for size)
#   /home/user/workspace/hawkes-universality/data/solar_v15_idata.pkl  (dev fallback)
# Requires hurricane v3 bootstrap params + event CSV at one of:
#   universality/data/hurricane_v3_bootstrap_params.npy
#   universality/data/hurricane_cat3plus_events.csv
PYTENSOR_FLAGS="cxx=,mode=FAST_RUN" python scripts/build_universality_table.py
python scripts/scaling_test.py
```

Outputs land in `universality/data/universality_summary.json`, `universality/data/scaling_test.json`, `universality/figures/01_universality_three_panel.png`, `universality/figures/02_universality_scaling.png`.

To re-derive the solar v15 idata.pkl: see `README.md` on `main` — run `python solar_flare_hawkes_v15.py`. To re-derive the hurricane v3 bootstrap: not currently public.

---

*— Diatom Sky R&D, exploratory work. Comments welcome via GitHub issues.*
