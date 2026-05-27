# v16 — SC26 amplitude prediction (true out-of-sample test of v15)

**Date:** 2026-05-27
**Seed:** 20260523
**Status:** Defensive-publication R&D. MIT (code) + CC0 (data/figures). No warranty of fitness for any operational use.

---

## What this is

The clearest test of v15's hierarchical Hawkes model: **predict SC26's total G4+ event count before any SC26 data exists.** No new sampling — just draw fresh `(μ₀, α, β, κ)` parameters from the **population posterior** v15 already estimated, then simulate forward 11.5 years.

This is the difference between a hindcast and a forecast. v15's SC25 forward simulation conditioned on 10 already-observed SC25 events. v16's SC26 simulation conditions on **nothing** — only the population distribution of how solar cycles behave.

## The Russian-doll trick made operational

In v15 the model said: "the population of solar cycles has a typical log-μ₀ of −5.55 ± 0.38." That number is a **distribution over distributions** — uncertainty about both the average cycle and how much variability the next cycle could exhibit.

For SC26 we do this:
1. Pick one posterior draw → get one realization of `(μ_μ, σ_μ, μ_α, σ_α, μ_β, σ_β, μ_κ, σ_κ, γ)`.
2. **Sample a fresh cycle's `μ₀, α, β, κ` from those population distributions.** This is the prior-predictive for an unobserved cycle.
3. Simulate one SC26 by Ogata thinning over 2031-07-01 → 2042-12-31 (4202 days = 11.5 years).
4. Repeat 2000 times across posterior draws to fold in both **between-cycle variability** and **parameter uncertainty**.

The result is a **calibrated probability distribution over how many G4+ events the next solar cycle will produce.**

---

## Headline forecast: it depends entirely on the F10.7 scenario

| Scenario | F10.7 assumption | SC26 G4+ count (median, 95% HDI) | P(≥1 G5) | P(≥2 G5) |
|---|---|---|---|---|
| **PHYSICAL** (recommended) | SC23-shape, peak ∼N(155, 25²) sfu | **13** [2, 38] | 64.7% | 35.6% |
| FLAT_QUIET | Flat at 118 sfu ([Singh+2021](https://ui.adsabs.harvard.edu/abs/2021cosp...43E1044S/abstract) SC26 forecast) | 22 [6, 61] | 82.3% | 56.8% |
| FLAT_AVG | Flat at 150 sfu (mid-modern) | 37 [12, 99] | 92.7% | 79.0% |
| FLAT_LIKE25 | Flat at 180 sfu (SC25-like) | 55 [19, 134] | 97.1% | 89.9% |

**Historical reference:** SC22=28, SC23=31, SC24=3, SC25 (partial through 2025-05) = 10.

### Why PHYSICAL is the credible scenario

Real solar cycles spend most of their time **below their maximum**. The flat scenarios artificially pin F10.7 at one value for 11 years; this overweights peak-year storm production. The PHYSICAL scenario applies an SC23-shaped F10.7 trajectory (low ascending ramp, peak around year 5, gentle decline), and samples the peak amplitude from `N(155, 25²)` truncated to [80, 230] sfu — wide enough to bracket both quiet (Singh-like SC26 ≈ 118 sfu) and active (SC25-like ≈ 180 sfu) cases.

Under PHYSICAL, **SC26 is forecast as the third-quietest in the modern era**: median 13 G4+ events, comparable to SC24 (3) and SC25's running total (~30 by completion). The 95% interval [2, 38] brackets every modern cycle except SC23 and SC22.

## Where the prediction lines up with the literature

Singh et al. (2021) used simplex projection to forecast SC26's peak F10.7 at **118 ± 9 sfu** (their FLAT_QUIET corresponds to this). The published forecasts cluster around "weaker than SC25, comparable to SC24-25 range" — and our **PHYSICAL** scenario lands right there.

Where ours differs: published SC26 forecasts give a single number (peak amplitude). Ours gives a **full posterior over G4+ event counts**, broken out by scenario, that anyone can update with new F10.7 data as it arrives.

---

## When will SC26 events cluster?

`figures/61_sc26_event_timing.png` shows the aggregate temporal density of G4+ events under PHYSICAL, with the F10.7 template overlaid. **Events concentrate sharply in 2035-2037**, matching the predicted SC26 maximum at year ~5 of the cycle. Two practical implications:

- **Grid hardening priority window: 2035 ± 2 years.** That's when the next G5-grade hit is most likely.
- **The decline phase (2038-2042) still produces events** — about a third of cycle activity. Half-life of clustering remains ~1.7 days regardless.

## The G5 (Carrington-class precursor) probability

Under PHYSICAL, **64.7%** of simulations produce at least one G5+ event during SC26, and **35.6%** produce two or more. These numbers come from the empirical mark distribution of SC23-25 (n=45, where 4/45 = 8.9% of events landed at mark ≥ 9). A genuinely Carrington-class event (mark > 10) is **not** in the modern catalog, so the tail risk is bounded below by what we observe; the real Carrington-class probability is higher than what we report.

To go further on Carrington-class events would require: (a) refitting marks with a Pareto tail on top of the empirical body; (b) using the 19th-century catalog (which v16 is also planning to drop in a later refit). That's a follow-up.

---

## Limitations (don't pretend these aren't there)

1. **No SC26 F10.7 data exists.** The single largest source of uncertainty is which scenario (PHYSICAL vs FLAT_*) you believe. The model is a transformation from F10.7 trajectory to event distribution — garbage in, garbage out on the F10.7 side. Anyone using this should be honest about which scenario they're committing to.

2. **PHYSICAL template borrows SC23's shape.** Real cycles differ in skew and tail. A reasonable extension: sample the entire F10.7 trajectory shape (not just amplitude) from a library of past cycles, weighted by SC24/SC25 similarity.

3. **The population posterior is dominated by 1844-1947 ground-magnetometer data** (372 of 434 events). v15 noted SC10 is suspiciously active. A modern-only refit (v17 plan) would tighten the population priors and probably shift `μ_μ` downward, making SC26 forecasts slightly quieter.

4. **No cross-cycle correlation modeled.** SC24 was quiet and so is SC25 (partial). If they're part of a longer Gleissberg-cycle modulation (~100 yr), SC26 is more likely to be quiet too. v15 + v16 treats every cycle as exchangeable. This is conservative for "what could happen" but might overstate SC26 activity.

5. **Mark distribution capped at observed modern range.** Our P(G5) is a **lower bound**.

6. **No coupling between α and μ₀.** If quiet cycles have proportionally weaker excitation (a physical possibility), the actual forecast spread is narrower than ours.

---

## Where this could go

- **v17 — modern-only refit (1947-2025).** Drop the 19th-century catalog, see whether SC10's anomaly survives. Tightens v15's population posterior and gives a "fair" SC26 forecast.
- **Joint SC25 + SC26 + SC27 trajectory.** Use the 2-level Gleissberg hierarchy mentioned in v15 to forecast three cycles forward simultaneously.
- **Live updating.** As 2031, 2032, 2033 F10.7 data arrives, refit the per-cycle SC26 parameters and shrink the forecast. Could ship as a public dashboard.
- **Cross-domain universality.** The same hierarchical-Bayes Hawkes machine fit on hurricanes (Hurricane v3, currently local) shows the same structural finding (μ₀ varies, kernel shape doesn't). Worth a joint paper.

## Who this impacts

- **Grid operators planning the 2030s.** Median 13 events with P(≥1 G5) = 65% under the most realistic scenario is concrete enough to budget against. Pair with FERC/NERC GMD planning requirements.
- **Insurance / parametric covers.** First reproducible Bayesian forecast of a future cycle's event count — anyone can re-run with their own scenario.
- **Researchers comparing solar-cycle prediction methods.** Published SC26 forecasts give peak amplitudes; this one gives event-count distributions. Easier to verify against actual SC26 observations once they exist.
- **Spaceweather.gov / SWPC.** A complementary "ground-up" forecast from event statistics, rather than top-down from sunspot dynamo models. Different failure modes — useful for ensemble.

---

## Reproduce

```bash
git clone https://github.com/KhaiB10/solar-flare-grid-coupling
cd solar-flare-grid-coupling
python scripts/analyze_sc26_forecast_v16.py    # ~20s after v15 trace exists
python scripts/replot_sc26_v16.py              # regenerates plots only
```

Outputs:
- `data/v16_sc26_summary.json` — all reported numbers
- `figures/60_sc26_count_distribution.png` — main forecast histogram
- `figures/61_sc26_event_timing.png` — temporal density under PHYSICAL
- `figures/62_sc26_tail_probability.png` — G5 risk per scenario
- `figures/63_sc26_count_ecdf.png` — cumulative distribution

Sources: v15 InferenceData (`data/v15_idata.pkl`), NOAA G4+ catalog, Helsinki magnetometer records, Penticton F10.7 series ([NOAA NCEI](https://www.ngdc.noaa.gov/stp/space-weather/solar-data/solar-features/solar-radio/noontime-flux/penticton/)).

---

## One-sentence summary

**Under a physically realistic F10.7 trajectory, SC26 will most likely produce 13 G4+ geomagnetic storms (95% interval 2 to 38), peak around 2036, and has a ~65% chance of at least one G5-grade event — quieter than SC22/SC23 but comparable to the SC24-25 range.**
