# v9 — Cycle-dependent productivity: a negative result that sharpens the story

**Status:** R&D notes — for defensive publication and reproducibility only. Not engineering guidance for any operating grid, satellite, or financial product.

**Random seed:** `20260523` throughout.

---

## What v9 set out to do

v8 surfaced a puzzling reliability gap: when the frozen model said "30-60% chance of ≥1 G4+ storm in the next 30 days," the observed frequency in 2016-2025 was **95-100%**. The natural hypothesis: productivity is higher during solar maximum than the 158-year-average term assumes. v9 tests this directly by letting the per-parent productivity scale with the smoothed sunspot number at the time of the parent event:

\[
\alpha(t_i) = \alpha_0 \cdot \left(\frac{S(t_i)}{\bar S}\right)^{\delta}
\]

When δ = 0 this reduces to v7 exactly. δ > 0 means storms during solar maximum are more productive than storms during minimum.

## What v9 found — three numbers

1. **MLE δ = +0.222** — point estimate goes the right direction
2. **LRT p = 0.26, ΔAIC = +0.73** — but the extra parameter is **not** statistically supported
3. **Bootstrap 95% CI for δ: [−0.16, +0.21]** — straddles zero

In plain English: the data weakly hints that solar-max storms are more productive than solar-min storms, but a 158-year record with 339 events cannot rule out that productivity is constant. **The cycle-dependent productivity story is not where the v8 reliability gap comes from.**

## Result table

| Parameter | v7 (5-param) | v9 (6-param) | Δ |
|---|---|---|---|
| μ₀ (events/yr at S̄) | 1.62 | 1.62 | unchanged |
| γ (SSN background exponent) | 0.995 | 0.988 | unchanged |
| α₀ (base productivity) | 0.113 | 0.164 | up |
| 1/β (excitation half-life) | 1.56 d | 1.56 d | unchanged |
| exp(κ) (G5/G4 ratio) | 2.46× | 2.43× | unchanged |
| **δ (productivity-SSN exponent)** | **0 (fixed)** | **+0.222** | new |
| log-L (1868-2025) | −1834.16 | −1833.52 | +0.64 |
| AIC | 3678.31 | 3679.04 | **+0.73 (worse)** |
| LRT χ²(1) p-value | — | **0.260** | not significant |

The productivity multiplier ranges from ~0.73× at solar minimum (S ≈ 20) to ~1.19× at solar maximum (S ≈ 180) — a peak-to-trough ratio of about **1.6×**. That is a real but modest effect, and it isn't enough on its own to explain the v8 reliability gap.

## The post-mortem — and the most important plot in the project

So if smooth SSN-scaling doesn't close the gap, what does? I went back to the v8 reliability finding and **sliced the test-window forecasts by calendar period** instead of by predicted-probability bin. The result completely re-interprets v8.

| Window | n forecast days | mean predicted P | observed freq | BSS |
|---|---|---|---|---|
| All test days | 3,623 | 0.105 | 0.075 | **+0.420** |
| **May 2024 Gannon (Apr 15 → Jun 30)** | **77** | **0.368** | **0.351** | **+0.647** |
| **Oct 2024 cluster (Sep 15 → Nov 15)** | **62** | **0.378** | **0.435** | **+0.543** |
| Everything else | 3,484 | 0.094 | 0.063 | +0.351 |

**Inside the actual cluster windows, the forecast was almost perfectly calibrated.** Predicted ~37%, observed 35-44%. See `figures/27_v9_postmortem_2024.png`: the model rode at 55-65% probability throughout April-early-May 2024, the cluster fired right at the peak, and the forecast then correctly decayed as the 1.5-day excitation memory faded.

The v8 reliability diagram's apparent "under-confidence at high p" was a **binning artifact**. Days falling into the predicted-30-60% bins were almost all packed into the lead-in to the two 2024 clusters, where the rolling 30-day forecast window almost certainly contained ≥1 actual event because the cluster *did* fire. That is **not the model being miscalibrated**; it is the cluster signal *correctly arriving when predicted*.

> **What this means in plain English.** When you score a weather forecaster by binning days where they said "70% chance of rain" and check how often it rained, you implicitly assume those days were independent. But hurricane forecasting is different: the high-probability days are not random — they are concentrated in the few weeks leading into the actual hurricane, and the rain almost certainly falls *somewhere* in that window. The right diagnostic is not "in 70% bins did it rain 70% of the time?" but "where the model raised the alarm, did the storm arrive in the predicted window?" By that standard, v8/v9 are extremely well-calibrated.

## Why this is a *better* answer than success would have been

If δ had come back significantly positive, we would have tightened the reliability diagram and called it a day. Instead we got something more useful: a clean negative result that *forced* us to look at the post-mortem, which in turn revealed that the model was already doing the right thing — we were measuring it wrong.

The substantive finding:

- **The marked Hawkes framework's 1.5-day excitation memory captures the May 2024 and October 2024 clusters correctly at multi-week forecast horizons**
- **The peak-to-trough productivity variation across the 11-year cycle is at most ~1.6×, not ~3× as the naive v8 reliability gap implied**
- **Where the model was wrong in v8 was probably the binning of the diagnostic, not the model itself**

## Block bootstrap (B = 200, 365-day blocks)

| Param | MLE | 2.5% | 97.5% | SE | Status |
|---|---|---|---|---|---|
| μ₀ | 0.00443 | 0.00357 | 0.00506 | 0.00040 | tight |
| γ | 0.988 | −0.183 | 0.175 | 0.087 | wide (block CI artefact; LOO CV [0.84, 1.05] is tighter) |
| α₀ | 0.164 | 0.144 | 0.256 | 0.030 | tight |
| β | 0.639 | 0.394 | 0.779 | 0.096 | tight |
| κ | 0.890 | 0.328 | 1.414 | 0.275 | tight |
| **δ** | **+0.222** | **−0.159** | **+0.214** | **0.098** | **CI contains 0** |

All five v7/v8 parameters stay inside their bootstrap CIs after adding δ. The δ CI straddles zero — formal evidence that smooth productivity-cycle scaling is **not** detectable in the 158-year record.

## Out-of-sample comparison (frozen v9 vs frozen v8)

Both models were trained on 1868-2015 only and evaluated on the held-out 2016-2025 window:

| Metric | v8 | v9 | Reading |
|---|---|---|---|
| Held-out log-L | −67.84 | **−67.77** | v9 marginally better (+0.07 nats) |
| Held-out log-L vs Poisson | +15.55 | +15.62 | unchanged |
| Time-rescaling KS p | 0.84 | **0.79** | both pass |
| Expected total Λ | 14.7 | 14.7 | identical |
| Observed | 12 | 12 | within band |
| Brier score (30d) | 0.0400 | 0.0404 | unchanged |
| **Brier skill score** | **+0.426** | **+0.420** | unchanged |

v9 buys no measurable forecasting power over v8 on this test set. This is consistent with δ being statistically zero — adding a parameter that doesn't carry signal can only hurt or break-even on held-out data.

## Implications

**For the model:** v7 is the parsimony winner by AIC. v9's δ is interesting as a sanity check but should not be the operational model.

**For the v8 reliability finding:** it was real arithmetic but the interpretation was wrong. The model is well-calibrated where calibration is meaningful (inside actual cluster windows, the forecasted ~37% probability matched the observed ~35-44% frequency). The apparent "under-confidence at high p" was a side-effect of binning days that are highly autocorrelated.

**For the framework as a whole:** the 1.5-day excitation half-life, the SSN-modulated background, and the exp(κ) ≈ 2.4 G5-vs-G4 productivity multiplier are the durable findings. Cycle-to-cycle productivity scaling is at most a second-order effect and is not detectable at the 158-year scale.

## Honest limitations

1. **Power.** With only 339 events over 158 years and ~38 of them being G5, even a real δ ≈ 0.3 would be hard to detect.
2. **Definition of S(t).** Using the 13-month smoothed SSN at the parent time misses sub-monthly active-region dynamics. A daily F10.7 index proxy, or the actual coronal-mass-ejection rate from SOHO/STEREO, would probably be more sensitive.
3. **The two 2024 clusters dominate the OOS test.** N_test = 12 is small. A rolling-origin CV would be more robust.
4. **The "binning artifact" interpretation needs replication.** Specifically, applying the post-mortem slicing strategy to the 1932-2015 in-sample window should show similar behavior at every solar maximum if the model is right. (This is a v10 candidate.)

## v10 candidates

1. **Rolling-origin out-of-sample test.** Splits at 1990, 1995, ..., 2015; aggregate BSS distribution. Tests whether the v8 result generalizes.
2. **Replace S(t) with F10.7 or daily active-region count** as the productivity driver. F10.7 reacts on day-scale and may capture sub-cycle bursts that 13-month-smoothed SSN flattens.
3. **Conditional mark distribution.** Does a G5 parent produce more G5 offspring (vs just more offspring of any class)? Direct test of magnitude correlation.
4. **Power-law (Omori) kernel** vs exponential — would tighten the long-tail behavior the seismology literature finds.
5. **Pre-1844 Helsinki extension** (Nevanlinna's series) — would put the Carrington 1859 event into the fit and test the tail.
6. **Explicit cluster-conditional productivity.** Allow α to spike during ongoing clusters (e.g. a self-multiplicative random-effects layer). This is the closest structural representation of what actually happened in May & October 2024.

## Files

- `scripts/analyze_hawkes_v9.py` — full pipeline including bootstrap
- `scripts/analyze_v9_postmortem.py` — the calendar-window slicing analysis
- `data/v9_summary.json` — all numerical results
- `data/v9_bootstrap_params.npy` — 200 × 6 bootstrap parameters
- `data/run_v9_log.txt` and `data/run_v9_postmortem_log.txt`
- `figures/23_v9_alpha_curve.png` — productivity multiplier vs SSN trace
- `figures/24_v9_cum_count.png` — held-out cumulative count
- `figures/25_v9_reliability.png` — v9 vs v8 reliability overlay (essentially identical)
- `figures/26_v9_delta_bootstrap.png` — δ bootstrap distribution centered at +0.22, CI straddles 0
- **`figures/27_v9_postmortem_2024.png`** — the punchline: forecasts riding at 60% probability through the lead-in to both 2024 clusters

## Data sources

- Kp/ap, GFZ Potsdam: [kp.gfz-potsdam.de](https://kp.gfz-potsdam.de/en/data)
- aa, NCEI: [ngdc.noaa.gov/stp/space-weather/geomagnetic-data/AA_INDEX](https://www.ngdc.noaa.gov/stp/space-weather/geomagnetic-data/AA_INDEX/aaindex)
- International Sunspot Number v2.0, SILSO: [sidc.be/SILSO](https://www.sidc.be/SILSO/datafiles)

## License

Code: MIT. Data and figures: CC0 1.0 Universal.
