# FINDINGS — Solar Flare → Grid Coupling

**Diatom Sky R&D · Defensive Publication**
**Author:** KhaiB10
**Date:** 2026-05-23
**Status:** Open-data exploratory analysis. Not an operational risk model.
**License:** CC0 1.0 / MIT (see LICENSE)

---

## Summary in one paragraph

Using 94 years of open three-hour geomagnetic data from the GFZ Potsdam Kp/ap record (1932–2025, n=274,672 observations), we fit a Peaks-Over-Threshold Generalized Pareto Distribution (GPD) to daily-maximum ap values above the 95th percentile and ran a 20,000-trial Monte Carlo simulation of decade-scale worst-storm magnitudes. The fitted tail is **near-exponential (shape ξ = -0.024)** and yields a per-decade probability of **~58.5%** of at least one "Carrington-class" event (daily ap ≥ 400, equivalent to a Kp = 9 / G5 day). Cross-referencing the seven best-documented modern geomagnetic grid impact events shows that every documented G5 day in the modern record coincided with measurable grid-side impact, ranging from the 1989 Quebec collapse to the 2024 Gannon storm's ~$500M US precision-agriculture GPS disruption.

## Why publish this

The hazard rate is not the hard part — NOAA, NERC, and academic groups have published similar tail estimates. The contribution here is:

1. **Fully reproducible from a single ~16 MB open file** (`Kp_ap_since_1932.txt` from [GFZ Potsdam](https://kp.gfz.de/)) with one Python script.
2. **The 2024 Gannon storm is now in the historical record**, and the Monte Carlo trial described here was re-baselined to include it. This is one of the first open replications that bakes 2024 into the decadal probability.
3. **Documented grid-side outcomes** are tabulated alongside the modeled hazard to make the conditional impact discussion concrete rather than abstract.

We make no claims about utility-specific risk, transformer reliability, or required mitigation spend. Those questions require utility-side data we do not have.

## Method (short version)

| Step | What we did |
|------|-------------|
| 1 | Load GFZ Potsdam `Kp_ap_since_1932.txt` (3-hour Kp & linear `ap` index) |
| 2 | Filter to 1932-01-01 → 2025-12-31, drop placeholder records (Kp < 0) |
| 3 | Aggregate to daily statistics (max/mean Kp and ap) |
| 4 | Classify each day on the NOAA G-scale (G0–G5) |
| 5 | Fit `scipy.stats.genpareto` to ap-max excesses above the 95th percentile (ap = 80) |
| 6 | Monte Carlo: draw Poisson(λ·10) exceedances per decade; sample magnitudes from the fitted GPD; record worst per trial; repeat 20,000× |
| 7 | Overlay seven documented modern GIC grid events onto the Kp/ap timeline |

Full code: [`scripts/analyze.py`](scripts/analyze.py). Re-running on a fresh checkout produces identical numbers (seed = 20260523).

## Headline numbers

| Quantity | Value |
|---|---|
| Records analyzed | 274,672 three-hour Kp values |
| Days analyzed | 34,334 |
| Years covered | 94 (1932–2025) |
| GPD threshold (95th percentile of daily ap) | **ap = 80** |
| GPD shape ξ | **−0.024** (near-exponential) |
| GPD scale σ | 65.3 |
| Exceedance rate above threshold | **λ = 16.0 / year** |
| P(≥1 ap ≥ 400 in a decade) | **0.585** |
| P(≥1 ap ≥ 600 in a decade) | 0.025 |
| Median worst-decadal ap | 413 |
| 95th-percentile worst-decadal ap | 561 |
| 99th-percentile worst-decadal ap | 650 |

The near-zero shape parameter is notable: it implies the tail is exponential rather than heavy. Geomagnetic storms appear to have a soft ceiling within the observed range, but the 1989 and 2003 events sit on the right edge of that tail — and **the model puts a non-trivial mass beyond observed**.

## The seven documented modern GIC grid events (overlay)

| Date | Event | Kp_max | ap_max | Min Dst (nT) | Grid impact |
|------|-------|--------|--------|--------------|-------------|
| 1989-03-13 | Quebec Blackout | 9.0 | 400 | -589 | Full Hydro-Québec collapse; 6M people, 9h outage; transformer damage |
| 2003-10-29 | Halloween Storm | 9.0 | 400 | -383 | Malmö, Sweden 50-min blackout; Eskom transformer damage in S. Africa |
| 2003-10-30 | Halloween (continued) | 9.0 | 400 | -401 | Continued GIC; multi-utility HV transformer degradation |
| 2015-03-17 | St. Patrick's Day | 7.7 | 179 | -223 | GIC observed in Finnish/US pipelines & transformers; no outages |
| 2017-09-08 | September 2017 | 8.3 | 236 | -142 | GIC during hurricane season; minor transformer heating |
| 2024-05-10 | **Gannon (Mother's Day)** | 8.7 | 300 | -412 | First G5 since 2003; ~$500M US precision-ag GPS disruption; no major blackouts |
| 2024-10-10 | October 2024 G4 | 8.7 | 300 | -335 | Strong GIC observed; minor utility-side anomalies |

Sources: [NASA / Boteler 2019](https://doi.org/10.1029/2019SW002278), [Pulkkinen et al. 2005](https://doi.org/10.1029/2004SW000123), [Kappenman 2010 (Metatech / NASA)](https://www.swpc.noaa.gov/sites/default/files/images/u33/Geomagnetic%20Storms.pdf), [NOAA SWPC archives](https://www.swpc.noaa.gov/).

## Interpretation

- **Decadal G5 odds are roughly coin-flip-favorable.** The 58.5% figure should not feel surprising — in the 94-year record there are several G5 days, and the fitted Poisson rate above the 95th percentile is ~16/year (almost all of those are G1–G3 noise; the heavy tail makes G5 plausible per decade).
- **The shape parameter ξ ≈ 0 is consistent with prior literature** (e.g. [Tsubouchi & Omura 2007](https://doi.org/10.1029/2007SW000329) found ξ in the same neighborhood). Our 2025-updated fit does not push toward a heavier tail.
- **Observed grid impact scales sharply non-linearly with ap.** A Kp = 7.7 storm (St. Patrick's 2015) produced GIC observations but no outages. A Kp = 9 storm (Quebec 1989) collapsed an entire interconnection. There is no public dataset rich enough to fit a smooth dose-response curve on this gap; that remains an open question for utility-side researchers.
- **The 2024 Gannon storm landed inside the predicted "expected once per decade" window** and largely confirms the conventional decadal hazard estimate while exposing a new vector (GPS / precision agriculture) that pre-2024 risk frameworks had not emphasized.

## What this analysis does NOT do

- Does not quantify per-utility GIC exposure (which depends on ground conductivity, transmission topology, and transformer design).
- Does not assert any specific mitigation policy, spending requirement, or compliance posture for any operator.
- Does not project solar-cycle-conditional rates (i.e. odds within Cycle 25 vs. Cycle 26 specifically); the Poisson rate is treated as homogeneous in time.
- Does not extrapolate beyond the 1859 Carrington event magnitude (commonly estimated at ap ≳ 1000–1700) — that lives in a regime where 94 years of data cannot constrain the tail.

## Reproducibility

```bash
git clone https://github.com/KhaiB10/solar-flare-grid-coupling
cd solar-flare-grid-coupling
pip install numpy pandas matplotlib scipy
curl -L -o data/Kp_ap_since_1932.txt https://kp.gfz.de/app/files/Kp_ap_since_1932.txt
python scripts/analyze.py
```

Outputs land in `figures/` and `data/derived_*.csv`. Seed is fixed at `20260523`.

## Open questions worth follow-up

1. Can the GFZ Kp record be extended pre-1932 using auroral observations + magnetometer fragments to better constrain the upper tail?
2. Does conditioning the Poisson rate on the 11-year solar cycle phase (rising, max, declining, min) sharpen the decadal estimate?
3. What is the right open-data proxy for grid-side severity that survives NERC's CEII redaction? OE-417 is partial; EAGLE-I outage data is post-hoc and rarely tagged "space weather."

## Citation

If you use this analysis or its derived tables, please cite:

> KhaiB10 (2026). *Solar Flare → Grid Coupling: a 94-year open replication.* Diatom Sky R&D. https://github.com/KhaiB10/solar-flare-grid-coupling

## Data & references

- [GFZ Potsdam Kp/ap historical record (1932–present)](https://kp.gfz.de/) — primary dataset
- [WDC Kyoto Dst index](https://wdc.kugi.kyoto-u.ac.jp/dstdir/) — used for event-table cross-validation
- [NOAA SWPC](https://www.swpc.noaa.gov/) — G-scale definitions, event archives
- [DOE OE-417 Annual Summaries](https://www.oe.netl.doe.gov/OE417_annual_summary.aspx) — referenced; not used directly (sparse space-weather tagging)
- Boteler, D.H. (2019). *A 21st Century View of the March 1989 Magnetic Storm.* Space Weather, 17(10).
- Pulkkinen et al. (2005). *Geomagnetic storm of 29–31 October 2003: Geomagnetically induced currents and their relation to problems in the Swedish high-voltage power transmission system.* Space Weather, 3(8).
- Kappenman, J. (2010). *Geomagnetic Storms and Their Impacts on the U.S. Power Grid.* Metatech Corp. / NASA.
- Tsubouchi & Omura (2007). *Long-term occurrence probabilities of intense geomagnetic storm events.* Space Weather, 5(12).
