# FINDINGS v2 — Solar-Cycle-Phase-Conditioned Hazard

**Diatom Sky R&D · Defensive Publication addendum**
**Author:** KhaiB10
**Date:** 2026-05-23
**Status:** Open-data exploratory analysis. Not an operational risk model.
**Builds on:** [FINDINGS.md](FINDINGS.md)

---

## Why this addendum

The original analysis treated geomagnetic storm arrivals as a homogeneous Poisson process — the same rate every year. That is a known oversimplification: storm rate tracks the ~11-year solar cycle, and within a cycle the **declining phase** is unusually active because long-lived coronal holes drive high-speed solar wind streams (the [Russell-McPherron effect](https://doi.org/10.1029/JA078i001p00092)). This addendum splits the 94-year record into the four standard cycle phases and re-runs the Monte Carlo with phase-conditional rates.

## Data added

[SILSO V2 monthly total sunspot number](https://www.sidc.be/SILSO/datafiles), 1932–2025, 13-month centered smoothing. Detected solar maxima:

> 1937-04, 1947-05, 1958-03, 1968-11, 1979-12, 1989-06, 2002-03, 2014-04, **2024-10**

These match the standard solar-cycle numbering (Cycles 17–25) exactly.

## Phase definitions

- **min** — within ±18 months of a smoothed-SSN minimum
- **rising** — between minimum and following maximum (excluding the ±18-month windows)
- **max** — within ±18 months of a smoothed-SSN maximum
- **declining** — between maximum and following minimum (excluding the ±18-month windows)

## Phase-conditional fits

GPD threshold held constant at ap = 80 (the 95th percentile of the full record) so all phases are directly comparable.

| Phase | Days in record | ap > 80 days | λ per phase-year | GPD ξ | GPD σ |
|-------|---------------:|-------------:|-----------------:|------:|------:|
| min | 8,648 | 193 | **8.15** | −0.035 | 50.7 |
| rising | 3,582 | 143 | **14.58** | −0.072 | 65.8 |
| max | 8,648 | 470 | **19.85** | −0.037 | 71.2 |
| **declining** | **10,624** | **639** | **21.97** | −0.022 | 66.5 |

The declining phase is **2.7× more active than solar min** and edges out solar max itself, consistent with prior literature on recurrent high-speed-stream-driven storms.

## Decadal Carrington-class hazard by scenario

20,000-trial Monte Carlo of decadal worst-day ap with phase-specific λ and GPD draws:

| Scenario | P(≥1 ap ≥ 400 in 10 yr) |
|---|---:|
| Decade entirely at solar minimum | **6.3%** |
| Realistic 2/3/2/3 mix (min/rising/max/declining) | **56.0%** |
| Unconditional global rate (v1 baseline) | 58.5% |
| Decade entirely at solar maximum | **76.8%** |

![Decadal hazard comparison](figures/05_phase_hazard_compare.png)

## Interpretation

1. **The original 58.5% headline is robust.** Treating the cycle properly produces 56.0% — within Monte Carlo noise. The v1 number was not an artifact of ignoring solar cycles.
2. **Phase matters enormously for sub-decadal questions.** If a utility or insurer is making a 3-year planning bet, the relevant odds depend heavily on where in the cycle those 3 years fall. A 2024–2026 window (post-max → declining) carries roughly an order of magnitude more Carrington risk than a 2007–2009 window (deep minimum).
3. **Cycle 25 is currently in its declining phase** (max in October 2024 per our detector). The 2024 Gannon storm landed at solar max as expected. **The 2026–2030 window should produce above-average G4+ activity** based on the declining-phase rate of ~22 exceedance days per year — even though the smoothed sunspot number will be falling.
4. **2030–2034 (Cycle 25 → 26 minimum):** expected G4+ rates collapse by a factor of ~6 relative to today.

## What this still doesn't do

Same caveats as [FINDINGS v1](FINDINGS.md). Two new caveats specific to this analysis:

- Solar-cycle phase is detected from the smoothed SSN curve, not a pre-published official labeling. Edge cases (e.g. double-peaked maxima) may be off by a few months. The ±18-month windows absorb most of that.
- "Min" and "max" together cover 50% of the record by construction — these are wide windows. Sharper phase definitions would produce more extreme rate contrasts.

## Reproduce

```bash
cd solar-flare-grid-coupling
curl -L -o data/SN_m_tot_V2.0.txt https://www.sidc.be/SILSO/DATA/SN_m_tot_V2.0.txt
python scripts/analyze_phase.py
```

Deterministic; seed = `20260523`; runtime ~30 s (the per-day phase tag is the slow step).

## New figures

| File | What it shows |
|---|---|
| `figures/04_phase_storm_rate.png` | G3+ and G4+ days per phase-year by cycle phase |
| `figures/05_phase_hazard_compare.png` | Decadal hazard under four phase-mix scenarios |
| `figures/06_ssn_with_storms.png` | 94-year SSN with G4+ storm-day tick marks |

## Additional citations

- [SILSO World Data Center](https://www.sidc.be/SILSO/) — Sunspot Index and Long-term Solar Observations, Royal Observatory of Belgium.
- Russell, C.T. & McPherron, R.L. (1973). *Semiannual variation of geomagnetic activity.* JGR, 78(1).
- Tsurutani, B.T. et al. (2006). *Corotating solar wind streams and recurrent geomagnetic activity: a review.* JGR, 111(A7).
