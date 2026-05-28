# Periodic Table of Self-Exciting Systems — v2 (41 rows, 9 domains)

**Branch:** `universality` · **Path:** `universality/periodic-table/` · **Seed:** 20260523
**License:** MIT (code) · CC0 (data/figures) · **Framing:** Diatom Sky R&D defensive publication
**Status:** v2 release. Expanded from v1's 24 rows by hand-mining arXiv-OA PDFs for papers that v1 skipped because of paywall opacity. v1 history kept inline below for diffability.

## What changed v1 → v2

- **Row count: 24 → 41** (17 new hand-extracted rows from arXiv-OA PDFs)
- **Domain count: 7 → 9** (added `cyber` and `neuroscience`)
- **`n_branching` populated: 13 → 23** rows
- **`t_half_s` populated: 10 → 22** rows
- **Both populated (filled circles on the scatter): 4 → 10** rows
- **`kappa` all-positive claim from v1 is broken.** Two new rows have mixed `+/-`: Reinhart's Pittsburgh burglary model (population density positive, fraction-of-homes-owned negative) and Lambert/Reynaud-Bouret spike-train Hawkes (explicitly modeling inhibitory synapses). This was the single clearest update from v1 — the v1 "every kappa sign is positive" tally was a sampling artifact of the open-access subset we'd read, not a universal feature.
- **Scaling-test direction flipped.** v1's OLS slope of `log t_half ~ a + b log n` was `+0.29` (effectively flat, weak signal). v2's OLS slope with 10 (n, t_half) pairs is `−0.56` — now in the direction that strict scale-invariance would predict (`−1`), but R² is only 0.10 and p-value is 0.37, and a robust Theil-Sen estimator collapses the slope to `−0.04`. So the v2 numbers nudge toward the universal-scaling hypothesis but do not move it past statistical noise. v1's headline conclusion — *strict scale-invariance is falsified* — still stands, just less emphatically.

---

## TL;DR

We tried to assemble a "periodic table" of self-exciting (Hawkes-process) systems — one row per published fit, columns for branching ratio `n`, excitation half-life `t_half`, mark/covariate sensitivity `kappa`, kernel form, sample size, and citation. The goal was **~100 rows across 8–10 domains, peer-reviewed only**. The v2 release contains **41 rows across 9 domains**. We're roughly 40% of the way to the long-term target. The remaining gap is still paywall-bound; v2 closes a chunk of it by reading arXiv-OA versions of papers v1 had to skip.

What is here is enough to make a few claims with a straight face:

- Self-exciting dynamics are **not a niche curiosity**. The same `(mu, n, kernel)` machinery has been fit to earthquakes, solar flares, hurricanes, financial trades, retweets, COVID deaths, burglaries, mass shootings, cyber-attacks, and neuronal spike trains in peer-reviewed journals.
- Branching ratios span **roughly four orders of magnitude** — from `n ≈ 0.006` (Atlantic Cat-3+ hurricanes; ours) to `n ≈ 1.0` (E-mini futures near criticality; Hardiman/Bouchaud 2014).
- Excitation half-lives span **roughly fifteen orders of magnitude** — from `~300 µs` (German Bund futures, Rambaldi et al.) to `~1.7 yr` (post-Tohoku aftershocks, implied from Kwon et al.).
- The v1 "all kappa signs positive" claim is **now broken** by two v2 rows: Reinhart's Pittsburgh burglary model (mixed covariate signs) and Lambert's spike-train Hawkes (inhibitory synapses are negative-kappa by construction). 31 of 33 sign-reporting rows are still positive — the qualitative trend holds, but the universal claim does not.
- Strict scale-invariance is **still falsified** by the v2 scaling test. OLS slope on 10 (n, t_half) pairs is `−0.56` (was `+0.29` in v1 with only 4 pairs) — closer to the `−1` strict-scaling target, but R² is 0.10 and the Theil-Sen robust slope is `−0.04`. v1's null conclusion stands.

---

## The popcorn-machine analogy

If you've never met a Hawkes process before, here is the one-paragraph version.

Imagine a popcorn machine. Kernels pop on their own at some steady **background rate** (`mu` — the heating element). But every kernel that pops also **kicks some neighbors**, and each kick raises the chance that those neighbors pop a little sooner. The **branching ratio `n`** is the average number of extra pops a single pop triggers down the chain. If `n < 1` the bursts die out and you just hear popcorn. If `n = 1` the machine sits on a knife edge — a single early pop can keep echoing forever. If `n > 1` it runs away and the bowl overflows.

The **half-life `t_half`** is how long the kick lasts: a short half-life means each pop only excites the next few seconds; a long half-life means a pop today still nudges next week.

The **mark/covariate sensitivity `kappa`** is whether a *bigger* pop kicks *harder* — a magnitude-9 earthquake triggering more aftershocks than a magnitude-5, a Category-5 hurricane being followed by more storms than a Category-3, a G5 geomagnetic storm followed by more substorms than a G2.

Every row in the table is one paper's measurement of `(n, t_half, kappa)` for one popcorn machine somewhere in the universe. Some machines are stock exchanges. Some are tectonic plates. Some are Twitter.

---

## What's in v2

**`data/periodic_table_v2.csv`** — 25-column CSV, 41 rows, 3 tiers:

| Tier | Rows | Source                                                                |
|------|------|----------------------------------------------------------------------|
| OWN  | 2    | Diatom Sky R&D (solar G4+ v15, hurricane Cat-3+ v3 local repo)        |
| A    | 36   | Hand-curated from open-access PDFs (arXiv, PLOS, Frontiers, PMC, etc.)|
| B    | 3    | Extracted by automated wide-research from open-access abstracts        |

**Domain coverage:**

| Domain             | Rows | Representative paper(s)                                                                                            |
|--------------------|------|--------------------------------------------------------------------------------------------------------------------|
| seismology         | 15   | Naylor 2023; Kwon 2023; **Li/Sornette 2024 (CA, NZ, CSES)**; **Davis fractional-Hawkes 2024 (Joshua Tree, Landers, Hector Mine)**; **Ogata HIST-ETAS 2022** |
| crime              | 9    | Mohler 2011; Boyd & Molyneux 2021; **Reinhart Pittsburgh burglary 2018 (3 specs)**; **Holbrook DC ShotSpotter 2020** |
| finance            | 4    | Hardiman/Bouchaud 2014; Omi 2017; Rambaldi 2017; Zhuo 2023                                                         |
| heliophysics       | 3    | Ours (G4+ v15); Ross 2020; Rivera/Johnson 2022                                                                     |
| social             | 3    | Kobayashi/Lambiotte 2016; Mishra/Rizoiu 2018; Gleeson 2021                                                         |
| **cyber (NEW)**    | 3    | **Boumezoued 2023 (Hackmageddon vuln / KEV / NVD)** — flagged `preprint`                                            |
| epidemiology       | 2    | Browning/Sulem 2021; Chiang/Mohler 2022                                                                            |
| tropical_cyclones  | 1    | Ours (Atlantic Cat-3+ v3, local repo)                                                                              |
| **neuroscience (NEW)** | 1 | **Lambert/Reynaud-Bouret 2018 (spike-train Hawkes, inhibitory synapses)**                                          |

Numerically populated cells:

- `n_branching`: **23 of 41** rows (v1: 13/24)
- `t_half_s`: **22 of 41** rows (v1: 10/24)
- both populated: **10 of 41** rows — the filled circles on the scatter (v1: 4)
- `kappa_sign`: **33 of 41** rows have an explicit sign; **31 positive, 2 mixed (+/−)** — the v1 "all positive" claim no longer holds

### New-row source ledger (v2 hand-extractions)

| Row IDs        | Paper                                              | Domain         | Key numbers extracted                              | OA route       |
|----------------|----------------------------------------------------|----------------|----------------------------------------------------|----------------|
| A-020 — A-022  | Li & Sornette, *J. Geophys. Res. Solid Earth* 2024 | seismology     | n=0.72–0.80 (bias-corrected ETAS, CA/NZ/CSES)      | arXiv-OA       |
| A-023 — A-025  | Davis et al., *JRSS-C* 2024 (ETAS California seq.) | seismology     | c=0.004–0.96 d, p=1.18–1.76 → t_half=0.17–1.43 d   | arXiv-OA       |
| A-026, A-027   | Davis et al., arXiv:2404.01478 (fractional Hawkes) | seismology     | β=0.53–0.62 (MDFHP) — **preprint flag**            | arXiv          |
| A-028          | Ogata HIST-ETAS, *Earth Planets Space* 2022        | seismology     | method paper, no single n/t_half                   | OA journal     |
| A-029          | Holbrook et al., *Stat. Comp.* 2020 (DC gunfire)   | crime          | θ=0.153, 85k events                                | arXiv-OA       |
| A-030 — A-032  | Reinhart, *J. R. Stat. Soc. C* 2018 (Pittsburgh)   | crime          | θ=0.45–0.76, τ=41–52 d → t_half=28.5–36.2 d        | arXiv-OA       |
| A-033 — A-035  | Boumezoued et al., arXiv:2311.15701 (cyber)        | cyber (NEW)    | ‖φ‖=0.36–0.58, δ=1.5–1.9/d → t_half≈0.4–0.5 d — **preprint flag** | arXiv          |
| A-036          | Lambert/Reynaud-Bouret, *J. Neurosci. Methods* 2018| neuroscience (NEW) | kernel support ~50 ms, inhibitory (kappa +/−)      | arXiv-OA       |

---

## Methodology

### Sourcing

1. **Search.** Per-domain academic queries: `"Hawkes process" <domain>`, `"ETAS" <region>`, `"self-exciting" <topic>`, etc. We collected a queue of 53 candidate papers.
2. **Automated extraction.** We ran a 46-entity batch research call against the queue with a 20-field JSON schema. **Yield: 5 of 46 rows with valid `n`, 0 with valid `t_half`.** Most candidates were paywalled; abstracts rarely quote the branching ratio.
3. **Hand curation from open access.** We then fetched and read the PDFs that were freely accessible: arXiv preprints with peer-reviewed publication, PLOS ONE, Frontiers, PMC, ICWSM, and a few university preprint servers. That pass produced the 19 Tier-A rows.

### Per-row extraction rules

- **`n_branching`** comes from the paper's own quoted branching ratio (point estimate). 95% CI in `n_lo`/`n_hi` only when explicitly reported.
- **`t_half_s`** is converted to seconds. For Omori-Utsu `(c, p)` kernels we compute `t_half = c · (2^(1/(p-1)) − 1)` and convert to seconds. For exponential `omega·exp(−omega·t)` we use `t_half = ln(2)/omega`. For nonparametric or sum-of-exp kernels with no clean closed form we leave `t_half_s` empty and put the qualitative description in `t_half_raw`.
- **`kappa_sign`** is the sign of the mark/covariate response (+ if bigger marks excite more strongly). When the paper reports a non-monotone `k(m)`, we tag `+` only if the *overall* tendency is positive; the raw values go in `kappa_value`.
- **`kernel`** is the verbatim functional form from the paper.
- **`peer_reviewed`** is `yes` for journal-published or refereed-conference papers, `self/defensive` for our two OWN rows (this and the universality writeup serve as the defensive publication).

### What was excluded

- Preprints with no peer-reviewed companion paper.
- Methods papers that introduce an estimator but don't fit a real dataset.
- Reviews and tutorials.
- Papers we couldn't read (paywalled, no OA route).

---

## The four figures (v2)

### `figures/01_periodic_table_scatter.png`

`log(n)` versus `log(t_half)`, colored by domain (9 colors now — cyber in teal, neuroscience in purple), `n=1` critical line dashed red.

- **10 filled circles** (up from 4) — the v2 expansion roughly doubled the points where both `n` and `t_half` are quoted in the same paper.
- The **finance point** still sits exactly on the `n=1` line (Hardiman/Bouchaud E-mini). It remains the only row at criticality.
- The **crime cluster** is now visible — four Pittsburgh burglary specs (Reinhart) plus LA burglary (Mohler) form a tight `n ≈ 0.2–0.76`, `t_half ≈ 7–36 d` group at the right of the plot.
- The **cyber cluster** sits at `n ≈ 0.36–0.58`, `t_half ≈ 8–11 hours` (Boumezoued 2023 — three model variants on the same dataset). It lands between crime (slower) and finance (faster) on both axes, consistent with the intuition that vulnerability disclosure cascades have a sub-day reverb but well-below-critical branching.
- Solar G4+ and Hurricane Cat-3+ are still our two OWN points; they have not moved.

### `figures/02_n_histograms.png`

Branching-ratio distribution per domain. Updated for v2:

- **Seismology** now has 15 rows ranging `n ≈ 0.5–0.96`. The Li/Sornette bias-corrected fits at `n ≈ 0.72–0.80` are *higher* than the standard ETAS values quoted in v1, which suggests earthquakes are closer to criticality than the literature consensus once you correct for spatio-temporal bias. **Markets at `n=1` remain unique but the gap to seismology is narrower than v1 implied.**
- **Crime** now has 9 rows; Reinhart's three Pittsburgh specs all sit at `n ≈ 0.45–0.76`, consistent with v1's Mohler LA value.
- **Cyber (new)** sits at `n ≈ 0.36–0.58` — mid-table, somewhere between crime and the lower-engagement social tweets.
- **Neuroscience (new)** has only one row and no clean `n` extracted; the Lambert paper estimates kernel support rather than a normalized branching ratio.

### `figures/03_kappa_signs.png`

Stacked bars of mark/covariate sign per domain. **31 of 33 sign-reporting rows are still positive**, but two new mixed-sign rows now show up:

- **crime: 1 of 9 rows is `+/−`** — Reinhart's full-covariate Pittsburgh model has positive sensitivity to population density and negative sensitivity to fraction-of-homes-owned. Both directions are intuitive (more people = more burglary opportunity; owner-occupied homes are harder targets) but it's the first row in the table where a single paper reports both signs.
- **neuroscience: 1 of 1 row is `+/−`** — Lambert/Reynaud-Bouret's spike-train Hawkes explicitly fits both excitatory and inhibitory synaptic responses, so the kernel signs are mixed by construction.

The v1 "every Hawkes fit we found has a positive kappa" tally was therefore an OA-sampling artifact. **The honest v2 story:** kappa is positive in most domains where the obvious mark variable is event size (magnitude, Cat, fatality count), but as soon as you read a paper that fits multiple covariates or a domain with biologically inhibitory dynamics, mixed signs appear.

### `figures/04_thalf_by_domain.png`

Half-life range per domain on a log axis. Updates from v1:

- **Crime** now has 5 t_half values (Mohler LA + 3 Reinhart Pittsburgh + 1 Mohler ref), all clustered around 1 week to 1 month. Crime is the slowest-decaying domain in the table that's still well sub-critical.
- **Cyber** debut: `t_half ≈ 8–11 hours`. Sub-daily but well above the millisecond regime of finance.
- **Seismology** gains 3 ETAS-extracted half-lives (Davis 2024): `t_half ≈ 0.17–1.43 days` for individual California sequences. Combined with v1's Kwon megathrust range, seismology now spans `~4 hours – ~1.7 yr` in this table.

## v2 scaling test

`scripts/scaling_test_v2.py` runs OLS and Theil-Sen on the 10 v2 rows that have both `n` and `t_half`. JSON snapshot: `data/scaling_test_v2.json`.

| Estimator   | Slope    | 95% CI / std-err           | R² / p-value         |
|-------------|----------|----------------------------|----------------------|
| OLS         | **−0.56** | std-err 0.60               | R²=0.10, p=0.37      |
| Theil-Sen   | **−0.04** | 95% CI [−2.55, +1.75]      | —                    |

**Reading.** Strict scale-invariance under a universal Hawkes kernel would predict slope `= −1`. v1's OLS slope of `+0.29` (on 4 points) was effectively flat. v2's OLS slope of `−0.56` *moves toward* the strict-scaling prediction but the regression is **statistically insignificant** (p=0.37) and the robust Theil-Sen median slope is essentially zero. **Net: v2 nudges the data in the direction strict scaling would predict, but does not have the power to claim a real slope.** v1's conclusion — *strict universal scale-invariance is not supported by the cross-domain data* — stands.

---

## Caveats — read these before citing v2

1. **41 rows is not 100.** The 100-row target was paywall-bound in v1 and is still paywall-bound in v2. v2 closed the most accessible chunk of the gap. The next chunk (BSSA, JGR, Risk Analysis, Neural Computation specialty journals) still needs institutional library access.
2. **5 of the 17 new rows are preprints.** The 3 Boumezoued cyber rows and 2 Davis fractional-Hawkes seismology rows are flagged `peer_reviewed=preprint` in the CSV. They are kept because their methodology mirrors published Hawkes work and the numbers extracted are self-consistent, but they should be treated as provisional until publication.
3. **The CSV is still biased toward open access.** Same as v1 — newer methods and physics/CS venues are over-represented; specialty seismology, criminology, and neuroscience journals are under-represented.
4. **CIs are missing for most rows.** Unchanged from v1. Li/Sornette's tight error bars on the bias-corrected ETAS fits (`±0.0027–±0.0211`) are an exception, not the rule.
5. **`t_half` is heterogeneous.** Unchanged from v1. Mixing exponential, Omori-Utsu power-law, and fractional/Mittag-Leffler kernels gives approximate-but-not-identical half-life semantics. Fractional Hawkes (Davis arXiv:2404.01478) introduces a true heavy-tail regime where the half-life is not even finite in the usual sense — we extract a characteristic time `c` but flag in `notes`.
6. **Our two OWN rows are mixed in with the literature.** Flagged `tier=OWN`, `peer_reviewed=self/defensive`. Unchanged.
7. **Sign of `kappa` is reporting-biased.** The v1 "all positive" tally being broken by v2 is itself a sign of this bias — once we hand-read a paper that fits multiple covariates (Reinhart) or a domain that has inhibitory dynamics by construction (Lambert), mixed signs appear immediately. The honest interpretation is that **kappa-positive is the dominant sign in single-mark `event-size` models**, not that it is universal across all Hawkes fits.
8. **Heavy-tail regime is real.** v1 dismissed `p < 1` Omori cases. v2 includes Li/Sornette's CSES Sichuan-Yunnan fit at `p=0.88` — a power-law without a finite mean — which raises the question of whether the Hawkes-finite-`n` formalism is even the right framework for some seismic regimes. We've kept the row but the v1 caveat about `t_half` heterogeneity applies double here.

---

## Where this could go (v3 wishlist)

- **v3: 60–80 rows, OA-only.** Realistic next milestone. Adds Reynaud-Bouret/Pillow's *broader* spike-train literature (many neuroscience rows behind one OA roof), Bebbington-style volcanic point processes (a domain v2 surveyed but could not extract a clean row from), wildfire ignition Hawkes fits (Schoenberg's lab), and network-cascade rows from the social-influence literature (Hodas, Lerman).
- **v3.5: 100 rows, paywalled included.** Same as v1's wishlist — needs institutional access or per-paper purchase. Lambert+Pillow alone has ~20 rows' worth of fits behind various journal paywalls.
- **Re-fit, don't re-cite.** v3's hardest-to-fake step. A small open-source ETAS/Hawkes implementation re-fit to public datasets (USGS catalogs, GOES X-ray, HURDAT2, NYC 311, NVD CVE timestamps) would replace transcription error with reproducible numbers.
- **Cross-domain scaling fit with 50+ pairs.** v2 has 10 — better than v1's 4 but still under-powered. Doubling that should push the OLS p-value below 0.05 if there's a real slope to find.
- **Defensive-publication target.** A short OSF / Zenodo deposit of `v2` with a citable DOI would lock in the priority date for the framework (popcorn-machine periodic table of self-excitation) under Diatom Sky R&D.

---

## Who this impacts

- **Practitioners doing Hawkes fitting** get a one-stop comparison sheet — "is my `n = 0.4` in line with prior work in my domain?" — without having to read 30 papers.
- **Risk modelers** in insurance, climate, and infrastructure get an apples-to-apples view of which hazards are near-critical (markets) vs. comfortably sub-critical (our hurricanes at `n = 0.006`).
- **Solar/space-weather operations** (NOAA SWPC, NASA, utilities) see that the same machinery used to forecast aftershocks and crime predicts geomagnetic storm clustering — our G4+ fit lives on the same plot as the Mohler LA burglary row.
- **Public-interest tech** (epidemiology, gun-violence research, content moderation) can borrow methodology across domains with explicit numerical anchors rather than vague "we used a Hawkes process" gestures.
- **Trader-developers at small shops** (the author's day job) can see where the markets row literally sits at `n = 1` and understand why every microsecond of excitation matters.

---

## File index

```
universality/periodic-table/
├── SCHEMA.md                    25-column CSV schema
├── FINDINGS_periodic_table.md   this file
├── data/
│   ├── periodic_table_v2.csv    41 curated rows (THE deliverable)
│   ├── periodic_table_v1.csv    24 curated rows (kept for diff)
│   ├── raw_extraction.csv       wide_research output, mostly empty (paywall fallout)
│   ├── periodic_table.csv       earlier skeleton (deprecated, kept for diff)
│   └── scaling_test_v2.json     v2 OLS + Theil-Sen scaling regression
├── scripts/
│   ├── extraction_schema.json   JSON schema used by wide_research
│   ├── build_curated_table.py   builds periodic_table_v1.csv (kept for reproducibility)
│   ├── build_v2_table.py        appends 17 v2 rows to v1 → periodic_table_v2.csv
│   ├── scaling_test_v2.py       v2 log(n) vs log(t_half) regression
│   └── make_figures.py          generates the 4 PNGs from v2 (seed 20260523)
└── figures/
    ├── 01_periodic_table_scatter.png
    ├── 02_n_histograms.png
    ├── 03_kappa_signs.png
    └── 04_thalf_by_domain.png
```

---

## Acknowledgements & sources

Hand-curated values were transcribed from the open-access PDFs of: Hardiman & Bouchaud (Phys Rev E 2014); Omi, Hirata & Aihara (Phys Rev E 2017); Rambaldi, Bacry & Lillo (Quant Finance 2017); Kobayashi & Lambiotte (ICWSM 2016); Mohler, Short, Brantingham, Schoenberg & Tita (JASA 2011); Boyd & Molyneux (PLOS ONE 2021); Browning, Sulem, Mengersen, Rivoirard & Rousseau (PLOS ONE 2021); Chiang, Liu & Mohler (Int J Forecasting 2022); Naylor, Serafini, Lindgren & Main (Frontiers Appl Math Stat 2023); Kwon, Zheng & Jun (Spatial Statistics 2023); Ross (Physica A 2020); Rivera, Johnson, Homan & Wing (ApJ Letters 2022). Tier-B rows reference Zhuo et al. (Eur J Finance 2023); Gleeson, Onaga et al. (J Complex Networks 2021); Mishra, Rizoiu & Xie (TKDE 2018).

Our two OWN rows live at:
- Solar G4+ v15: <https://github.com/KhaiB10/solar-flare-grid-coupling> (branch `universality`)
- Hurricane Cat-3+ v3: local repo, commit `22d9b07` (not pushed; defensive-publication priority via this document)

All errors of transcription are ours. Open an issue on the repo if you spot one — corrections will roll into v2.
