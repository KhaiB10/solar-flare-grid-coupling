# `periodic_table.csv` schema

One row = one published Hawkes fit (one dataset, one model spec, one set of reported parameters).

## Columns

| col | type | required | notes |
|---|---|---|---|
| `row_id` | int | yes | running counter, unique |
| `domain` | str | yes | one of: `seismology`, `heliophysics`, `severe_weather`, `social`, `finance`, `crime`, `neuroscience`, `epidemiology`, `infrastructure`, `ecology` |
| `subdomain` | str | yes | freeform, e.g. "Twitter retweet cascade — 2014 World Cup" |
| `n_branching` | float | yes | branching ratio, dimensionless. If only reported as range, store midpoint. |
| `n_branching_lo` | float | no | lower bound if CI/range reported |
| `n_branching_hi` | float | no | upper bound if CI/range reported |
| `t_half_seconds` | float | yes | kernel half-life converted to seconds. For exp kernel: ln(2)/β. For power-law (Omori): c·(2^(1/p)−1). Must specify natural unit in `notes`. |
| `t_half_seconds_lo` | float | no | lower bound |
| `t_half_seconds_hi` | float | no | upper bound |
| `kappa_sign` | str | yes | one of: `+`, `-`, `0`, `unknown`. The sign of mark productivity. `0` = mark-insensitive (no marks in model). |
| `kappa_value` | float | no | natural-log productivity exponent if reported in a comparable form |
| `kernel_form` | str | yes | one of: `exp`, `power_law`, `omori`, `rayleigh`, `weibull`, `nonparametric`, `other` |
| `n_events` | int | yes | catalog size used for the fit |
| `obs_window_days` | float | no | observation window in days |
| `mu_background` | float | no | background intensity (per day) if reported |
| `paper_title` | str | yes | full paper title |
| `authors_short` | str | yes | "First Author et al." or "Author1 & Author2" |
| `year` | int | yes | publication year |
| `venue` | str | yes | journal / conference name |
| `peer_reviewed` | bool | yes | TRUE only if journal article or top-tier peer-reviewed conference (NeurIPS, ICML, KDD, etc.). FALSE for arXiv preprints with no published venue. |
| `doi_or_url` | str | yes | DOI link preferred, stable URL otherwise |
| `notes` | str | no | mark units, stratification, anything weird ("Mmin=4.0", "Cat-3 threshold", "scaled to 1/day") |
| `extracted_by` | str | yes | `manual`, `wide_research`, `direct_fetch` |
| `extraction_date` | str | yes | ISO date string |

## Conversion conventions

- **Time units → seconds.** This is the universal axis. Days → ×86400. Hours → ×3600. Minutes → ×60.
- **Branching ratio.** Always dimensionless. If a paper reports α and β separately, n = α/β only for univariate exponential Hawkes — for marked or multivariate, compute n = α·E[g(m)] if mark dist available, else leave blank and add a note.
- **κ (mark productivity).** Convert to natural-log scale. If paper reports `productivity ∝ 10^(α·M)`, then `kappa_value = α · ln(10) ≈ α · 2.3`. If paper reports `exp(κ·m)`, kappa_value = κ directly. If mark units differ from natural mark choice for that domain (e.g. paper used log10(intensity) instead of category), flag in notes.
- **Half-life.** Always preferred over decay rate for cross-domain comparability. If only β reported, t_half = ln(2)/β. If only Omori (p, c) reported, t_half = c·(2^(1/p) − 1).

## Quality gate

A row is admitted iff:
1. The paper is peer-reviewed (journal or top-tier conference).
2. Both `n_branching` and `t_half_seconds` are extractable from the paper (or computable from reported α, β, c, p).
3. `kernel_form` is identifiable.

Rows missing any of these go to `papers_excluded.csv` with a reason.
