# `universality/` — Hawkes universality across solar, hurricanes, earthquakes

Cross-domain comparison branch off the main solar Hawkes work. **Not** loaded by default from `main` — this lives on the `universality` branch.

See **[FINDINGS_universality.md](FINDINGS_universality.md)** for the 3-domain writeup with analogy, math, results, and caveats.

For the broader 24-row peer-reviewed cross-domain table, see **[periodic-table/FINDINGS_periodic_table.md](periodic-table/FINDINGS_periodic_table.md)**.

## Contents

```
universality/
├── FINDINGS_universality.md      ~250-line writeup (start here)
├── README.md                      this file
├── scripts/
│   ├── build_universality_table.py   main analysis (3 domains, 2 figures, 2 JSONs)
│   └── scaling_test.py               OLS log-log slope test
├── data/
│   ├── universality_summary.json     n, t_half, kappa per domain with CIs
│   ├── universality_samples.npz      raw posterior arrays
│   └── scaling_test.json             OLS slope = 0.24, R² = 0.73
├── figures/
│   ├── 01_universality_three_panel.png   branching ratio, half-life, κ
│   └── 02_universality_scaling.png       t_half vs τ_forcing on log-log
└── periodic-table/                       24-row peer-reviewed Hawkes table (v1)
    ├── FINDINGS_periodic_table.md        ~190-line writeup (popcorn-machine framing)
    ├── SCHEMA.md                         25-column CSV schema
    ├── data/periodic_table_v1.csv        24 curated rows, 7 domains
    ├── papers/                           candidate paper queues
    ├── scripts/                          builder + figure scripts (seed 20260523)
    └── figures/                          4 PNGs (scatter, n histograms, κ signs, t_half range)
```

## Quick numbers

| Domain | n (branching ratio) | t_half | κ (mark) |
|---|---|---|---|
| Solar G4+ (v15)        | 0.134 [0.078, 0.241]    | 1.16 d    | +1.06 / Kp |
| Hurricanes Cat-3+ (v3) | 0.006 [0.002, 0.012]    | 17.4 d    | +0.17 / cat |
| ETAS quakes (lit.)     | 0.5 – 0.8               | 3 min – 7 h | implicit |

All three: sub-critical (n < 1), monotonically-decaying kernel, mark-positive productivity. The *structural* form universalizes. The *dimensionless* ratio t_half / τ_forcing does **not** (OLS slope = 0.24, far from the slope=1 scale-invariance reference).

## Inputs not committed here

- `solar_v15_idata.pkl` (10 MB PyMC trace) — produced by `solar_flare_hawkes_v15.py` on `main`. Place at `universality/data/solar_v15_idata.pkl` to re-run.
- `hurricane_v3_*` (bootstrap params + event CSV) — lives in a separate local-only repo, not on GitHub. Summary JSON (`hurricane_v3_summary.json`) is sufficient for the headline numbers but not for the per-bootstrap κ histogram in panel C.

## License

MIT (code) · CC0 (data + figures).
