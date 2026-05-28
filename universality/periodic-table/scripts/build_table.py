"""
Build the curated Periodic Table of Self-Exciting Systems.

Sources:
1. Our own v15/v16 fits (solar, hurricane, ETAS reference) — 3 rows
2. wide_research extraction from 46 peer-reviewed papers — filter to rows
   with at least branching ratio OR kernel form recovered
3. Hand-curated annotations from open-access abstracts for headline numbers
   that the LLM extractor missed but are visible in published abstracts

Output: periodic_table.csv
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw_extraction.csv"
OUT = ROOT / "data" / "periodic_table.csv"

# ----- Our own three rows (from v15/v16/hurricane) -----
OWN_ROWS = [
    {
        "row_id": "OWN-001",
        "domain": "heliophysics",
        "subdomain": "Solar G4+ geomagnetic storms",
        "n_branching": "0.134",
        "n_branching_lo": "0.078",
        "n_branching_hi": "0.241",
        "t_half_seconds": "100224",  # 1.16 d
        "t_half_raw": "1.16 d",
        "kernel_form": "Omori-Utsu power-law",
        "kappa_sign": "+",
        "kappa_value": "+1.06 per Kp index",
        "n_events": "1,156",
        "obs_window_days": "~22,000",
        "mu_background": "see v15",
        "peer_reviewed": "self (defensive publication)",
        "paper_title": "v15 hierarchical Bayesian Hawkes per solar cycle",
        "authors_short": "Diatom Sky R&D",
        "year": "2026",
        "venue": "GitHub: solar-flare-grid-coupling",
        "doi_or_url": "https://github.com/KhaiB10/solar-flare-grid-coupling",
        "notes": "Per-cycle hierarchical fit; SC25 forecast through 2030 published.",
    },
    {
        "row_id": "OWN-002",
        "domain": "tropical_cyclones",
        "subdomain": "North Atlantic Cat-3+ hurricanes",
        "n_branching": "0.006",
        "n_branching_lo": "0.002",
        "n_branching_hi": "0.012",
        "t_half_seconds": "1503360",  # 17.4 d
        "t_half_raw": "17.4 d",
        "kernel_form": "Omori-Utsu power-law",
        "kappa_sign": "+",
        "kappa_value": "+0.17 per Saffir-Simpson category",
        "n_events": "263",
        "obs_window_days": "~62,000",
        "mu_background": "see hurricane v3",
        "peer_reviewed": "self (defensive publication)",
        "paper_title": "Hurricane Hawkes clustering v3",
        "authors_short": "Diatom Sky R&D",
        "year": "2026",
        "venue": "GitHub: hurricane-hawkes-clustering (local)",
        "doi_or_url": "local repo",
        "notes": "Saffir-Simpson category as mark; column-pick bug caught pre-publication.",
    },
    {
        "row_id": "OWN-003",
        "domain": "seismology",
        "subdomain": "Reference ETAS literature range",
        "n_branching": "0.65",
        "n_branching_lo": "0.50",
        "n_branching_hi": "0.80",
        "t_half_seconds": "900",  # ~15 min, mid of 3min-7hr range
        "t_half_raw": "3 min to 7 h",
        "kernel_form": "Omori-Utsu power-law",
        "kappa_sign": "+",
        "kappa_value": "alpha ~ 1.5-2.3 (magnitude productivity)",
        "n_events": "varies",
        "obs_window_days": "varies",
        "mu_background": "varies",
        "peer_reviewed": "yes (literature consensus)",
        "paper_title": "ETAS branching ratio summary (consensus from multiple regional fits)",
        "authors_short": "Ogata, Zhuang, Helmstetter (compilation)",
        "year": "2000-2024",
        "venue": "multiple",
        "doi_or_url": "see seismology rows below",
        "notes": "Headline range from California, Japan, Italy ETAS fits.",
    },
]


def slugify_kernel(k):
    k = (k or "").lower()
    if "expon" in k:
        return "exponential"
    if "omori" in k or "power" in k:
        return "Omori-Utsu power-law"
    if "mittag" in k:
        return "Mittag-Leffler"
    if "nonparam" in k:
        return "nonparametric"
    if "param" in k:
        return "parametric"
    if k in ("na", "", "other"):
        return "NA"
    return k


def main():
    # Read raw_extraction
    with open(RAW) as f:
        rows = list(csv.DictReader(f))

    cols = [
        "row_id", "domain", "subdomain", "n_branching",
        "n_branching_lo", "n_branching_hi", "t_half_seconds", "t_half_raw",
        "kernel_form", "kappa_sign", "kappa_value", "n_events",
        "obs_window_days", "mu_background", "peer_reviewed",
        "paper_title", "authors_short", "year", "venue", "doi_or_url",
        "notes",
    ]

    out_rows = list(OWN_ROWS)
    counter = 1

    for r in rows:
        n_str = r["Branching ratio n (point estimate)"]
        kernel = r["Kernel form (exp/Omori/power-law/other)"]
        # Keep if we have ANY of: n, kernel (non-other), mark
        has_n = n_str not in ("NA", "", "N/A", "n/a")
        has_kernel = kernel not in ("NA", "", "N/A", "n/a", "other")
        mark = r["Mark/covariate magnitude"]
        has_mark = mark not in ("NA", "", "N/A", "n/a")
        if not (has_n or has_kernel or has_mark):
            continue

        # Convert raw t_half to seconds if possible
        traw = r["Half-life as reported (raw units)"]
        tsec = r["Half-life in seconds (converted)"]
        if tsec in ("NA", "", "N/A"):
            tsec = "NA"

        row = {
            "row_id": f"LIT-{counter:03d}",
            "domain": r["Domain"],
            "subdomain": r["Subdomain / Dataset"][:120],
            "n_branching": n_str if has_n else "NA",
            "n_branching_lo": "NA",
            "n_branching_hi": "NA",
            "t_half_seconds": tsec,
            "t_half_raw": traw,
            "kernel_form": slugify_kernel(kernel),
            "kappa_sign": r["Mark/covariate effect sign (+/-/0/NA)"],
            "kappa_value": mark,
            "n_events": r["Sample size (#events)"][:50],
            "obs_window_days": r["Observation window"][:50],
            "mu_background": r["Background rate (with units)"][:80],
            "peer_reviewed": r["Peer Reviewed (yes/no/preprint)"],
            "paper_title": r["Paper Title"][:120],
            "authors_short": r["Authors (short)"],
            "year": r["Year"],
            "venue": r["Venue"][:60],
            "doi_or_url": r["DOI"],
            "notes": r["Notes / Caveats"][:200],
        }
        out_rows.append(row)
        counter += 1

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in out_rows:
            w.writerow({k: r.get(k, "") for k in cols})

    print(f"Wrote {len(out_rows)} rows to {OUT}")

    # Summary
    from collections import Counter
    domain_counts = Counter(r["domain"] for r in out_rows)
    print("\nDomain breakdown:")
    for d, c in sorted(domain_counts.items(), key=lambda x: -x[1]):
        print(f"  {d:25s} {c}")
    n_valid = sum(1 for r in out_rows if r["n_branching"] not in ("NA", ""))
    t_valid = sum(1 for r in out_rows if r["t_half_seconds"] not in ("NA", ""))
    print(f"\nRows with n_branching: {n_valid}")
    print(f"Rows with t_half:      {t_valid}")


if __name__ == "__main__":
    main()
