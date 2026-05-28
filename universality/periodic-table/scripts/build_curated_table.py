"""
Build the curated Periodic Table of Self-Exciting Systems v1.

Two tiers of data quality:
  TIER A: hand-verified from open-access PDFs (arXiv, PMC, PLOS, journal-OA)
  TIER B: from wide_research extraction with high confidence (n reported)
  OWN:    our own v15/v16 fits

t_half conversion: 1 d = 86400 s, 1 h = 3600 s, 1 min = 60 s.

Random seed: 20260523 (per user instruction)
"""
import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "periodic_table_v1.csv"

DAY = 86400.0
HOUR = 3600.0
MIN_S = 60.0

def half_life_exp(beta_per_unit, unit_seconds):
    """For exp(-beta*t), t_half = ln(2)/beta in seconds."""
    return math.log(2) / beta_per_unit * unit_seconds if beta_per_unit > 0 else None

def half_life_omori(c_days, p):
    """For c*(1+t/c)^-p, t_half_in_days = c*(2^(1/(p-1)) - 1)."""
    if p <= 1:
        return None
    return c_days * (2**(1.0/(p-1)) - 1) * DAY  # convert days to seconds

# ============================================================
# OWN ROWS (our work)
# ============================================================
OWN = [
    {
        "row_id": "OWN-001", "domain": "heliophysics",
        "subdomain": "Solar G4+ geomagnetic storms (SC23-SC25)",
        "n_branching": 0.134, "n_lo": 0.078, "n_hi": 0.241,
        "t_half_s": 1.16 * DAY, "t_half_raw": "1.16 d",
        "kernel": "Omori-Utsu power-law",
        "kappa_sign": "+", "kappa_value": "+1.06 per Kp index",
        "n_events": 1156, "obs_window": "1968-2025 (~57 yr)",
        "mu_background": "hierarchical per-cycle",
        "venue": "GitHub: solar-flare-grid-coupling",
        "authors": "Diatom Sky R&D", "year": 2026,
        "peer_reviewed": "self/defensive",
        "doi_or_url": "https://github.com/KhaiB10/solar-flare-grid-coupling",
        "notes": "v15 hierarchical Bayesian; SC25 forecast through 2030",
        "tier": "OWN",
    },
    {
        "row_id": "OWN-002", "domain": "tropical_cyclones",
        "subdomain": "North Atlantic Cat-3+ hurricanes",
        "n_branching": 0.006, "n_lo": 0.002, "n_hi": 0.012,
        "t_half_s": 17.4 * DAY, "t_half_raw": "17.4 d",
        "kernel": "Omori-Utsu power-law",
        "kappa_sign": "+", "kappa_value": "+0.17 per Saffir-Simpson cat",
        "n_events": 263, "obs_window": "1851-2024 (~173 yr)",
        "mu_background": "see hurricane v3",
        "venue": "GitHub: hurricane-hawkes-clustering (local)",
        "authors": "Diatom Sky R&D", "year": 2026,
        "peer_reviewed": "self/defensive",
        "doi_or_url": "local repo",
        "notes": "Saffir-Simpson cat mark; col-pick bug caught pre-publication",
        "tier": "OWN",
    },
]

# ============================================================
# TIER A (hand-verified from OA PDFs)
# ============================================================
TIER_A = [
    # --- Finance ---
    {
        "row_id": "A-001", "domain": "finance",
        "subdomain": "E-mini S&P 500 futures mid-price (1998-2011)",
        "n_branching": 1.00, "n_lo": 0.95, "n_hi": 1.05,
        "t_half_s": 1000.0, "t_half_raw": "~1000 s (crossover regime)",
        "kernel": "power-law (eps=0.15 short, 0.45 long)",
        "kappa_sign": "NA", "kappa_value": "NA",
        "n_events": "~10^6 per 2-month window", "obs_window": "1998-2011",
        "mu_background": "spikes during crisis",
        "venue": "Phys Rev E", "authors": "Hardiman, Bercot, Bouchaud", "year": 2014,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1103/PhysRevE.90.062807",
        "notes": "Branching ratio fluctuates about n_c=1; markets near criticality",
        "tier": "A",
    },
    {
        "row_id": "A-002", "domain": "finance",
        "subdomain": "Nikkei 225 mini tick data (Jan-Jun 2016)",
        "n_branching": 0.41, "n_lo": None, "n_hi": None,
        "t_half_s": None, "t_half_raw": "NA (sum-of-exp kernel)",
        "kernel": "sum of exponentials",
        "kappa_sign": "NA", "kappa_value": "NA",
        "n_events": "avg 2090/session, range 465-18505", "obs_window": "Jan 4 - Jun 30 2016",
        "mu_background": "time-dependent (Bayesian change-point)",
        "venue": "Phys Rev E", "authors": "Omi, Hirata, Aihara", "year": 2017,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1103/PhysRevE.96.012303",
        "notes": "n=0.41 indicates relative importance of exogenous factors",
        "tier": "A",
    },
    {
        "row_id": "A-003", "domain": "finance",
        "subdomain": "Bund Future EUREX (Jul 2013 - Nov 2014)",
        "n_branching": None, "n_lo": None, "n_hi": None,
        "t_half_s": 3e-4, "t_half_raw": "~300 us reaction time",
        "kernel": "nonparametric multivariate (Wiener-Hopf)",
        "kappa_sign": "+", "kappa_value": "volume-dependent kernels",
        "n_events": "~150e6", "obs_window": "Jul 2013-Nov 2014",
        "mu_background": "varies by component",
        "venue": "Quant Finance", "authors": "Rambaldi, Bacry, Lillo", "year": 2017,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1080/14697688.2016.1260759",
        "notes": "576 kernels over 24-dim model; chars timescale ~300us",
        "tier": "A",
    },
    # --- Social ---
    {
        "row_id": "A-004", "domain": "social",
        "subdomain": "Twitter retweets (TiDeH, popular tweets >2000RT)",
        "n_branching": None, "n_lo": None, "n_hi": None,
        "t_half_s": 2.0 * DAY, "t_half_raw": "tau_m = 2.02 +/- 0.66 d (synth)",
        "kernel": "power-law with cutoff: c0=6.49e-4/s, s0=300s, theta=0.242",
        "kappa_sign": "NA", "kappa_value": "NA",
        "n_events": "738 cascades, 166076 tweets", "obs_window": "Oct 7 - Nov 7 2011",
        "mu_background": "circadian time-dependent",
        "venue": "ICWSM", "authors": "Kobayashi, Lambiotte", "year": 2016,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1609/icwsm.v10i1.14717",
        "notes": "TiDeH; characteristic popularity decay ~2 days",
        "tier": "A",
    },
    # --- Crime (Mohler 2011 LA burglary) ---
    {
        "row_id": "A-005", "domain": "crime",
        "subdomain": "San Fernando Valley LA residential burglary (2004-05)",
        "n_branching": 0.20, "n_lo": None, "n_hi": None,
        "t_half_s": math.log(2) * 10.0 * DAY, "t_half_raw": "omega^-1 = ~10 d",
        "kernel": "exponential (theta*omega*exp(-omega*t))",
        "kappa_sign": "NA", "kappa_value": "NA",
        "n_events": 5376, "obs_window": "2004-2005",
        "mu_background": "mu=5.71/(km^2 day) sim baseline",
        "venue": "JASA", "authors": "Mohler, Short, Brantingham, Schoenberg, Tita", "year": 2011,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1198/jasa.2011.ap09546",
        "notes": "Sim params from Table 1; n derived from theta=0.20 (per-event offspring)",
        "tier": "A",
    },
    # --- Crime (Boyd & Molyneux mass shootings - 4 datasets) ---
    {
        "row_id": "A-006", "domain": "crime",
        "subdomain": "US mass shootings (Brady catalog)",
        "n_branching": 0.90, "n_lo": None, "n_hi": None,
        "t_half_s": None, "t_half_raw": "NA (23% offspring in 14d, long tail)",
        "kernel": "nonparametric (histogram, 2wk/3mo/6mo/+1y bins)",
        "kappa_sign": "+", "kappa_value": "k(m): 0.69(1)->1.44(5)->0.97(6-8)->0.59(9+)",
        "n_events": 477, "obs_window": "Feb 2005 - Jan 2013",
        "mu_background": "0.007-0.016/d",
        "venue": "PLOS ONE", "authors": "Boyd, Molyneux", "year": 2021,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1371/journal.pone.0248437",
        "notes": "n=expected offspring per event; mark dependence non-monotonic",
        "tier": "A",
    },
    {
        "row_id": "A-007", "domain": "crime",
        "subdomain": "US mass shootings (Stanford MSA, >=Jan 1999)",
        "n_branching": 0.72, "n_lo": None, "n_hi": None,
        "t_half_s": None, "t_half_raw": "NA (60% offspring in 14d)",
        "kernel": "nonparametric (histogram)",
        "kappa_sign": "+", "kappa_value": "k(m): 0.92(1)->1.38(5)->0.41(6-7)->0.13(7+)",
        "n_events": 262, "obs_window": "Jan 1999 - Jun 2016",
        "mu_background": "0.007-0.016/d",
        "venue": "PLOS ONE", "authors": "Boyd, Molyneux", "year": 2021,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1371/journal.pone.0248437",
        "notes": "Multi-dataset analysis row 2/4",
        "tier": "A",
    },
    {
        "row_id": "A-008", "domain": "crime",
        "subdomain": "US mass shootings (Gun Violence Archive)",
        "n_branching": 0.66, "n_lo": None, "n_hi": None,
        "t_half_s": None, "t_half_raw": "NA (~90% offspring in 14d)",
        "kernel": "nonparametric (histogram)",
        "kappa_sign": "+", "kappa_value": "k(m): 0.58(1)->0.87(5)->0.72(6-9)->0.13(10+)",
        "n_events": 2669, "obs_window": "Jan 2012-2021",
        "mu_background": "0.33/d",
        "venue": "PLOS ONE", "authors": "Boyd, Molyneux", "year": 2021,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1371/journal.pone.0248437",
        "notes": "Multi-dataset analysis row 3/4; loosest definition",
        "tier": "A",
    },
    {
        "row_id": "A-009", "domain": "crime",
        "subdomain": "US mass shootings (Mother Jones, strict)",
        "n_branching": 0.46, "n_lo": None, "n_hi": None,
        "t_half_s": None, "t_half_raw": "NA (98% offspring AFTER 14d)",
        "kernel": "nonparametric (histogram)",
        "kappa_sign": "+", "kappa_value": "k(m): 1.23(1)->0.19(7-10)->0.004(11-17)->0.29(18+)",
        "n_events": 118, "obs_window": "Jan 1999-present",
        "mu_background": "0.007-0.016/d",
        "venue": "PLOS ONE", "authors": "Boyd, Molyneux", "year": 2021,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1371/journal.pone.0248437",
        "notes": "Multi-dataset analysis row 4/4; strictest definition (4+ killed)",
        "tier": "A",
    },
    # --- Epidemiology (Sulem/Browning COVID, 10 countries) ---
    {
        "row_id": "A-010", "domain": "epidemiology",
        "subdomain": "COVID-19 daily deaths (Browning et al, geometric kernel)",
        "n_branching": None, "n_lo": None, "n_hi": None,
        "t_half_s": 2.5 * DAY, "t_half_raw": "~2.5 d (beta=0.4) or ~1 d (beta=0.9)",
        "kernel": "discrete-time geometric: beta*(1-beta)^(t-1)",
        "kappa_sign": "NA", "kappa_value": "NA",
        "n_events": "daily counts 10 countries", "obs_window": "Dec 2019 - Feb 2021",
        "mu_background": "import rate",
        "venue": "PLOS ONE", "authors": "Browning, Sulem, Mengersen, Rivoirard, Rousseau", "year": 2021,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1371/journal.pone.0250015",
        "notes": "Brazil/China/France/Germany/India/Italy/Spain/Sweden/UK/US; multi-phase",
        "tier": "A",
    },
    # --- Epidemiology (Mohler/Chiang COVID US counties) ---
    {
        "row_id": "A-011", "domain": "epidemiology",
        "subdomain": "COVID-19 county-level US (Mohler/Chiang Hawkes-mobility)",
        "n_branching": 1.5, "n_lo": 0.5, "n_hi": 2.5,
        "t_half_s": None, "t_half_raw": "NA (Weibull alpha,beta)",
        "kernel": "Weibull (alpha shape, beta scale)",
        "kappa_sign": "+", "kappa_value": "mobility covariates (Google index)",
        "n_events": "3217 counties, 2824 with >=10 cases", "obs_window": "Feb 2020 - Mar 2021",
        "mu_background": "imported infections mu_c",
        "venue": "Int J Forecasting", "authors": "Chiang, Liu, Mohler", "year": 2022,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1016/j.ijforecast.2021.07.001",
        "notes": "Dynamic R(t)=exp(theta.x), R0=2.5 early -> 1 post stay-home",
        "tier": "A",
    },
    # --- Seismology (Naylor ETAS.inlabru — though synthetic, anchors ETAS scale) ---
    {
        "row_id": "A-012", "domain": "seismology",
        "subdomain": "ETAS synthetic anchor (Naylor: Apennines/Ridgecrest-like)",
        "n_branching": None, "n_lo": None, "n_hi": None,
        "t_half_s": half_life_omori(c_days=0.11, p=1.08), "t_half_raw": "c=0.11d, p=1.08",
        "kernel": "ETAS Omori (c=0.11, p=1.08, alpha=2.29, K=0.089)",
        "kappa_sign": "+", "kappa_value": "alpha = 2.29 (magnitude productivity)",
        "n_events": "varies 117-2530", "obs_window": "1000-day synthetic",
        "mu_background": "mu=0.1/d",
        "venue": "Frontiers Appl Math Stat", "authors": "Naylor, Serafini, Lindgren, Main", "year": 2023,
        "peer_reviewed": "yes",
        "doi_or_url": "10.3389/fams.2023.1126759",
        "notes": "Synthetic but anchored to Apennines/Ridgecrest scale; method paper",
        "tier": "A",
    },
    # --- Seismology (Kwon Chile/Japan, 5 catalogs) ---
    {
        "row_id": "A-013", "domain": "seismology",
        "subdomain": "Chile A subduction (2001-2006, M>=4.0)",
        "n_branching": None, "n_lo": None, "n_hi": None,
        "t_half_s": half_life_omori(c_days=0.0327, p=1.0947), "t_half_raw": "c=0.0327d, p=1.0947 (synth fit)",
        "kernel": "kernel-based nonparametric ETAS",
        "kappa_sign": "+", "kappa_value": "magnitude productivity",
        "n_events": 1569, "obs_window": "2001-2006",
        "mu_background": "NA",
        "venue": "Spatial Statistics", "authors": "Kwon, Zheng, Jun", "year": 2023,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1016/j.spasta.2023.100728",
        "notes": "Kwon multi-catalog row 1/5",
        "tier": "A",
    },
    {
        "row_id": "A-014", "domain": "seismology",
        "subdomain": "Chile B subduction (2010-2015 post-2010 M8.8)",
        "n_branching": None, "n_lo": None, "n_hi": None,
        "t_half_s": None, "t_half_raw": "Omori (region-specific)",
        "kernel": "kernel-based nonparametric ETAS",
        "kappa_sign": "+", "kappa_value": "magnitude productivity",
        "n_events": 3110, "obs_window": "Feb 2010 - Sep 2015",
        "mu_background": "NA",
        "venue": "Spatial Statistics", "authors": "Kwon, Zheng, Jun", "year": 2023,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1016/j.spasta.2023.100728",
        "notes": "Kwon multi-catalog row 2/5; post-M8.8 Maule sequence",
        "tier": "A",
    },
    {
        "row_id": "A-015", "domain": "seismology",
        "subdomain": "Chile C subduction (2015-2021 post-M8.3)",
        "n_branching": None, "n_lo": None, "n_hi": None,
        "t_half_s": None, "t_half_raw": "Omori (region-specific)",
        "kernel": "kernel-based nonparametric ETAS",
        "kappa_sign": "+", "kappa_value": "magnitude productivity",
        "n_events": 2552, "obs_window": "Sep 2015 - Sep 2021",
        "mu_background": "NA",
        "venue": "Spatial Statistics", "authors": "Kwon, Zheng, Jun", "year": 2023,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1016/j.spasta.2023.100728",
        "notes": "Kwon multi-catalog row 3/5",
        "tier": "A",
    },
    {
        "row_id": "A-016", "domain": "seismology",
        "subdomain": "Japan A region (2003-2008)",
        "n_branching": None, "n_lo": None, "n_hi": None,
        "t_half_s": None, "t_half_raw": "Omori (region-specific)",
        "kernel": "kernel-based nonparametric ETAS",
        "kappa_sign": "+", "kappa_value": "magnitude productivity",
        "n_events": 1327, "obs_window": "2003-2008",
        "mu_background": "NA",
        "venue": "Spatial Statistics", "authors": "Kwon, Zheng, Jun", "year": 2023,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1016/j.spasta.2023.100728",
        "notes": "Kwon multi-catalog row 4/5",
        "tier": "A",
    },
    {
        "row_id": "A-017", "domain": "seismology",
        "subdomain": "Japan B post-Tohoku M9.1 (2011-2017)",
        "n_branching": None, "n_lo": None, "n_hi": None,
        "t_half_s": None, "t_half_raw": "Omori (region-specific)",
        "kernel": "kernel-based nonparametric ETAS",
        "kappa_sign": "+", "kappa_value": "magnitude productivity",
        "n_events": 7420, "obs_window": "Mar 2011 - Mar 2017",
        "mu_background": "NA",
        "venue": "Spatial Statistics", "authors": "Kwon, Zheng, Jun", "year": 2023,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1016/j.spasta.2023.100728",
        "notes": "Kwon multi-catalog row 5/5; post-Tohoku aftershock sequence",
        "tier": "A",
    },
    # --- Heliophysics (Ross 2020 solar self-excitation) ---
    {
        "row_id": "A-018", "domain": "heliophysics",
        "subdomain": "Solar flares (GOES catalog, waiting times)",
        "n_branching": None, "n_lo": None, "n_hi": None,
        "t_half_s": None, "t_half_raw": "NA (qualitative)",
        "kernel": "self-exciting (parametric Hawkes)",
        "kappa_sign": "+", "kappa_value": "burstiness above Poisson baseline",
        "n_events": "GOES catalog (~10^4)", "obs_window": "multi-cycle",
        "mu_background": "solar-cycle modulated",
        "venue": "Physica A", "authors": "Ross", "year": 2020,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1016/j.physa.2020.124775",
        "notes": "Rejects pure non-stationary Poisson; significant self-excitation",
        "tier": "A",
    },
    {
        "row_id": "A-019", "domain": "heliophysics",
        "subdomain": "Solar flares (Rivera/Johnson clustering ~3h scale)",
        "n_branching": None, "n_lo": None, "n_hi": None,
        "t_half_s": 3.0 * HOUR, "t_half_raw": "~3 h characteristic recurrence",
        "kernel": "nonparametric clustering signature",
        "kappa_sign": "+", "kappa_value": "memory above Bayesian-block surrogate",
        "n_events": "~10^4 GOES flares", "obs_window": "SC23+SC24",
        "mu_background": "rate-variable Bayesian blocks",
        "venue": "ApJ Letters", "authors": "Rivera, Johnson, Homan, Wing", "year": 2022,
        "peer_reviewed": "yes",
        "doi_or_url": "10.3847/2041-8213/ac8de9",
        "notes": "Surrogate-based nonparametric test; ~3h recurrence excess",
        "tier": "A",
    },
]

# ============================================================
# TIER B (from wide_research; n value confirmed)
# ============================================================
TIER_B = [
    {
        "row_id": "B-001", "domain": "finance",
        "subdomain": "China A-share stocks (108 firms)",
        "n_branching": 0.81, "n_lo": None, "n_hi": None,
        "t_half_s": None, "t_half_raw": "NA",
        "kernel": "power-law (nonparametric)",
        "kappa_sign": "NA", "kappa_value": "NA",
        "n_events": "108 stocks", "obs_window": "NA",
        "mu_background": "NA",
        "venue": "Eur J Finance", "authors": "Zhuo et al.", "year": 2023,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1080/1351847X.2023.2251531",
        "notes": "Mean internal branching ratio ~0.81; BIC selected power-law",
        "tier": "B",
    },
    {
        "row_id": "B-002", "domain": "social",
        "subdomain": "Twitter URL cascades (Gleeson/Onaga 2020)",
        "n_branching": 0.90, "n_lo": None, "n_hi": None,
        "t_half_s": None, "t_half_raw": "NA",
        "kernel": "branching process (Hawkes-related)",
        "kappa_sign": "NA", "kappa_value": "NA",
        "n_events": "Twitter URL cascades", "obs_window": "NA",
        "mu_background": "NA",
        "venue": "J Complex Networks", "authors": "Gleeson, Onaga et al.", "year": 2021,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1093/comnet/cnab002",
        "notes": "Effective branching number 0.90; branching/Hawkes hybrid",
        "tier": "B",
    },
    {
        "row_id": "B-003", "domain": "social",
        "subdomain": "Twitter retweet popularity prediction (Mishra/Rizoiu)",
        "n_branching": 0.12, "n_lo": None, "n_hi": None,
        "t_half_s": None, "t_half_raw": "NA",
        "kernel": "marked Hawkes power-law: phi_m(tau)=kappa*m^beta*(tau+c)^-(1+theta)",
        "kappa_sign": "+", "kappa_value": "kappa=0.75 example cascade",
        "n_events": "Twitter cascades", "obs_window": "NA",
        "mu_background": "NA",
        "venue": "TKDE", "authors": "Mishra, Rizoiu, Xie", "year": 2018,
        "peer_reviewed": "yes",
        "doi_or_url": "10.1109/TKDE.2018.2885271",
        "notes": "Hybrid Hawkes+predictive overlay; n=0.12 small-tweet regime",
        "tier": "B",
    },
]

# ============================================================
# Assemble and write
# ============================================================
def main():
    cols = [
        "row_id", "tier", "domain", "subdomain", "n_branching", "n_lo", "n_hi",
        "t_half_s", "t_half_raw", "kernel", "kappa_sign", "kappa_value",
        "n_events", "obs_window", "mu_background",
        "venue", "authors", "year", "peer_reviewed", "doi_or_url", "notes",
    ]

    all_rows = OWN + TIER_A + TIER_B

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in all_rows:
            w.writerow({k: r.get(k, "") for k in cols})

    print(f"Wrote {len(all_rows)} rows to {OUT}\n")

    from collections import Counter
    domain_counts = Counter(r["domain"] for r in all_rows)
    tier_counts = Counter(r["tier"] for r in all_rows)
    print("Domain breakdown:")
    for d, c in sorted(domain_counts.items(), key=lambda x: -x[1]):
        print(f"  {d:25s} {c}")
    print("\nTier breakdown:")
    for t, c in tier_counts.items():
        print(f"  {t:5s} {c}")

    n_valid = sum(1 for r in all_rows if r.get("n_branching") is not None)
    t_valid = sum(1 for r in all_rows if r.get("t_half_s") is not None)
    both = sum(1 for r in all_rows if r.get("n_branching") is not None and r.get("t_half_s") is not None)
    print(f"\nRows with n_branching:   {n_valid}")
    print(f"Rows with t_half:        {t_valid}")
    print(f"Rows with BOTH:          {both}")
    print(f"\nDomains covered:         {len(domain_counts)}")

if __name__ == "__main__":
    main()
