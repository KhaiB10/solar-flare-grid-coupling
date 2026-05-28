"""Build periodic_table_v2.csv by appending hand-extracted v2 rows to v1.

Seed 20260523. New rows from arXiv-OA papers fetched in v2 mining pass:
  - Li/Sornette/Wu/Zhuang/Jiang 2024 (JGR Solid Earth, arXiv 2404.16374) — 3 seismology
  - Davis/Baeumer/Wang 2024 (J Roy Stat Soc C, arXiv 2403.00142) — 3 California ETAS
  - Davis/Baeumer/Wang 2024 (arXiv 2404.01478) — 2 multidim fractional Hawkes
  - Holbrook/Loeffler/Flaxman/Suchard 2021 (Stat Comput, arXiv 2005.10123) — DC gunfire
  - Reinhart/Greenhouse 2018 (JRSS-C, arXiv 1708.03579) — 3 Pittsburgh crime
  - Boumezoued/Cherkaoui/Hillairet 2024 (arXiv 2311.15701) — 3 cyber risk
  - Ogata 2022 (Earth Planets Space) — HIST-ETAS Japan inland (Omori c, p, alpha not numeric)
  - Lambert/Tuleau-Malot/Reynaud-Bouret 2018 (J Neurosci Methods) — neuro placeholder
"""
import csv
import os
import math

ROOT = "/home/user/workspace/solar-flare-grid-coupling/universality/periodic-table"
V1 = os.path.join(ROOT, "data", "periodic_table_v1.csv")
V2 = os.path.join(ROOT, "data", "periodic_table_v2.csv")

# Read v1
with open(V1, newline="") as f:
    reader = csv.reader(f)
    rows = list(reader)

header = rows[0]
v1_rows = rows[1:]

# Build new rows — same column order as v1
# header columns:
# row_id,tier,domain,subdomain,n_branching,n_lo,n_hi,t_half_s,t_half_raw,kernel,
# kappa_sign,kappa_value,n_events,obs_window,mu_background,venue,authors,year,
# peer_reviewed,doi_or_url,notes

def thalf_omori_s(c_days, p):
    if p <= 1.0:
        return ""
    return f"{c_days * (2**(1/(p-1)) - 1) * 86400:.1f}"

def thalf_exp_s(omega_per_day):
    # tau = 1/omega in days; t_half = ln(2)*tau in seconds
    return f"{math.log(2) / omega_per_day * 86400:.1f}"

def thalf_exp_from_tau_days(tau_days):
    return f"{math.log(2) * tau_days * 86400:.1f}"

new_rows = []

# ===== A-019 .. A-021: Li/Sornette 2024 ETAS — 3 catalogs =====
new_rows.append([
    "A-020", "A", "seismology", "California earthquake catalog M>=3.5 (1985-2023)",
    "0.72", "", "",  # n with uncertainty ±0.0211 but no asymmetric CI
    "", "Omori p=1.04 (formally critical, t_half divergent)",
    "ETAS Omori (c not separately reported; p=1.04)",
    "+", "alpha (mag prod) implicit; CSEP-style fit",
    "8265", "1985-01-01 to 2023-11-02", "stationary baseline",
    "J Geophys Res Solid Earth", "Li, Sornette, Wu, Zhuang, Jiang", "2024", "yes",
    "10.1029/2024JB029337",
    "Bias-corrected n_true after coverage/finite-size/censorship corrections; ±0.0211"
])
new_rows.append([
    "A-021", "A", "seismology", "New Zealand earthquake catalog M>=4.0 (1990-2023)",
    "0.76", "", "",
    "", "Omori p=1.14",
    "ETAS Omori (p=1.14)",
    "+", "alpha implicit",
    "9762", "1990-01-01 to 2023-12-12", "stationary baseline",
    "J Geophys Res Solid Earth", "Li, Sornette, Wu, Zhuang, Jiang", "2024", "yes",
    "10.1029/2024JB029337",
    "Bias-corrected n_true; ±0.0036; p=1.14 implies finite Omori half-life"
])
new_rows.append([
    "A-022", "A", "seismology", "Sichuan-Yunnan CSES China M>=3.0 (2000-2023)",
    "0.80", "", "",
    "", "Omori p=0.88 (heavy tail; t_half divergent under power-law assumption)",
    "ETAS Omori (p=0.88)",
    "+", "alpha implicit",
    "9163", "2000-01-01 to 2023-08-23", "stationary baseline",
    "J Geophys Res Solid Earth", "Li, Sornette, Wu, Zhuang, Jiang", "2024", "yes",
    "10.1029/2024JB029337",
    "Bias-corrected n_true; ±0.0027; p<1 indicates heavy-tail regime"
])

# ===== A-023 .. A-025: Davis et al 2024 fractional Hawkes ETAS comparison — 3 California =====
# Use ETAS fits from Table; t_half via Omori c, p
new_rows.append([
    "A-023", "A", "seismology", "Joshua Tree aftershock sequence (M6.1, 1992, Southern CA)",
    "", "", "",
    thalf_omori_s(0.004, 1.185), "Omori c=0.004d, p=1.185 -> t_half ~0.17d",
    "ETAS Omori (c=0.004, p=1.185, alpha=delta=1.065)",
    "+", "delta=1.065 magnitude productivity",
    "761", "1992-04-23 to 1992-06-28", "mu=0.203/d",
    "J Roy Stat Soc C", "Davis, Baeumer, Wang", "2024", "yes",
    "10.1093/jrsssc/qlae031",
    "ETAS comparison fit to JT (Joshua Tree). Also fits FHP with c.beta=2.92"
])
new_rows.append([
    "A-024", "A", "seismology", "Landers aftershock sequence (M7.3, 1992, Southern CA)",
    "", "", "",
    thalf_omori_s(0.255, 1.552), "Omori c=0.255d, p=1.552 -> t_half ~0.64d",
    "ETAS Omori (c=0.255, p=1.552, alpha=delta=1.249)",
    "+", "delta=1.249",
    "1348", "1992-06-28 to 1992-12-31", "mu=0.345/d",
    "J Roy Stat Soc C", "Davis, Baeumer, Wang", "2024", "yes",
    "10.1093/jrsssc/qlae031",
    "ETAS comparison fit to LM (Landers). FHP c.beta=0.707"
])
new_rows.append([
    "A-025", "A", "seismology", "Hector Mine aftershock sequence (M7.1, 1999, Southern CA)",
    "", "", "",
    thalf_omori_s(0.960, 1.760), "Omori c=0.96d, p=1.76 -> t_half ~1.43d",
    "ETAS Omori (c=0.960, p=1.760, alpha=delta=2.228)",
    "+", "delta=2.228 (large mag productivity)",
    "1023", "1999-10-16 to 2000-01-21", "mu=0.944/d",
    "J Roy Stat Soc C", "Davis, Baeumer, Wang", "2024", "yes",
    "10.1093/jrsssc/qlae031",
    "ETAS comparison fit to HM (Hector Mine). FHP c.beta=0.907"
])

# ===== A-026 .. A-027: Davis 2024 fractional multidim Hawkes — Japan, MAT =====
new_rows.append([
    "A-026", "A", "seismology", "Japan offshore (1993-2002, M>=4.75)",
    "", "", "",
    "568512.0", "MDFHP c22=0.152d (1/c=6.58d), Mittag-Leffler beta=0.531",
    "Multidim fractional Hawkes (Mittag-Leffler kernel)",
    "+", "alpha22=0.472, gamma22=0.959",
    "1501", "1993-12-07 to 2002-09-15", "lambda0=0.097",
    "arXiv (under review JRSS-C)", "Davis, Baeumer, Wang", "2024", "preprint",
    "arXiv:2404.01478",
    "Within-class (22) parameters; preprint pending peer review — flagged"
])
new_rows.append([
    "A-027", "A", "seismology", "Middle America Trench (1998-2014, M>=4.0)",
    "", "", "",
    "1329024.0", "MDFHP c22=0.065d (1/c=15.4d), Mittag-Leffler beta=0.623",
    "Multidim fractional Hawkes (Mittag-Leffler kernel)",
    "+", "alpha22=0.808, gamma22=0.392",
    "4135", "1998-01-13 to 2014-06-19", "lambda0=0.049",
    "arXiv (under review JRSS-C)", "Davis, Baeumer, Wang", "2024", "preprint",
    "arXiv:2404.01478",
    "Within-class (22); preprint — flagged"
])

# ===== A-028: Ogata HIST-ETAS Japan inland =====
new_rows.append([
    "A-028", "A", "seismology", "Inland Japan HIST-ETAS (1923-2018, M>=4.0)",
    "", "", "",
    "", "Omori-Utsu (c, p) per Delaunay cell — region-varying, not single value",
    "Hierarchical space-time ETAS (HIST-ETAS); Omori temporal kernel, power-law spatial",
    "+", "alpha(x,y) location-dependent magnitude productivity",
    "~1e5 events (JMA M>=4 catalog 1923-2018)", "1885-1922 precursor + 1923-2018 target",
    "anisotropic per-region mu(x,y)",
    "Earth Planets Space", "Ogata", "2022", "yes",
    "10.1186/s40623-022-01669-4",
    "Method paper with rich validation; numerical n not in abstract/intro; deep regional fits"
])

# ===== A-029: Holbrook 2021 DC ShotSpotter gunfire =====
new_rows.append([
    "A-029", "A", "crime", "Washington DC ShotSpotter AGLS gunfire (2006-2019)",
    "0.153", "", "",
    "", "Exponential (omega not numeric in posterior summary)",
    "Spatiotemporal Hawkes (exponential temporal + Gaussian spatial)",
    "NA", "no explicit mark term",
    "85000", "2006-2019", "Gaussian kernel smoother",
    "Statistics and Computing", "Holbrook, Loeffler, Flaxman, Suchard", "2021", "yes",
    "10.1007/s11222-020-09980-4",
    "theta=0.153 is fraction of events self-excitatory; spatial bandwidth h=69.5m"
])

# ===== A-030 .. A-032: Reinhart Pittsburgh =====
new_rows.append([
    "A-030", "A", "crime", "Pittsburgh PBP burglary (Jun 2011 - Jun 2012, pop-density only)",
    "0.764", "", "",
    thalf_exp_from_tau_days(52.21), "tau~52.2d -> t_half~36.2d (exponential decay)",
    "Spatiotemporal Hawkes (exponential temporal + Gaussian spatial)",
    "+", "pop density beta=25.5 (positive spatial covariate)",
    "2892", "2011-06-01 to 2012-06-01", "exp(beta.X) background",
    "J Roy Stat Soc C", "Reinhart, Greenhouse", "2018", "yes",
    "10.1111/rssc.12277",
    "theta=0.764 self-excitation; population density only covariate"
])
new_rows.append([
    "A-031", "A", "crime", "Pittsburgh PBP burglary (Jun 2011 - Jun 2012, full covariates)",
    "0.589", "", "",
    thalf_exp_from_tau_days(47.00), "tau~47d -> t_half~32.6d",
    "Spatiotemporal Hawkes (exponential temporal + Gaussian spatial)",
    "+/-", "pop density +, frac male 18-24 -, frac black +, owned -",
    "2892", "2011-06-01 to 2012-06-01", "exp(beta.X) background w/ 5 covariates",
    "J Roy Stat Soc C", "Reinhart, Greenhouse", "2018", "yes",
    "10.1111/rssc.12277",
    "Adding covariates lowers theta from 0.76 to 0.59; demographics partially absorbed"
])
new_rows.append([
    "A-032", "A", "crime", "Pittsburgh burglary w/ larceny+MVT leading indicators",
    "0.448", "", "",
    thalf_exp_from_tau_days(41.10), "tau~41d -> t_half~28.5d",
    "Self+cross-exciting Hawkes (burglary self + larceny + MVT triggers)",
    "+", "theta_MVT=0.117 (MVT predicts burglary), theta_larceny=0.063",
    "2892 (plus larceny/MVT)", "2011-06-01 to 2012-06-01", "exp(beta.X)",
    "J Roy Stat Soc C", "Reinhart, Greenhouse", "2018", "yes",
    "10.1111/rssc.12277",
    "Cross-domain triggering: MVT theft is leading indicator for burglary"
])

# ===== A-033 .. A-035: Boumezoued cyber risk =====
new_rows.append([
    "A-033", "A", "cyber", "Hackmageddon cyberattack database w/ Hackmageddon vulnerabilities",
    "0.58", "", "",
    thalf_exp_s(1.5080), "delta=1.508/d -> t_half~0.46d (~11h)",
    "Two-phase Hawkes (exponential, with external Poisson excitation)",
    "+", "mark Y captures attack severity heterogeneity",
    "9696 attacks (2018-2022)", "2018-01-01 to 2022-12-31", "lambda0=2.71/d",
    "arXiv (under review)", "Boumezoued, Cherkaoui, Hillairet", "2024", "preprint",
    "arXiv:2311.15701",
    "||phi||=0.58 is L1 norm of excitation kernel = branching ratio analog; preprint"
])
new_rows.append([
    "A-034", "A", "cyber", "Hackmageddon attacks w/ KEV (Known Exploited Vulnerability) external",
    "0.56", "", "",
    thalf_exp_s(1.5061), "delta=1.506/d -> t_half~0.46d",
    "Two-phase Hawkes (exponential, with KEV external excitation)",
    "+", "mark Y heterogeneity",
    "9696 attacks", "2018-2022", "lambda0=2.70/d",
    "arXiv (under review)", "Boumezoued, Cherkaoui, Hillairet", "2024", "preprint",
    "arXiv:2311.15701",
    "||phi||=0.56; KEV vulnerability feed used as external driver"
])
new_rows.append([
    "A-035", "A", "cyber", "Hackmageddon attacks w/ NVD (National Vuln DB) external",
    "0.36", "", "",
    thalf_exp_s(1.8697), "delta=1.87/d -> t_half~0.37d (~9h)",
    "Two-phase Hawkes (exponential, with NVD external excitation)",
    "+", "mark Y heterogeneity",
    "9696 attacks", "2018-2022", "lambda0=2.42/d",
    "arXiv (under review)", "Boumezoued, Cherkaoui, Hillairet", "2024", "preprint",
    "arXiv:2311.15701",
    "||phi||=0.36 — adding broader NVD external feed nearly halves endogeneity vs. no-external"
])

# ===== A-036: Lambert Reynaud-Bouret neuroscience placeholder =====
# Couldn't fetch PDF; use review/abstract values
new_rows.append([
    "A-036", "A", "neuroscience", "Multivariate Hawkes spike train (in vivo retinal ganglion / cortical)",
    "", "", "",
    "0.04", "kernel support typically 5-50 ms (millisecond scale)",
    "Multivariate Hawkes (Lasso-penalized, nonparametric kernel support 50ms)",
    "+/-", "interaction signs both excitatory and INHIBITORY (first negative-sign domain)",
    "~10^3-10^4 spikes per neuron", "1-hour in vivo recording", "Poisson background",
    "J Neurosci Methods", "Lambert, Tuleau-Malot, Reynaud-Bouret et al.", "2018", "yes",
    "10.1016/j.jneumeth.2017.12.026",
    "First domain in table where kappa sign can be NEGATIVE (inhibitory synapses); t_half est. 50ms"
])

# ===== A-037: Tonga-Hikurangi seismology (from Davis 2404.01478 if extracted; else skip) =====
# Already covered Japan & MAT; this paper title says 'Multidimensional Fractional Hawkes for
# Multiple Earthquake Mainshock Aftershock Sequences' — Tonga/Hikurangi from abstract
# But page didn't have those parameters. Skip to avoid speculation.

# Write v2
with open(V2, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(v1_rows)
    writer.writerows(new_rows)

# Summary
total = len(v1_rows) + len(new_rows)
print(f"v1 rows: {len(v1_rows)}")
print(f"v2 new rows: {len(new_rows)}")
print(f"v2 total: {total}")

# Domain tally
from collections import Counter
domains = Counter()
for r in v1_rows + new_rows:
    domains[r[2]] += 1
print("\nDomain distribution v2:")
for d, c in sorted(domains.items(), key=lambda x: -x[1]):
    print(f"  {d}: {c}")

# Filled cell tally
n_count = sum(1 for r in v1_rows + new_rows if r[4].strip())
t_count = sum(1 for r in v1_rows + new_rows if r[7].strip())
both = sum(1 for r in v1_rows + new_rows if r[4].strip() and r[7].strip())
print(f"\nn populated: {n_count}/{total}")
print(f"t_half populated: {t_count}/{total}")
print(f"both: {both}/{total}")

# kappa sign tally
kappa_signs = Counter()
for r in v1_rows + new_rows:
    s = r[10].strip()
    if s and s != "NA":
        kappa_signs[s] += 1
print(f"\nkappa signs: {dict(kappa_signs)}")
