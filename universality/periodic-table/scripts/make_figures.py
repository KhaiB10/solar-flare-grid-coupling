"""
Periodic Table of Self-Exciting Systems — figures.

Random seed: 20260523.
"""
import csv
import math
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from collections import Counter

np.random.seed(20260523)

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "data" / "periodic_table_v1.csv"
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(exist_ok=True)


def load_rows():
    rows = []
    with open(CSV_PATH) as f:
        for r in csv.DictReader(f):
            try:
                r["n_branching_f"] = float(r["n_branching"]) if r["n_branching"] else None
            except ValueError:
                r["n_branching_f"] = None
            try:
                r["t_half_s_f"] = float(r["t_half_s"]) if r["t_half_s"] else None
            except ValueError:
                r["t_half_s_f"] = None
            rows.append(r)
    return rows


DOMAIN_COLOR = {
    "heliophysics": "#f59f00",
    "tropical_cyclones": "#1c7ed6",
    "seismology": "#7048e8",
    "finance": "#2f9e44",
    "social": "#e64980",
    "crime": "#fa5252",
    "epidemiology": "#15aabf",
}


# ----- Figure 1: log(n) vs log(t_half) scatter, colored by domain -----
def figure_scatter():
    rows = load_rows()
    fig, ax = plt.subplots(figsize=(10, 7))

    has_both = [r for r in rows if r["n_branching_f"] is not None and r["t_half_s_f"] is not None]
    has_n_only = [r for r in rows if r["n_branching_f"] is not None and r["t_half_s_f"] is None]
    has_t_only = [r for r in rows if r["n_branching_f"] is None and r["t_half_s_f"] is not None]

    for r in has_both:
        x = r["t_half_s_f"]
        y = r["n_branching_f"]
        c = DOMAIN_COLOR.get(r["domain"], "#888")
        ax.scatter(x, y, s=160, c=c, edgecolors="black", linewidths=1.3, alpha=0.9, zorder=3)
        # Label our own rows
        if r["row_id"].startswith("OWN"):
            label = "Solar G4+" if "solar" in r["subdomain"].lower() else "Hurricane Cat3+"
            ax.annotate(label, (x, y),
                        xytext=(8, 8), textcoords="offset points",
                        fontsize=10, fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=c, lw=1.5))

    # Show n-only rows as horizontal bars near a fake low t_half
    for r in has_n_only:
        c = DOMAIN_COLOR.get(r["domain"], "#888")
        ax.scatter(1e-5, r["n_branching_f"], s=80, c=c, marker=">",
                   edgecolors="black", linewidths=0.8, alpha=0.5, zorder=2)

    # Show t-only as vertical bars near fake n
    for r in has_t_only:
        c = DOMAIN_COLOR.get(r["domain"], "#888")
        ax.scatter(r["t_half_s_f"], 0.001, s=80, c=c, marker="^",
                   edgecolors="black", linewidths=0.8, alpha=0.5, zorder=2)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Excitation half-life $t_{1/2}$ (seconds)", fontsize=13)
    ax.set_ylabel(r"Branching ratio $n$", fontsize=13)
    ax.set_title("Periodic Table of Self-Exciting Systems\nlog(n) vs log(half-life), curated peer-reviewed Hawkes fits",
                 fontsize=13, pad=14)

    # Reference lines: criticality, stability
    ax.axhline(1.0, color="red", linestyle="--", alpha=0.6, label=r"$n=1$ (critical)")
    ax.axhline(0.5, color="grey", linestyle=":", alpha=0.5, label=r"$n=0.5$ (half-supercritical)")

    # Human time-scale anchors
    for t_s, label, col in [(1, "1 s", "#ddd"), (60, "1 min", "#ddd"), (3600, "1 h", "#ddd"),
                            (86400, "1 d", "#ddd"), (86400*30, "1 mo", "#ddd"),
                            (86400*365, "1 yr", "#ddd")]:
        ax.axvline(t_s, color=col, linestyle="-", alpha=0.4, lw=0.8, zorder=0)
        ax.text(t_s, 1.4e-3, label, rotation=90, fontsize=8, ha="right", va="bottom", color="#666", alpha=0.7)

    # Domain legend
    handles = []
    for d, c in DOMAIN_COLOR.items():
        if any(r["domain"] == d for r in rows):
            handles.append(plt.Line2D([0], [0], marker="o", color="w",
                                      markerfacecolor=c, markeredgecolor="black",
                                      markersize=10, label=d.replace("_", " ")))
    handles.append(plt.Line2D([0], [0], linestyle="--", color="red", label=r"$n=1$ critical"))
    ax.legend(handles=handles, loc="lower left", fontsize=9, framealpha=0.9)

    ax.grid(True, which="both", alpha=0.3)
    ax.set_xlim(1e-5, 1e9)
    ax.set_ylim(1e-3, 2)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "01_periodic_table_scatter.png", dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Wrote {FIG_DIR / '01_periodic_table_scatter.png'}")


# ----- Figure 2: domain histograms of branching ratio -----
def figure_n_hist():
    rows = load_rows()
    # Collect n values per domain
    by_domain = {}
    for r in rows:
        if r["n_branching_f"] is None:
            continue
        by_domain.setdefault(r["domain"], []).append(r["n_branching_f"])

    domains = sorted(by_domain.keys(), key=lambda d: -len(by_domain[d]))
    n_cols = 2
    n_rows = (len(domains) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(11, 2.6 * n_rows), squeeze=False)

    bins = np.linspace(0, 1.1, 12)
    for i, d in enumerate(domains):
        ax = axes[i // n_cols][i % n_cols]
        vals = by_domain[d]
        ax.hist(vals, bins=bins, color=DOMAIN_COLOR.get(d, "#888"),
                edgecolor="black", linewidth=0.5, alpha=0.85)
        ax.axvline(1.0, color="red", linestyle="--", alpha=0.7)
        ax.set_title(f"{d.replace('_', ' ')} — N={len(vals)}", fontsize=11)
        ax.set_xlim(0, 1.15)
        ax.set_xlabel("branching ratio n")
        ax.set_ylabel("count")
        ax.grid(True, alpha=0.3)

    # Hide unused axes
    for j in range(len(domains), n_rows * n_cols):
        axes[j // n_cols][j % n_cols].axis("off")

    fig.suptitle("Branching ratio distribution by domain (curated peer-reviewed fits)", fontsize=13, y=1.0)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "02_n_histograms.png", dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Wrote {FIG_DIR / '02_n_histograms.png'}")


# ----- Figure 3: kappa-sign tally -----
def figure_kappa():
    rows = load_rows()
    counts = Counter()
    for r in rows:
        sign = (r.get("kappa_sign") or "").strip()
        if sign in ("+", "-", "0", "NA"):
            counts[(r["domain"], sign)] += 1

    domains = sorted({d for d, _ in counts.keys()})
    signs = ["+", "-", "0", "NA"]
    sign_color = {"+": "#37b24d", "-": "#f03e3e", "0": "#888", "NA": "#dee2e6"}

    fig, ax = plt.subplots(figsize=(10, max(4, 0.6 * len(domains) + 2)))
    y = np.arange(len(domains))
    bottoms = np.zeros(len(domains))

    for sign in signs:
        widths = np.array([counts.get((d, sign), 0) for d in domains])
        ax.barh(y, widths, left=bottoms, color=sign_color[sign],
                edgecolor="black", linewidth=0.7, label=f"kappa {sign}")
        # Annotate
        for yi, w, b in zip(y, widths, bottoms):
            if w > 0:
                ax.text(b + w / 2, yi, str(w), ha="center", va="center",
                        fontsize=10, fontweight="bold",
                        color="white" if sign in ("+", "-") else "black")
        bottoms += widths

    ax.set_yticks(y)
    ax.set_yticklabels([d.replace("_", " ") for d in domains])
    ax.set_xlabel("number of rows")
    ax.set_title("Mark/covariate effect sign (kappa) by domain", fontsize=13)
    ax.legend(loc="lower right", fontsize=10, framealpha=0.95)
    ax.grid(True, alpha=0.3, axis="x")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "03_kappa_signs.png", dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Wrote {FIG_DIR / '03_kappa_signs.png'}")


# ----- Figure 4: half-life range by domain -----
def figure_thalf():
    rows = load_rows()
    by_domain = {}
    for r in rows:
        if r["t_half_s_f"] is None:
            continue
        by_domain.setdefault(r["domain"], []).append((r["t_half_s_f"], r["subdomain"][:30]))

    domains = sorted(by_domain.keys(), key=lambda d: np.median([v[0] for v in by_domain[d]]))
    fig, ax = plt.subplots(figsize=(10, max(4, 0.6 * len(domains) + 2)))

    for i, d in enumerate(domains):
        vals = [v[0] for v in by_domain[d]]
        ax.scatter(vals, [i] * len(vals), s=140, c=DOMAIN_COLOR.get(d, "#888"),
                   edgecolors="black", linewidths=1.1, alpha=0.9, zorder=3)
        if len(vals) > 1:
            ax.plot([min(vals), max(vals)], [i, i], color="black", alpha=0.4, zorder=2)

    ax.set_xscale("log")
    ax.set_yticks(range(len(domains)))
    ax.set_yticklabels([d.replace("_", " ") for d in domains])
    ax.set_xlabel(r"$t_{1/2}$ (seconds)")
    ax.set_title("Excitation half-life range by domain (log scale)", fontsize=13)

    for t_s, label in [(1, "1 s"), (60, "1 min"), (3600, "1 h"),
                        (86400, "1 d"), (86400*30, "1 mo"), (86400*365, "1 yr")]:
        ax.axvline(t_s, color="#ddd", lw=0.8, zorder=0)
        ax.text(t_s, len(domains) - 0.4, label, rotation=90, fontsize=8, ha="right", va="top", color="#666")

    ax.grid(True, alpha=0.3, axis="x")
    ax.set_xlim(1e-5, 1e9)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "04_thalf_by_domain.png", dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Wrote {FIG_DIR / '04_thalf_by_domain.png'}")


if __name__ == "__main__":
    figure_scatter()
    figure_n_hist()
    figure_kappa()
    figure_thalf()
