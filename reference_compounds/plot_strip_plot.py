"""Strip plot comparing the experimental hits of each AL round with the reference compounds.

Rows: targets. Columns: KD, ligand efficiency, molecular weight and linearity ratio.
The reference compounds are shown as a shaded band (range and interquartile range)
plus their median.

Run compute_properties.py first.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.transforms as mtransforms
from matplotlib.lines import Line2D

OUT = Path(__file__).resolve().parent / "results"

# Define the targets, the colors for each round, and the metrics to plot
TARGETS = {"IL1b": "IL-1β", "IL2": "IL-2", "IL36g": "IL-36γ"}
ROUND_COLORS = {"IL1b": ["#F4A97F", "#D65F1A"], "IL2": ["#82C882"],
                "IL36g": ["#C299D6", "#9B5FBF", "#7A3BAA", "#551F8A"]}
METRICS = {"KD_uM": "K$_D$ (µM)", "LigandEfficiency": "Ligand Efficiency (LE)",
           "MW": "Molecular Weight (Da)", "linearity_ratio": "Linearity Ratio"}

# Only active compounds are plotted: reference compounds with KD < 10 µM and
# experimental hits with KD < 1 mM.
REF_KD_CUTOFF = 10.0     # µM
HIT_KD_CUTOFF = 1000.0   # µM

# Define the plot styles
REF_COLOR, BAND_COLOR, BAND_ALPHA = "#4878CF", "#C8D8F0", 0.3

plt.rcParams.update({"font.size": 13, "axes.spines.top": False, "axes.spines.right": False})

# Plot
fig, axes = plt.subplots(3, 4, figsize=(17.5, 12), gridspec_kw={"hspace": 0.5, "wspace": 0.12})
fig.subplots_adjust(left=0.08, right=0.97, top=0.92, bottom=0.07)
for i, (target, target_label) in enumerate(TARGETS.items()):
    ref = pd.read_csv(OUT / f"{target}_reference_properties.csv").drop_duplicates("SMILES")
    ref = ref[ref["KD_uM"] < REF_KD_CUTOFF]

    rounds = {}
    for csv in sorted(OUT.glob(f"{target}_round*_properties.csv")):
        rnd = int(csv.stem.split("_round")[1].split("_")[0])
        hits = pd.read_csv(csv).drop_duplicates("SMILES")
        hits = hits[hits["KD_uM"] < HIT_KD_CUTOFF]
        if len(hits):
            rounds[rnd] = hits

    for j, (metric, xlabel) in enumerate(METRICS.items()):
        ax = axes[i, j]
        ref_vals = ref[metric].dropna()
        ax.axvspan(ref_vals.min(), ref_vals.max(), color=BAND_COLOR, alpha=BAND_ALPHA)
        ax.axvspan(ref_vals.quantile(0.25), ref_vals.quantile(0.75), color=BAND_COLOR, alpha=BAND_ALPHA + 0.15)
        ax.axvline(ref_vals.median(), color=REF_COLOR, lw=1.5, ls="--")
        rng = np.random.default_rng(42)
        ax.scatter(ref_vals, rng.uniform(-0.08, 0.08, len(ref_vals)), color=REF_COLOR, alpha=0.35, s=48, lw=0)

        groups = [(0, ref_vals, REF_COLOR)]
        for y, (rnd, hits) in enumerate(rounds.items(), start=1):
            vals, color = hits[metric].dropna(), ROUND_COLORS[target][y - 1]
            ax.axhline(y, color="#DDDDDD", lw=0.6)
            ax.scatter(vals, y + rng.uniform(-0.08, 0.08, len(vals)), color=color, alpha=0.85, s=80, lw=0, zorder=4)
            ax.scatter(vals.median(), y, color=color, s=176, marker="D", edgecolors="white", zorder=5)
            groups.append((y, vals, color))

        ax.set_ylim(-0.6, len(rounds) + 1.3)
        ax.set_yticks(range(len(rounds) + 1))
        ax.set_yticklabels(["Reference"] + [f"Round {r}" for r in rounds] if j == 0 else [])
        ax.set_xlabel(xlabel, fontsize=14)
        ax.grid(axis="x", ls="--", lw=0.4, alpha=0.5)
        ax.margins(x=0.06)
        if i == 0:
            ax.set_title(xlabel.split(" (")[0], fontsize=16, fontweight="bold")
        if j == 3:  # number of compounds per row, next to the last column
            trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
            for y, vals, color in groups:
                ax.text(1.02, y, f"n={len(vals)}", transform=trans, va="center", fontsize=11, color=color,
                        style="italic")
    axes[i, 0].annotate(target_label, xy=(-0.4, 0.5), xycoords="axes fraction", fontsize=18,
                        fontweight="bold", ha="center", va="center", rotation=90)

legend = [mpatches.Patch(color=BAND_COLOR, alpha=BAND_ALPHA + 0.15, label="Reference IQR"),
          mpatches.Patch(color=BAND_COLOR, alpha=BAND_ALPHA, label="Reference range (min–max)"),
          Line2D([0], [0], color=REF_COLOR, lw=1.5, ls="--", label="Reference median"),
          Line2D([0], [0], marker="D", color="gray", ls="none", ms=9, mec="white", label="Round median")]
fig.legend(handles=legend, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 0.995), frameon=False)

fig.savefig(OUT / "strip_plot_KD_MW_LE_LR.png", dpi=600, bbox_inches="tight")
fig.savefig(OUT / "strip_plot_KD_MW_LE_LR.pdf", bbox_inches="tight")
