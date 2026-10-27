"""Docking enrichment, chemotype diversity and scaffold diversity along the AL rounds.

Figure 3x3 grid (rows: metric, columns: target) with one line per bias condition.
Run compute_chemotypes.py first.
"""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib import ticker

# Define paths, targets, biases, rounds and metrics to plot
OUT = Path(__file__).resolve().parent / "results"
TARGETS = {"IL1b": "IL-1β", "IL2": "IL-2", "IL36g": "IL-36γ"}
BIASES = {"unbias": ("grey", "o", "Unbias"), "bias1": ("blue", "s", "Bias 1"), "bias2": ("red", "^", "Bias 2")}
ROUNDS = ["initial", "pred1", "pred2"]

METRICS = [("top1_glide", "Top-1 Glide Gscore\n(kcal/mol)", "Docking\nEnrichment", False, "%.1f"),
           ("chemotype_diversity", "Chemotype Diversity", "Chemotype\nDiversity", True, "%.2f"),
           ("scaffold_chemotype_ratio", "Chemotypes with a unique\nBemis-Murcko scaffold", "Scaffold\nDiversity", True, "%.2f")]

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "xtick.labelsize": 8, "ytick.labelsize": 8})

# Load data (previously computed with compute_chemotypes.py)
df = pd.read_csv(OUT / "summary_chemotypes.csv")
df["x"] = df["Iteration"].map({r: i for i, r in enumerate(ROUNDS)})

# Create figure
fig, axes = plt.subplots(3, 3, figsize=(9, 9), gridspec_kw={"hspace": 0.2, "wspace": 0.2})
for i, (metric, ylabel, row_title, sharey, fmt) in enumerate(METRICS):
    pad = (df[metric].max() - df[metric].min()) * 0.12
    for j, (target, target_label) in enumerate(TARGETS.items()):
        ax = axes[i, j]
        for bias, (color, marker, label) in BIASES.items():
            sub = df[(df["Target"] == target) & (df["Bias"] == bias)].sort_values("x")
            ax.plot(sub["x"], sub[metric], color=color, marker=marker, alpha=0.6, lw=1.8, ms=6.5, clip_on=False)

        ax.set_xlim(-0.25, 2.25)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["Initial", "Pred 1", "Pred 2"] if i == 2 else [])
        if i == 2:
            ax.set_xlabel("AL Cycle", fontsize=8)
        if sharey:
            ax.set_ylim(df[metric].min() - pad, df[metric].max() + pad)
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter(fmt))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4, prune="both"))
        ax.grid(axis="y", ls="--", lw=0.45, alpha=0.55)
        if j == 0:
            ax.set_ylabel(ylabel, fontsize=8)
        if i == 0:
            ax.set_title(target_label, fontsize=11, fontweight="bold")
        if j == 2:
            ax.annotate(row_title, xy=(1.1, 0.5), xycoords="axes fraction", va="center",
                        fontsize=8, color="#555555", rotation=270)

handles = [mlines.Line2D([], [], color=c, marker=m, alpha=0.6, lw=1.8, label=l) for c, m, l in BIASES.values()]
fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.975), frameon=False,
           title="Bias condition")

fig.savefig(OUT / "chemical_diversity.png", dpi=600, bbox_inches="tight")
fig.savefig(OUT / "chemical_diversity.pdf", bbox_inches="tight")
