"""Grid of Glide gscore histograms for every target (rows) and bias (columns).

Each panel overlays the initial random set and the two AL predictions, with the
docked reference compounds as dashed lines and an inset zooming into the
best-scoring tail.

Outputs (results/): score_histograms_grid.png/.pdf
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# Define paths, targets, labels, colors and plotting parameters
DATA = Path(__file__).resolve().parents[1] / "data" / "ALScreen"
OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)

TARGETS = {"IL1b": ["unbias", "bias1"], "IL2": ["unbias", "bias1", "bias2"], "IL36g": ["unbias", "bias1", "bias2"]}
TARGET_LABELS = {"IL1b": "IL-1β", "IL2": "IL-2", "IL36g": "IL-36γ"}
BIASES = {"unbias": "Unbias", "bias1": "Bias 1", "bias2": "Bias 2"}
ROUNDS = {"initial": ("#999999", "Initial 50K"), "pred1": ("#4878CF", "Prediction 1"), "pred2": ("#D65F5F", "Prediction 2")}
REF_COLOR = "#B22222"
INSET_XLIM = (-10, -6)
N_BINS = 100

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "xtick.labelsize": 7.5, "ytick.labelsize": 7.5})

# Helper
def load_gscores(path):
    gscores = pd.read_csv(path)["glide_gscore"].dropna()
    return gscores[gscores != 10000]

# Create Figure
fig, axes = plt.subplots(3, 3, figsize=(11, 9.6), gridspec_kw={"hspace": 0.6, "wspace": 0.42})
ymax = 0
for i, (target, biases) in enumerate(TARGETS.items()):
    for j, bias in enumerate(BIASES):
        ax = axes[i, j]
        if bias not in biases:
            ax.axis("off")
            continue
        folder = DATA / target / bias
        scores = {rnd: load_gscores(folder / f"{rnd}.csv") for rnd in ROUNDS}
        refs = load_gscores(folder / "reference.csv")
        all_scores = pd.concat(scores.values())
        bins = np.linspace(all_scores.min(), all_scores.max(), N_BINS)

        # This is the zoomed-in panel
        axins = ax.inset_axes([0.55, 0.55, 0.6, 0.6])
        for z, (rnd, (color, _)) in enumerate(ROUNDS.items()):
            counts, _, _ = ax.hist(scores[rnd], bins=bins, color=color, alpha=0.75, zorder=2 + z)
            axins.hist(scores[rnd], bins=np.linspace(*INSET_XLIM, N_BINS), color=color, alpha=0.75, zorder=2 + z)
            ymax = max(ymax, counts.max())
        for a in (ax, axins):
            for ref in refs:
                a.axvline(ref, color=REF_COLOR, ls="--", lw=0.8, alpha=0.8, zorder=5)
        axins.set_xlim(*INSET_XLIM)
        axins.tick_params(labelsize=5.5, length=2)
        _, connectors = ax.indicate_inset_zoom(axins, edgecolor="0.3", linestyle=":")
        for line in connectors:
            line.set_visible(False)

        # Formatting
        ax.set_xlabel("Glide gscore (kcal/mol)", fontsize=10)
        ax.grid(axis="y", ls="--", lw=0.4, alpha=0.4)
        if j == 0:
            ax.set_ylabel("Counts", fontsize=10)
            ax.annotate(TARGET_LABELS[target], xy=(-0.3, 0.5), xycoords="axes fraction", ha="right",
                        va="center", fontsize=12, fontweight="bold", rotation=90)
        if i == 0 or (i == 1 and bias == "bias2"):
            ax.set_title(BIASES[bias], fontsize=12, fontweight="bold", pad=35)

# Same y range in all panels so the distributions can be compared directly
for ax in axes.flat:
    if ax.axison:
        ax.set_ylim(0, ymax * 1.08)

# Add legend
handles = [mpatches.Patch(color=c, alpha=0.75, label=l) for c, l in ROUNDS.values()]
handles.append(Line2D([0], [0], color=REF_COLOR, ls="--", label="Reference compounds"))
axes[0, 2].legend(handles=handles, loc="center", fontsize=10, frameon=False, labelspacing=1.2)

# Save Figure
fig.savefig(OUT / "score_histograms_grid.png", dpi=600, bbox_inches="tight")
fig.savefig(OUT / "score_histograms_grid.pdf", bbox_inches="tight")
