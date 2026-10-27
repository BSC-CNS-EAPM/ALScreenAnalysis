"""Statistical comparison of the Glide gscore distributions between AL rounds.

For each target and bias, the initial random set and the two AL predictions are
compared pairwise with a Mann-Whitney U test, a Kolmogorov-Smirnov test and the
Cliff's delta effect size. P-values are corrected with Benjamini-Hochberg (FDR).

Outputs (results/): statistical_summary.csv, cliffs_delta_heatmap.png/.pdf
"""
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, ks_2samp
from statsmodels.stats.multitest import multipletests

# Define paths, targets and labels
DATA = Path(__file__).resolve().parents[1] / "data" / "ALScreen"
OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)

TARGETS = {"IL1b": ["unbias", "bias1"], "IL2": ["unbias", "bias1", "bias2"], "IL36g": ["unbias", "bias1", "bias2"]}
TARGET_LABELS = {"IL1b": "IL-1β", "IL2": "IL-2", "IL36g": "IL-36γ"}
BIAS_LABELS = {"unbias": "Unbias", "bias1": "Bias 1", "bias2": "Bias 2"}
COMPARISONS = [("initial", "pred1"), ("pred1", "pred2"), ("initial", "pred2")]

def load_gscores(target, bias, rnd):
    gscores = pd.read_csv(DATA / target / bias / f"{rnd}.csv")["glide_gscore"].dropna()
    return gscores[gscores != 10000].values

# Compute the statistical comparisons
rows = []
for target, biases in TARGETS.items():
    for bias in biases:
        scores = {rnd: load_gscores(target, bias, rnd) for rnd in ["initial", "pred1", "pred2"]}
        for a, b in COMPARISONS:
            x, y = scores[a], scores[b]
            u, mw_p = mannwhitneyu(x, y, alternative="two-sided")
            rows.append({"Target": target, "Bias": bias, "Comparison": f"{a} vs {b}",
                         "Mean_1": x.mean(), "Mean_2": y.mean(),
                         "Median_1": np.median(x), "Median_2": np.median(y),
                         "MW_pvalue": mw_p, "KS_pvalue": ks_2samp(x, y).pvalue,
                         # Cliff's delta = P(x > y) - P(x < y), obtained from the U statistic
                         "Cliffs_delta": 2 * u / (len(x) * len(y)) - 1})

# Save results and print the summary table
stats = pd.DataFrame(rows)
stats["MW_qvalue"] = multipletests(stats["MW_pvalue"], method="fdr_bh")[1]
stats["KS_qvalue"] = multipletests(stats["KS_pvalue"], method="fdr_bh")[1]
stats["Significance"] = pd.cut(stats["MW_qvalue"], [0, 0.001, 0.01, 0.05, 1], labels=["***", "**", "*", "ns"],
                               include_lowest=True)
stats.to_csv(OUT / "statistical_summary.csv", index=False)
print(stats[["Target", "Bias", "Comparison", "Cliffs_delta", "MW_qvalue", "Significance"]])

# Cliff's delta heatmap (rows: target and bias, columns: comparison)
order = [f"{a} vs {b}" for a, b in COMPARISONS]
delta = stats.pivot(index=["Target", "Bias"], columns="Comparison", values="Cliffs_delta")[order]
stars = stats.pivot(index=["Target", "Bias"], columns="Comparison", values="Significance")[order]
index = [(t, b) for t, biases in TARGETS.items() for b in biases]
delta, stars = delta.loc[index], stars.loc[index]
annot = delta.apply(lambda col: col.map("{:.3f}".format)) + "\n" + stars.astype(str)

fig, ax = plt.subplots(figsize=(9, 6))
sns.heatmap(delta, annot=annot, fmt="", cmap="coolwarm", vmin=0, vmax=1, linewidths=2, ax=ax,
            annot_kws={"fontsize": 11, "fontweight": "bold"},
            cbar_kws={"label": "Cliff's δ", "ticks": [0, 0.147, 0.330, 0.474, 1]})
ax.collections[0].colorbar.set_ticklabels(["0\nnegligible", "0.147\nsmall", "0.330\nmedium", "0.474\nlarge", "1.0"])
ax.set_xticklabels(["Initial vs Pred 1", "Pred 1 vs Pred 2", "Initial vs Pred 2"])
ax.set_yticklabels([f"{TARGET_LABELS[t]} · {BIAS_LABELS[b]}" for t, b in index], rotation=0)
ax.xaxis.tick_top()
ax.set_xlabel("")
ax.set_ylabel("")
for i in np.cumsum([len(b) for b in TARGETS.values()])[:-1]:
    ax.axhline(i, color="white", lw=6)
plt.tight_layout()

plt.savefig(OUT / "cliffs_delta_heatmap.png", dpi=600)
plt.savefig(OUT / "cliffs_delta_heatmap.pdf")
