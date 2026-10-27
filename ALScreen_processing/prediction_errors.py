"""Error of the ALScreen predictions against the Glide gscore of the docked compounds.

For every target, bias and AL prediction it computes the MAE, RMSE, Pearson and
Spearman correlation between the ALScreen predicted score (dockscore) and the Glide
gscore, and the fold improvement from prediction 1 to prediction 2.

Outputs (results/): prediction_errors.csv, prediction_errors.png/.pdf
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

# Define paths, targets and labels
DATA = Path(__file__).resolve().parents[1] / "data" / "ALScreen"
OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)

TARGETS = {"IL1b": ["unbias", "bias1"], "IL2": ["unbias", "bias1", "bias2"], "IL36g": ["unbias", "bias1", "bias2"]}
TARGET_LABELS = {"IL1b": "IL-1β", "IL2": "IL-2", "IL36g": "IL-36γ"}
BIAS_LABELS = {"unbias": "Unbias", "bias1": "Bias 1", "bias2": "Bias 2"}

# Compute the errors
rows = []
for target, biases in TARGETS.items():
    for bias in biases:
        for pred in ["pred1", "pred2"]:
            df = pd.read_csv(DATA / target / bias / f"{pred}_alscreen.csv").dropna(subset=["dockscore", "glide_gscore"])
            df = df[df["glide_gscore"] < 50]  # failed dockings are stored as 10000
            diff = df["glide_gscore"] - df["dockscore"]
            rows.append({"Target": target, "Bias": bias, "Prediction": pred, "n": len(df),
                         "MAE": diff.abs().mean(), "RMSE": np.sqrt((diff ** 2).mean()),
                         "Pearson": pearsonr(df["dockscore"], df["glide_gscore"])[0],
                         "Spearman": spearmanr(df["dockscore"], df["glide_gscore"])[0]})

errors = pd.DataFrame(rows)
wide = errors.pivot_table(index=["Target", "Bias"], columns="Prediction", values=["MAE", "RMSE"])
for metric in ["MAE", "RMSE"]:
    wide[(metric, "fold")] = wide[(metric, "pred1")] / wide[(metric, "pred2")]

# Compute fold improvement
folds = pd.DataFrame({f"{metric}_fold": wide[(metric, "fold")] for metric in ["MAE", "RMSE"]}).reset_index()
errors = errors.merge(folds, on=["Target", "Bias"])

# Save results
errors.to_csv(OUT / "prediction_errors.csv", index=False)

# Create Figure
# Grouped bar plot: one pair of bars (pred1, pred2) per target, metric and bias
bars = [(t, m, b) for t, biases in TARGETS.items() for m in ["MAE", "RMSE"] for b in biases]
x = np.arange(len(bars))
width = 0.35
fig, ax = plt.subplots(figsize=(18, 7))
for shift, pred, color, label in [(0, "pred1", "blue", "Prediction 1"), (width, "pred2", "red", "Prediction 2")]:
    ax.bar(x + shift, [wide.loc[(t, b), (m, pred)] for t, m, b in bars], width, color=color, alpha=0.5, label=label)
for xi, (t, m, b) in zip(x, bars):
    ax.text(xi + width, wide.loc[(t, b), (m, "pred2")] + 0.05, f"Fold: {wide.loc[(t, b), (m, 'fold')]:.2f}",
            ha="center", rotation=90, fontsize=14)
for xi in np.cumsum([2 * len(b) for b in TARGETS.values()])[:-1]:
    ax.axvline(xi - 0.5 + width / 2, ls="--", lw=2, color="gray", alpha=0.7)

ax.set_xticks(x + width / 2)
ax.set_xticklabels([f"{TARGET_LABELS[t]}\n{m}\n{BIAS_LABELS[b]}" for t, m, b in bars], fontsize=16)
ax.set_ylabel("Prediction Error", fontsize=16)
ax.tick_params(labelsize=16)
ax.legend(fontsize=16)
ax.grid(axis="y", alpha=0.4)
plt.tight_layout()

plt.savefig(OUT / "prediction_errors.png", dpi=300)
plt.savefig(OUT / "prediction_errors.pdf")
