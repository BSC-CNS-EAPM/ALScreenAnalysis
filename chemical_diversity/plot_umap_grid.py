"""UMAP of the chemical space explored in each AL round.
One panel per target and bias.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import umap
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

# Define paths, targets, biases and rounds to plot
DATA = Path(__file__).resolve().parents[1] / "data" / "ALScreen"
OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)

TARGETS = {"IL1b": ["unbias", "bias1"], "IL2": ["unbias", "bias1", "bias2"], "IL36g": ["unbias", "bias1", "bias2"]}
TARGET_LABELS = {"IL1b": "IL-1β", "IL2": "IL-2", "IL36g": "IL-36γ"}
BIAS_LABELS = {"unbias": "Unbias", "bias1": "Bias 1", "bias2": "Bias 2"}
ROUNDS = {"initial": ("grey", "Initial 50k"), "pred1": ("blue", "Prediction 1"), "pred2": ("red", "Prediction 2")}
N_SAMPLE = 50000
SEED = 42
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

# Functions to featurize molecules and compute UMAP embeddings
def featurize(df):
    """Morgan fingerprints as a numpy array, dropping invalid SMILES."""
    mols = df["SMILES"].astype(str).apply(Chem.MolFromSmiles)
    df = df[mols.notna()].reset_index(drop=True)
    return df, np.array([FPGEN.GetFingerprintAsNumPy(m) for m in mols.dropna()], dtype=np.uint8)

# All compounds are computed in the same UMAP space
compounds, references = [], []
for target, biases in TARGETS.items():
    for bias in biases:
        for rnd in ROUNDS:
            df = pd.read_csv(DATA / target / bias / f"{rnd}.csv")
            compounds.append(df.loc[df["glide_gscore"] != 10000, ["SMILES"]].assign(Target=target, Bias=bias, Iteration=rnd))
        references.append(pd.read_csv(DATA / target / bias / "reference.csv")[["SMILES"]].assign(Target=target, Bias=bias))

compounds = pd.concat(compounds, ignore_index=True)
sample = compounds.iloc[np.random.RandomState(SEED).choice(len(compounds), N_SAMPLE, replace=False)]
sample, fps = featurize(sample)
references, ref_fps = featurize(pd.concat(references, ignore_index=True))

reducer = umap.UMAP(n_neighbors=25, min_dist=0.1, metric="jaccard", random_state=SEED)
sample[["UMAP1", "UMAP2"]] = reducer.fit_transform(fps)
references[["UMAP1", "UMAP2"]] = reducer.transform(ref_fps)

# Altough they are in the same UMAP space, we plot one panel per target and bias to avoid overcrowding
fig, axes = plt.subplots(3, 3, figsize=(14, 12), sharex=True, sharey=True)
for i, (target, biases) in enumerate(TARGETS.items()):
    for j, bias in enumerate(BIAS_LABELS):
        ax = axes[i, j]
        if bias not in biases:
            ax.set_axis_off()
            continue
        panel = sample[(sample["Target"] == target) & (sample["Bias"] == bias)]
        sns.scatterplot(data=panel, x="UMAP1", y="UMAP2", hue="Iteration", hue_order=list(ROUNDS),
                        palette={r: c for r, (c, _) in ROUNDS.items()}, s=10, alpha=0.3, lw=0, marker=".",
                        legend=False, ax=ax)
        refs = references[(references["Target"] == target) & (references["Bias"] == bias)]
        ax.scatter(refs["UMAP1"], refs["UMAP2"], marker="*", color="black", s=50, lw=0)
        if i == 0 or (i == 1 and bias == "bias2"):
            ax.set_title(BIAS_LABELS[bias], fontsize=18, fontweight="bold")
        if j == 2 or (i == 0 and j == 1):
            ax.annotate(TARGET_LABELS[target], xy=(1.03, 0.5), xycoords="axes fraction", va="center",
                        fontsize=18, fontweight="bold", rotation=270)

handles = [mlines.Line2D([], [], color=c, marker="o", ls="None", label=l) for c, l in ROUNDS.values()]
handles.append(mlines.Line2D([], [], color="black", marker="*", ls="None", label="Reference compounds"))
axes[0, 2].legend(handles=handles, loc="center", frameon=False, prop={"size": 14})
plt.tight_layout()

plt.savefig(OUT / "umap_grid.png", dpi=600)
plt.savefig(OUT / "umap_grid.pdf")
