"""Cluster the interaction fingerprints of one target and plot the interaction frequencies per cluster.

Compounds are clustered with DBSCAN on the Tanimoto distance between their interaction
fingerprints, loaded from the pickle written by compute_fingerprints.py. For each cluster,
a heatmap shows the fraction of compounds making each interaction with each residue.

    python analyse_fingerprints.py --target IL1b --ifp results/IL1b_ifp.pkl

"""
import argparse
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN

OUT = Path(__file__).resolve().parent / "results"

# DBSCAN parameters (eps in Tanimoto distance) used in the paper
DBSCAN_PARAMS = {"IL1b": (0.14, 75), "IL2": (0.16, 25), "IL36g": (0.15, 100)}


def tanimoto_matrix(fps):
    """Pairwise Tanimoto similarity between the rows of a boolean matrix."""
    fps = fps.astype(np.float32)
    common = fps @ fps.T
    n_bits = fps.sum(axis=1)
    with np.errstate(invalid="ignore"):
        sim = common / (n_bits[:, None] + n_bits[None, :] - common)
    return np.nan_to_num(sim)


def plot_interaction_heatmap(df, outname):
    """Heatmap of interaction frequencies (residues x interaction types)."""
    freq = df.mean(axis=0).unstack(level="interaction").fillna(0)
    freq = freq.loc[(freq != 0).any(axis=1), (freq != 0).any(axis=0)]
    plt.figure(figsize=(10, 8), dpi=500)
    sns.heatmap(freq, annot=freq.round(4), cmap="Blues")
    plt.xticks(rotation=90, fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("Interactions", fontsize=14)
    plt.ylabel("Protein residue", fontsize=14)
    plt.tight_layout()
    plt.savefig(outname)
    plt.close()


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--target", required=True, help="target name, used to name the output files")
parser.add_argument("--ifp", required=True, type=Path, help="fingerprint pickle from compute_fingerprints.py")
parser.add_argument("--eps", type=float, help="DBSCAN eps in Tanimoto distance")
parser.add_argument("--min_samples", type=int, help="DBSCAN min_samples")
parser.add_argument("--output_dir", type=Path, default=OUT)
args = parser.parse_args()

if not args.ifp.is_file():
    raise FileNotFoundError(f"{args.ifp} does not exist; run compute_fingerprints.py first")

if args.eps is None or args.min_samples is None:
    if args.target not in DBSCAN_PARAMS:
        raise ValueError(f"No DBSCAN parameters stored for {args.target}; pass --eps and --min_samples")
    eps, min_samples = DBSCAN_PARAMS[args.target]
    eps = args.eps if args.eps is not None else eps
    min_samples = args.min_samples if args.min_samples is not None else min_samples
else:
    eps, min_samples = args.eps, args.min_samples

args.output_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_pickle(args.ifp).droplevel("ligand", axis=1)
distances = 1 - tanimoto_matrix(df.values)
labels = DBSCAN(metric="precomputed", eps=eps, min_samples=min_samples).fit(distances).labels_
print(f"{args.target} (eps={eps}, min_samples={min_samples}): {dict(Counter(labels))}")

pd.DataFrame({"ID": df.index, "cluster": labels}).to_csv(args.output_dir / f"{args.target}_clusters.csv", index=False)
for cluster in sorted(set(labels) - {-1}):
    plot_interaction_heatmap(df[labels == cluster],
                             args.output_dir / f"{args.target}_cluster{cluster}_int_perc_heatmap.pdf")
