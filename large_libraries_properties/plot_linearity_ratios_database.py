#!/usr/bin/env python3
"""Linearity ratio distribution of two libraries, with the reference compounds as dashed lines.

Reads the reports written by calculate_linearity_ratios_database.py chunk by chunk,
so the full libraries never need to fit in memory.

Reference compounds with KD < 10 µM are read from ../reference_compounds/results/ (run compute_properties.py there first).

    python plot_linearity_ratios_database.py results/linearity/enamine_5B results/linearity/enamine_65B
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
REFERENCES = HERE.parent / "reference_compounds" / "results"
TARGETS = [("IL1b", "IL-1β", "green"), ("IL2", "IL-2", "red"), ("IL36g", "IL-36γ", "teal")]
KD_CUTOFF_uM = 10


def linearity_histogram(results_dir, edges):
    """Histogram of the linearity ratios of every report in results_dir."""
    files = sorted(results_dir.glob("*_linearity_report.parquet"))
    if not files:
        raise FileNotFoundError(f"No linearity reports found in {results_dir}")
    counts = np.zeros(len(edges) - 1, dtype=np.int64)
    for f in files:
        reader = pq.ParquetFile(f)
        for i in range(reader.num_row_groups):
            ratios = reader.read_row_group(i, columns=["linearity_ratio"]).column(0).to_pandas().dropna()
            counts += np.histogram(ratios, bins=edges)[0]
    print(f"{results_dir}: {counts.sum():,} molecules from {len(files)} files")
    return counts


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("results_dirs", nargs=2, type=Path, help="output folders of calculate_linearity_ratios_database.py")
parser.add_argument("--labels", nargs=2, default=["Enamine REAL 5B", "Enamine REAL 65B"])
parser.add_argument("--bins", type=int, default=100)
parser.add_argument("--output", type=Path, default=HERE / "results" / "linearity_ratio_distribution.png")
args = parser.parse_args()

edges = np.linspace(0, 1, args.bins + 1)
fig, ax = plt.subplots(figsize=(8, 5))
for results_dir, label, color in zip(args.results_dirs, args.labels, ["steelblue", "darkorange"]):
    counts = linearity_histogram(results_dir, edges)
    ax.stairs(counts / (counts.sum() * np.diff(edges)), edges, fill=True, alpha=0.45, color=color,
              label=f"{label} (n={counts.sum():,})")

for target, label, color in TARGETS:
    ref = pd.read_csv(REFERENCES / f"{target}_reference_properties.csv")
    ratios = ref.loc[ref["KD_uM"] < KD_CUTOFF_uM, "linearity_ratio"]
    for ratio in ratios:
        ax.axvline(ratio, color=color, ls="--", lw=1.2, alpha=0.9)
    ax.plot([], [], color=color, ls="--", lw=1.2, label=f"{label} reference (<10 µM, n={len(ratios)})")

ax.set(title="Linearity-ratio distributions", xlabel="Linearity ratio", ylabel="Density", xlim=(0, 1))
ax.grid(axis="y", alpha=0.25)
ax.legend()
fig.tight_layout()
args.output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.output, dpi=300)
