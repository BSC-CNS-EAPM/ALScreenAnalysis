#!/usr/bin/env python3
"""Molecular weight distribution of two libraries, with the reference compounds as dashed lines.

Bins the per-molecule weights written by calculate_mw_database.py. The files are read one
at a time and in chunks, so the full libraries never have to fit in memory.

Reference compounds with KD < 10 µM are read from ../reference_compounds/results/ (run
compute_properties.py there first).

    python plot_mw_database.py results/mw/enamine_5B results/mw/enamine_65B
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REFERENCES = HERE.parent / "reference_compounds" / "results"
TARGETS = [("IL1b", "IL-1β", "forestgreen"), ("IL2", "IL-2", "mediumvioletred"), ("IL36g", "IL-36γ", "teal")]
KD_CUTOFF_uM = 10
MW_MAX = 2000.0       # upper edge of the histogram; heavier molecules only count towards the total
CHUNK = 10_000_000    # molecules loaded at a time from each file


def read_mw(results_dir, bin_width):
    """Histogram every *_mw.npy in results_dir, one chunk at a time."""
    files = sorted(results_dir.glob("*_mw.npy"))
    if not files:
        raise FileNotFoundError(f"No MW files found in {results_dir}")
    edges = np.arange(0, MW_MAX + bin_width, bin_width)
    counts = np.zeros(len(edges) - 1, dtype=np.int64)
    n_valid, mw_sum = 0, 0.0
    for f in files:
        mws = np.load(f, mmap_mode="r")
        n_valid += mws.size
        for start in range(0, mws.size, CHUNK):
            chunk = np.asarray(mws[start:start + CHUNK], dtype=np.float64)
            counts += np.histogram(chunk, bins=edges)[0]
            mw_sum += float(chunk.sum())
    print(f"{results_dir}: {n_valid:,} molecules from {len(files)} files, mean MW {mw_sum / n_valid:.1f} Da")
    return counts, edges, n_valid


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("results_dirs", nargs=2, type=Path, help="output folders of calculate_mw_database.py")
parser.add_argument("--labels", nargs=2, default=["Enamine REAL 5B", "Enamine REAL 65B"])
parser.add_argument("--bin_width", type=float, default=5.0, help="bin width in Da")
parser.add_argument("--xlim", type=float, nargs=2, default=(100, 800))
parser.add_argument("--output", type=Path, default=HERE / "results" / "mw_distribution.png")
args = parser.parse_args()

fig, ax = plt.subplots(figsize=(8, 5))
for results_dir, label, color in zip(args.results_dirs, args.labels, ["steelblue", "darkorange"]):
    counts, edges, n_valid = read_mw(results_dir, args.bin_width)
    # Density over all valid molecules (including those above 2000 Da), so areas are comparable
    ax.stairs(counts / (n_valid * args.bin_width), edges,
              fill=True, alpha=0.45, color=color, label=f"{label} (n={n_valid:,})")

for target, label, color in TARGETS:
    ref = pd.read_csv(REFERENCES / f"{target}_reference_properties.csv")
    mws = ref.loc[ref["KD_uM"] < KD_CUTOFF_uM, "MW"]
    for mw in mws:
        ax.axvline(mw, color=color, ls="--", lw=1.2, alpha=0.9)
    ax.plot([], [], color=color, ls="--", lw=1.2, label=f"{label} reference (<10 µM, n={len(mws)})")

ax.set(title="Molecular-weight distributions", xlabel="Molecular weight (Da)", ylabel="Density", xlim=args.xlim)
ax.grid(axis="y", alpha=0.25)
ax.legend()
fig.tight_layout()
args.output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.output, dpi=300)
