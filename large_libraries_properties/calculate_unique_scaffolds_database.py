#!/usr/bin/env python3
"""Count the unique Bemis-Murcko scaffolds of a molecular library stored as Parquet files.

Each partition is reduced to its set of unique scaffolds and saved as
<stem>_scaffolds.txt. Finished partitions are skipped, so a job can be resubmitted.
Once all partitions are done, the per-file sets are merged into unique_scaffolds.txt
and the total is written to summary.txt.

    python calculate_unique_scaffolds_database.py --input_dir library_parquet/ --output_dir scaffolds/ --workers 48
"""
import argparse
import os
from multiprocessing import Pool
from pathlib import Path

import pyarrow.parquet as pq
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def get_scaffold(smiles):
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return None
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(scaffold) if scaffold.GetNumAtoms() else None


def process_file(job):
    input_file, output_dir, smiles_column = job
    output_file = output_dir / f"{input_file.stem}_scaffolds.txt"
    smiles = pq.read_table(input_file, columns=[smiles_column]).column(0).to_pylist()
    scaffolds = {s for s in map(get_scaffold, smiles) if s}

    # Write to a temp name first so a killed job never leaves a truncated "finished" file
    tmp_file = output_file.with_suffix(".tmp")
    tmp_file.write_text("\n".join(sorted(scaffolds)) + "\n")
    os.replace(tmp_file, output_file)
    return input_file.name, len(scaffolds)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input_dir", required=True, type=Path)
    parser.add_argument("--output_dir", required=True, type=Path)
    parser.add_argument("--smiles_column", default="SMILES")
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    args = parser.parse_args()

    partial_dir = args.output_dir / "per_file"
    partial_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(args.input_dir.glob("*.parquet"))
    jobs = [(f, partial_dir, args.smiles_column) for f in files
            if not (partial_dir / f"{f.stem}_scaffolds.txt").exists()]
    print(f"{len(files)} files found, {len(jobs)} remaining", flush=True)

    if jobs:
        with Pool(min(args.workers, len(jobs))) as pool:
            for i, (name, n) in enumerate(pool.imap_unordered(process_file, jobs), start=1):
                print(f"[{i}/{len(jobs)}] {name}: {n:,} unique scaffolds", flush=True)

    all_scaffolds = set()
    for f in sorted(partial_dir.glob("*_scaffolds.txt")):
        all_scaffolds.update(f.read_text().split())
    (args.output_dir / "unique_scaffolds.txt").write_text("\n".join(sorted(all_scaffolds)) + "\n")
    (args.output_dir / "summary.txt").write_text(f"Total unique scaffolds: {len(all_scaffolds)}\n")
    print(f"Total unique scaffolds: {len(all_scaffolds):,}")


if __name__ == "__main__":
    main()
