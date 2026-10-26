#!/usr/bin/env python3
"""Compute the linearity ratio and branching descriptors of every molecule in a library.

linearity_ratio = graph diameter / (heavy atoms - 1), i.e. 1 for a fully linear chain.
For each input file (Parquet or csv) a report <stem>_linearity_report.<ext> is written
with one row per molecule. Finished files are skipped, so a killed job can be resubmitted.

    python calculate_linearity_ratios_database.py --input_dir library_parquet/ --output_dir linearity/ --workers 48
"""
import argparse
import os
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdmolops


def compute_branching_metrics(smiles):
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is not None:
        mol = rdmolops.RemoveHs(mol)
    if mol is None or mol.GetNumAtoms() == 0:
        return {"valid": False, "smiles": smiles}
    n_atoms = mol.GetNumAtoms()
    ring_info = mol.GetRingInfo()

    terminal_atoms = sum(atom.GetDegree() == 1 for atom in mol.GetAtoms())
    # Branch atoms: degree >= 3, excluding ring fusion atoms without substituents
    branch_atoms = sum(atom.GetDegree() >= 3 and (ring_info.NumAtomRings(atom.GetIdx()) == 0
                                                 or any(not b.IsInRing() for b in atom.GetBonds()))
                       for atom in mol.GetAtoms())
    diameter = int(rdmolops.GetDistanceMatrix(mol).max()) if n_atoms > 1 else 0

    return {"valid": True, "smiles": smiles, "num_heavy_atoms": n_atoms, "num_rings": ring_info.NumRings(),
            "num_branch_atoms": branch_atoms, "num_terminal_atoms": terminal_atoms, "graph_diameter": diameter,
            "linearity_ratio": round(diameter / (n_atoms - 1), 3) if n_atoms > 1 else 1.0,
            "branch_density": round(branch_atoms / n_atoms, 3)}


def classify_branching(row):
    if not row["valid"]:
        return "invalid"
    if row["num_heavy_atoms"] <= 2:
        return "trivial (too small to branch)"
    if row["num_branch_atoms"] == 0 and row["linearity_ratio"] >= 0.5:
        return "linear"
    if row["branch_density"] >= 0.50:
        return "highly branched"
    if row["branch_density"] >= 0.30:
        return "branched"
    return "slightly branched"


def output_path(input_file, output_dir, output_format):
    ext = input_file.suffix.lower() if output_format == "same" else f".{output_format}"
    return output_dir / f"{input_file.stem}_linearity_report{ext}"


def process_file(job):
    input_file, output_file, smiles_column, id_column = job
    columns = [smiles_column, id_column]
    df = pd.read_csv(input_file, usecols=columns) if input_file.suffix == ".csv" else pd.read_parquet(input_file, columns=columns)

    report = pd.DataFrame([{**compute_branching_metrics(smi), id_column: mol_id, smiles_column: smi}
                           for smi, mol_id in zip(df[smiles_column], df[id_column])])
    if not report.empty:
        report["classification"] = report.apply(classify_branching, axis=1)

    tmp_file = output_file.with_name("tmp_" + output_file.name)
    report.to_csv(tmp_file, index=False) if output_file.suffix == ".csv" else report.to_parquet(tmp_file, index=False)
    os.replace(tmp_file, output_file)
    return input_file.name


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input_dir", required=True, type=Path, help="directory with parquet/csv files")
    parser.add_argument("--output_dir", required=True, type=Path)
    parser.add_argument("--smiles_column", default="SMILES")
    parser.add_argument("--id_column", default="ID")
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--output_format", choices=["same", "csv", "parquet"], default="same")
    parser.add_argument("--overwrite", action="store_true", help="reprocess files that already have a report")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    files = sorted([*args.input_dir.glob("*.parquet"), *args.input_dir.glob("*.csv")])
    jobs = [(f, output_path(f, args.output_dir, args.output_format), args.smiles_column, args.id_column)
            for f in files]
    if not args.overwrite:
        jobs = [job for job in jobs if not job[1].exists()]
    print(f"{len(files)} files found, {len(jobs)} remaining", flush=True)
    if not jobs:
        return

    with Pool(min(args.workers, len(jobs))) as pool:
        for i, name in enumerate(pool.imap_unordered(process_file, jobs), start=1):
            print(f"[{i}/{len(jobs)}] {name} done", flush=True)


if __name__ == "__main__":
    main()
