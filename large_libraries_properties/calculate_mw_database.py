#!/usr/bin/env python3
"""Compute the molecular weight of every molecule in a Parquet library.

MW is RDKit ``Descriptors.MolWt`` (average MW incl. implicit H)

Designed to run as a SLURM array: task ``--task_id`` of ``--num_tasks`` takes every
``num_tasks``-th file. Finished files are skipped, so resubmitting resumes.
"""

import argparse
import os
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from rdkit import Chem
from rdkit.Chem import Descriptors


def process_file(job):
    input_file, output_dir, smiles_column = job
    output_file = Path(output_dir) / f"{input_file.stem}_mw.npy"
    if output_file.exists():
        return input_file.name, "skip"

    chunks = []
    n_invalid = 0

    reader = pq.ParquetFile(input_file)
    for row_group_index in range(reader.num_row_groups):
        smiles = reader.read_row_group(row_group_index, columns=[smiles_column]).column(0).to_pylist()
        mws = np.empty(len(smiles), dtype=np.float32)
        n = 0
        for smi in smiles:
            mol = Chem.MolFromSmiles(smi) if smi else None
            if mol is None:
                n_invalid += 1
                continue
            mws[n] = Descriptors.MolWt(mol)
            n += 1
        if n:
            chunks.append(mws[:n])

    mws = np.concatenate(chunks) if chunks else np.empty(0, dtype=np.float32)

    # Write to a temp name first so a killed job never leaves a truncated "finished" file.
    tmp_file = output_file.with_name(output_file.name + ".tmp.npy")
    np.save(tmp_file, mws)
    os.replace(tmp_file, output_file)
    mean = f" mean={mws.mean():.1f}" if mws.size else ""
    return input_file.name, f"done valid={mws.size:,} invalid={n_invalid:,}{mean}"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input_dir", required=True, type=Path)
    parser.add_argument("--output_dir", required=True, type=Path)
    parser.add_argument("--smiles_column", default="SMILES")
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--task_id", type=int, default=0)
    parser.add_argument("--num_tasks", type=int, default=1)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(args.input_dir.glob("*.parquet"))[args.task_id::args.num_tasks]
    jobs = [
        (f, args.output_dir, args.smiles_column)
        for f in files
        if not (args.output_dir / f"{f.stem}_mw.npy").exists()
    ]
    print(f"Task {args.task_id}/{args.num_tasks}: {len(files)} files assigned, {len(jobs)} remaining", flush=True)
    if not jobs:
        return

    with Pool(min(args.workers, len(jobs))) as pool:
        for i, (name, status) in enumerate(pool.imap_unordered(process_file, jobs, chunksize=1), start=1):
            print(f"[{i}/{len(jobs)}] {name}: {status}", flush=True)


if __name__ == "__main__":
    main()
