import json
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import matplotlib.pyplot as plt
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import Descriptors, rdMolDescriptors
from multiprocessing import Pool, cpu_count

from utils.comp import calculate_isim

print("========== SAMPLE ANALYSIS ==========")


def scaffold_worker(smile):
    mol = Chem.MolFromSmiles(smile)

    if mol is None:
        return None
    
    return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))


def process_batch(smiles_batch):
    return [scaffold_worker(s) for s in smiles_batch]


def get_scaffolds(sample, outname, processors, batch_size):
    if processors is None:
        processors = cpu_count()

    print(f"\nUsing {processors} CPUs")
    print(f"Batch size: {batch_size}")

    pf = pq.ParquetFile(sample)

    writer = None
    total_processed = 0

    with Pool(processors) as pool:
        for batch in pf.iter_batches(batch_size=batch_size):
            
            smiles = batch.column("SMILES").to_pylist()

            chunk_size = len(smiles) // processors + 1
            chunks = [smiles[i:i + chunk_size]
                      for i in range(0, len(smiles), chunk_size)]

            results = pool.map(process_batch, chunks)

            scaffolds = []
            for res in results:
                scaffolds.extend(res)

            scaffolds = [s if s is not None else "NA" for s in scaffolds]

            table = pa.Table.from_batches([batch])
            scaffold_array = pa.array(scaffolds)
            table = table.append_column("scaffolds", scaffold_array)

            if writer is None:
                writer = pq.ParquetWriter(outname, table.schema)

            writer.write_table(table)

            total_processed += len(smiles)
            print(f"Processed: {total_processed:,} molecules")

    if writer:
        writer.close()

    print("\n Finished writing parquet with scaffolds!")


def scaffold_uniqueness(file, batch_size):
    """ Calculates scaffold uniqueness fraction from the sampled molecules

    Args:
        file: Dataframe file belonging to the sampled molecules with computed scaffolds
        batch_size: Batch size to stream parquet file

    Returns:
        Metrics relative to scaffold uniqueness in our sampled set
    """
    pf = pq.ParquetFile(file)

    unique_scaffolds = set()
    total = 0

    for batch in pf.iter_batches(batch_size=batch_size, columns=["scaffolds"]):
        scaffolds = batch.column("scaffolds").to_pylist()

        total += len(scaffolds)

        for s in scaffolds:
            if s is not None and s != "NA":
                unique_scaffolds.add(s)

    unique_count = len(unique_scaffolds)
    frac = unique_count / total

    print("----- Scaffold Statistics -----")
    print(f"Total molecules: {total}")
    print(f"Unique scaffolds: {unique_count}")
    print(f"Scaffold uniqueness: {frac:.4f}")

    return unique_scaffolds


def filter_unique_scaffolds(file, output, unique_scaffolds, batch_size):
    pf = pq.ParquetFile(file)

    seen = set()
    writer = None

    for batch in pf.iter_batches(batch_size=batch_size):
        table = pa.Table.from_batches([batch])

        scaffolds = table.column("scaffolds").to_pylist()

        mask = []
        for s in scaffolds:
            if s not in seen and s is not None and s != "NA":
                seen.add(s)
                mask.append(True)
            else:
                mask.append(False)

        filtered_table = table.filter(pa.array(mask))

        if filtered_table.num_rows == 0:
            continue

        if writer is None:
            writer = pq.ParquetWriter(output, filtered_table.schema)

        writer.write_table(filtered_table)

    if writer:
        writer.close()


if __name__ == "__main__":
    import argparse

    print("\n========== PARSING ARGUMENTS ==========")
    
    parser = argparse.ArgumentParser(description="Script to perform stratified sampling on large DBs")
    
    parser.add_argument("--sample",
                        help="Sampled parquet file containing molecules",
                        required=True)
    parser.add_argument("--processors",
                        help="Numbers of processors to parallelise on",
                        type=int,
                        required=True)
    parser.add_argument("--batch_size",
                        help="Batch size to stream parquet file",
                        type=int,
                        required=True)
    
    args = parser.parse_args()

    print("\n========== PARAMETERS ==========")
    print(args)
    sample = args.sample
    processors = args.processors
    batch_size = args.batch_size
    outname = sample.replace(".parquet", "_refine.parquet")
   
    print("\n========== STEP 1: COMPUTE SCAFFOLDS ==========")
    get_scaffolds(sample, outname, processors, batch_size)

    print("\n========== STEP 2: SCAFFOLD UNIQUENESS ==========")
    unique_scaffolds = scaffold_uniqueness(outname, batch_size)
    
    print("\n========== STEP 3: FILTER UNIQUE SCAFFOLDS ==========")
    unique_out = outname.replace(".parquet", "_unique.parquet") 
    filter_unique_scaffolds(outname, unique_out, unique_scaffolds, batch_size) 
