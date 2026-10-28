import os
import json
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import matplotlib.pyplot as plt
from pathlib import Path
from itertools import product
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import Descriptors, rdMolDescriptors
from multiprocessing import Pool, cpu_count

from utils.comp import calculate_isim

print("========== SAMPLE BIN DIVISION ==========")


def assign_quantile_bin(value, q1, q2, q3):
    if value <= q1:
        return "q1"
    elif value <= q2:
        return "q2"
    elif value <= q3:
        return "q3"
    else:
        return "q4"


def generate_bins(base_dir="bins"):
    props = ["MW", "LogP", "TPSA"]
    quantiles = ["q1", "q2", "q3", "q4"]

    os.makedirs(base_dir, exist_ok=True)

    bin_names = []
    for combo in product(quantiles, repeat=3):
        name = f"MW{combo[0]}_LogP{combo[1]}_TPSA{combo[2]}"
        path = os.path.join(base_dir, name)
        os.makedirs(path, exist_ok=True)
        bin_names.append(name)

    return bin_names


def property_worker(smile):
    
    mol = Chem.MolFromSmiles(smile)

    if mol is None:
        return None
    
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)

    return mw, logp, tpsa


def process_batch(smiles_batch):
    return [property_worker(s) for s in smiles_batch]


def get_properties(sample, outname, processors, batch_size, bins_dir="bins"):

    print("\n================ START PIPELINE ================")
    print(f"Input file: {sample}")
    print(f"Processors: {processors}")
    print(f"Batch size: {batch_size}")
    print(f"Bins dir: {bins_dir}")

    if processors is None:
        processors = cpu_count()

    pf = pq.ParquetFile(sample)
    print("\nParquet metadata:")
    print(f"Num row groups: {pf.num_row_groups}")
    print(f"Schema: {pf.schema}")

    print("\n=== PASS 1: COMPUTING PROPERTIES FOR QUANTILES ===")

    mws, logps, tpsas = [], [], []

    with Pool(processors) as pool:
        for batch in pf.iter_batches(batch_size=batch_size):

            smiles = batch.column("SMILES").to_pylist()
            print(f"Batch size: {len(smiles)}")

            chunk_size = len(smiles) // processors + 1
            chunks = [smiles[i:i + chunk_size]
                      for i in range(0, len(smiles), chunk_size)]
            print(f"Total chunks: {len(chunks)}")

            results = pool.map(process_batch, chunks)
            print("Collected results from workers")

            for chunk in results:
                for res in chunk:
                    if res is None:
                        continue
                    mw, logp, tpsa = res
                    mws.append(mw)
                    logps.append(logp)
                    tpsas.append(tpsa)

    print("\nFinished PASS 1")
    
    print("\n========== COMPUTING QUANTILES ==========")
    q1_mw, q2_mw, q3_mw = np.quantile(mws, [0.25, 0.5, 0.75])
    q1_logp, q2_logp, q3_logp = np.quantile(logps, [0.25, 0.5, 0.75])
    q1_tpsa, q2_tpsa, q3_tpsa = np.quantile(tpsas, [0.25, 0.5, 0.75])

    print("\nQuantiles computed!")
    print("MW quantiles:", q1_mw, q2_mw, q3_mw)
    print("LogP quantiles:", q1_logp, q2_logp, q3_logp)
    print("TPSA quantiles:", q1_tpsa, q2_tpsa, q3_tpsa)

    print("\n=== PASS 2: BINNING MOLECULES ===")

    pf = pq.ParquetFile(sample)

    writers = {}
    total_written = 0

    with Pool(processors) as pool:

        for batch in pf.iter_batches(batch_size=batch_size):

            batch_df = batch.to_pandas()
            
            smiles = batch_df["SMILES"].tolist()

            chunk_size = len(smiles) // processors + 1

            chunks = [
                smiles[i:i + chunk_size]
                for i in range(0, len(smiles), chunk_size)]

            results = pool.map(process_batch, chunks)
            flat_results = [res for chunk in results for res in chunk]

            mw_vals = []
            logp_vals = []
            tpsa_vals = []

            for res in flat_results:

                if res is None:

                    mw_vals.append(None)
                    logp_vals.append(None)
                    tpsa_vals.append(None)

                else:

                    mw, logp, tpsa = res

                    mw_vals.append(mw)
                    logp_vals.append(logp)
                    tpsa_vals.append(tpsa)

            batch_df["MW"] = mw_vals
            batch_df["LogP"] = logp_vals
            batch_df["TPSA"] = tpsa_vals

            batch_df["MW_bin"] = batch_df["MW"].apply(
                lambda x: assign_quantile_bin(x, q1_mw, q2_mw, q3_mw)
            )

            batch_df["LogP_bin"] = batch_df["LogP"].apply(
                lambda x: assign_quantile_bin(x, q1_logp, q2_logp, q3_logp)
            )

            batch_df["TPSA_bin"] = batch_df["TPSA"].apply(
                lambda x: assign_quantile_bin(x, q1_tpsa, q2_tpsa, q3_tpsa)
            )

            batch_df["bin_name"] = (
                "MW" + batch_df["MW_bin"]
                + "_LogP" + batch_df["LogP_bin"]
                + "_TPSA" + batch_df["TPSA_bin"])

            grouped = batch_df.groupby("bin_name")

            for bin_name, sub_df in grouped:
                sub_df = sub_df.drop(columns=["MW_bin", "LogP_bin", "TPSA_bin", "bin_name"])

                table = pa.Table.from_pandas(sub_df)
                outpath = os.path.join(bins_dir, bin_name, "data.parquet")

                if bin_name not in writers:
                    writers[bin_name] = pq.ParquetWriter(outpath, table.schema)

                writers[bin_name].write_table(table)

                total_written += len(sub_df)

            print(f"Total rows written so far: {total_written}")

    for writer in writers.values():
        writer.close()

    print("\n========== DONE ==========")
    print(f"Total rows written: {total_written}")


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
    parser.add_argument("--bins_dir",
                        help="Directory to divide the sampled molecules in",
                        default="bins")
    
    args = parser.parse_args()

    print("\n========== PARAMETERS ==========")
    print(args)
    sample = args.sample
    processors = args.processors
    batch_size = args.batch_size
    bins_dir = args.bins_dir
    outname = sample.replace(".parquet", "_prop.parquet")
   
    print("\n========== STEP 1: COMPUTE PROPERTY BINS ==========")
    generate_bins(base_dir=bins_dir)

    print("\n========== STEP 2: COMPUTE PROPERTIES ==========")
    get_properties(sample, outname, processors, batch_size, bins_dir) 
