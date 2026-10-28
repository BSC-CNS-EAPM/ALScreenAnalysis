import os
import subprocess
import numpy as np
import pyarrow.parquet as pq
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from multiprocessing import Pool, cpu_count

print("RUNNING:", os.path.abspath(__file__))

_FPGEN = None

def init_worker(radius, fpsize):
    global _FPGEN
    _FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=radius,
                                                       fpSize=fpsize)


def find_parquet_files(root_dir):
    print(f"[INFO] Scanning directory: {root_dir}")
    parquet_files = []

    for subdir, _, files in os.walk(root_dir):
        print(f"[DEBUG] Entering: {subdir} ({len(files)} files)")
        for f in files:
            if f.endswith(".parquet") and not f.endswith("_filtered.parquet"):
                full_path = os.path.join(subdir, f)
                print(f"[FOUND] Parquet: {full_path}")
                parquet_files.append(full_path)

    parquet_files.sort()
    print(f"[INFO] Total parquet files found: {len(parquet_files)}")
    return parquet_files


def fingerprint_worker(args):
    smile = args
    if not smile:
        print("[WARN] Received None SMILES")
        return None

    mol = Chem.MolFromSmiles(smile)
    if mol is None:
        print(f"[WARN] Invalid SMILES: {smile}")
        return None

    fp = _FPGEN.GetFingerprint(mol)
    
    return list(fp.GetOnBits())


def get_sparse_fps(infile, outfile, processors, batch_size, radius, fpsize):
    
    print("\n==============================")
    print(f"[START] Processing file: {infile}")
    print(f"[INFO] Output sparse file: {outfile}")
    print(f"[INFO] Batch size: {batch_size}")
    print(f"[INFO] Radius: {radius}, FP size: {fpsize}")

    if processors is None:
        processors = cpu_count()
    
    print(f"[INFO] Using {processors} processors")

    pf = pq.ParquetFile(infile)
    print(f"[INFO] Num row groups: {pf.num_row_groups}")
    print(f"[INFO] Schema: {pf.schema}")

    total_rows = 0
    written_docs = 0

    with open(outfile, "w") as f, Pool(processors, initializer=init_worker, initargs=(radius, fpsize)) as pool:
        for batch_id, batch in enumerate(pf.iter_batches(batch_size=batch_size, columns=["SMILES"]), start=1):
            print(f"\n[BATCH {batch_id}] Reading batch")
            smiles = batch.column(0).to_pylist()
            total_rows += len(smiles)
            print(f"[BATCH {batch_id}] Rows: {len(smiles)}")

            results = pool.map(fingerprint_worker, smiles)

            for onbits in results:
                if onbits is None:
                    continue

                doc_id = written_docs + 1
                for bit in onbits:
                    f.write(f"{doc_id} {bit} 1\n")
                written_docs += 1

            print(f"[BATCH {batch_id}] Valid fingerprints written so far: {written_docs}")

    print("\n[END] Finished processing file")
    print(f"[SUMMARY] Total rows processed: {total_rows}")
    print(f"[SUMMARY] Total fingerprints written: {written_docs}")
    print(f"[OUTPUT] Sparse file saved at: {outfile}")
    
    return outfile


def calc_tapnn(input_file, outname, root_tapnn, t_sim):

    tapnn_path = os.path.join(root_tapnn, "build/apss")

    print("\n==============================")
    print("[START] TAPNN execution")
    print(f"[INFO] Executable: {tapnn_path}")
    print(f"[INFO] Input: {input_file}")
    print(f"[INFO] Output: {outname}")
    print(f"[INFO] Threshold: {t_sim}")

    if not os.path.exists(tapnn_path):
        print("[ERROR] TAPNN binary not found!")
        raise FileNotFoundError(tapnn_path)
    
    cmd = [
        tapnn_path,
        "tapnn",
        input_file,
        outname,
        "-fmtRead=ijv",
        "-readZidx",
        "-t",
        str(t_sim),
        "-sim",
        "tan"
    ]
    
    print("REAL CMD:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    print("[END] TAPNN finished successfully")
    print(f"[OUTPUT] Similarity file: {outname}")


def plot_similarity_histogram(sim_file, outpng):

    print("\n==============================")
    print("[START] Histogram generation")
    print(f"[INFO] Input similarity file: {sim_file}")

    sims = []

    with open(sim_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            try:
                sims.append(float(parts[2]))
            except ValueError:
                continue
    
    if not sims:
        print("[WARN] No similarity values found, skipping plot")
        return None

    sims = np.asarray(sims, dtype=np.float32)
    
    plt.figure()
    plt.hist(sims, bins=100)
    plt.xlabel("Tanimoto similarity")
    plt.ylabel("Frequency")

    output_png = outpng or f"{sim_file}.png"
    plt.savefig(output_png, dpi=1000)
    plt.close()

    print(f"[OUTPUT] Histogram saved at: {output_png}")
    print("[END] Histogram generation")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Script that converts molecules to sparse fingerprints to calculate pairwise similarity")
    
    parser.add_argument("--samples",
                        help="Path to directory containing parquet files",
                        required=True)
    parser.add_argument("--outname",
                        help="Output name for pairwise similarity results",
                        required=True)
    parser.add_argument("--batch_size",
                        help="Batch size for parquet file optimised reading",
                        type=int,
                        default=50000)
    parser.add_argument("--radius",
                        help="Radius for the fps generator",
                        type=int,
                        default=4)
    parser.add_argument("--fpsize",
                        help="Size of the fingerprint vector",
                        type=int,
                        default=1024) 
    parser.add_argument("--root_tapnn",
                        help="Root directory for the tapnn module",
                        required=True)
    parser.add_argument("--t_sim",
                        help="Tanimoto similarity threshold",
                        type=float,
                        required=True)
    parser.add_argument("--processors",
                        help="Number of processors to parallelise on",
                        type=int,
                        default=None)
    parser.add_argument("--outpng",
                        help="Output name for the distribution histogram")
    parser.add_argument("--plot_hist",
                        help="Flag to generate similarity hist or not",
                        action="store_true")
    
    args = parser.parse_args()
    
    samples = args.samples
    outname = args.outname
    batch_size = args.batch_size
    radius = args.radius
    fpsize = args.fpsize
    root_tapnn = args.root_tapnn
    processors = args.processors
    t_sim = args.t_sim
    outpng = args.outpng
    plot_hist = args.plot_hist
    
    print("\n========== CONFIG ==========")
    print(f"[ARGS] samples: {args.samples}")
    print(f"[ARGS] outname: {args.outname}")
    print(f"[ARGS] batch_size: {args.batch_size}")
    print(f"[ARGS] radius: {args.radius}")
    print(f"[ARGS] fpsize: {args.fpsize}")
    print(f"[ARGS] root_tapnn: {args.root_tapnn}")
    print(f"[ARGS] t_sim: {args.t_sim}")
    print(f"[ARGS] processors: {args.processors}")
    print(f"[ARGS] outpng: {args.outpng}")
    print(f"[ARGS] plot_hist: {args.plot_hist}")
    print("============================\n")

    parquet_files = find_parquet_files(args.samples)
    print(f"\nFound {len(parquet_files)} parquet files")
    
    for i, infile in enumerate(parquet_files, 1):
        
        print(f"\n######## FILE {i}/{len(parquet_files)} ########")

        base = os.path.splitext(os.path.basename(infile))[0] 
        outdir = os.path.dirname(infile)
    
        sparse_file = os.path.join(outdir, f"{base}.txt")
        output_file = os.path.join(outdir, f"{base}_sim.txt")

        print(f"[INFO] Input file: {infile}")
        print(f"[INFO] Sparse output: {sparse_file}")
        print(f"[INFO] TAPNN output: {output_file}")

        # Step 1: fingerprints to sparse file
        get_sparse_fps(infile, sparse_file, processors,
                       batch_size, radius, fpsize)

        # Step 2: TAPNN similarity
        calc_tapnn(sparse_file, output_file, root_tapnn, t_sim)

        # Setp 3: Generate similarity histogram (optional)
        if plot_hist:
            plot_similarity_histogram(output_file, outpng)
