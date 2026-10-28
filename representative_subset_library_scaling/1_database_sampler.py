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

print("========== DATABASE SAMPLER START ==========")

rng = np.random.default_rng(1234)

BASE_DB = {
    "001": "Enamine_REAL_5B",
}


def list_hac_files(mol_dir, hac, extension="*.parquet"):
    """ Retrieve list of files within the specified HAC directory.

    Args:
        mol_dir: Root path to the processed Composite DB
        hac: HAC specification to access only that specified directory
        extension: File type to extract

    Return:
        List of files from the specified HAC
    """

    print(f"\n[LIST] Looking for HAC directory: HAC_{hac}")

    hac_dir = Path(mol_dir) / f"HAC_{hac}"

    if not hac_dir.exists():
        print(f"[WARNING] HAC directory does not exist: {hac_dir}")
        return []

    return list(hac_dir.glob(extension))


def get_sample_quota(met_file, size, hac_min, hac_max, frac, min_sample, uniform, u_lambda):
    """ Given a metrics file with the counts per HAC and per DB calculate the quota.

    Args:
        met_file: Metrics files containing count information
        size: Target size to sample
        hac_min: Minimum HAC to sample from
        hac_max: Maximum HAC to sample from
        frac: Minimum representativity percentage per DB
        min_sample: Minimum number to sample per DB
        uniform: Flag to apply uniformity to the weights
        u_lambda: Weight normalisation factor

    Returns:
        Quota dictionary organised as:
        {hac: {db_id: target_samples}}
    """
    print("\n========== COMPUTING SAMPLE QUOTA ==========")
    print(f"[QUOTA] Metrics file: {met_file}")

    df = pd.read_csv(met_file)

    print(f"[QUOTA] Metrics rows loaded: {len(df)}")

    pop = df.groupby(["db_id", "hac_avg"])["count"].sum().reset_index()

    print(f"[QUOTA] Unique DB x HAC combinations: {len(pop)}")

    pop = pop[(pop["hac_avg"] >= hac_min) &
              (pop["hac_avg"] <= hac_max)]

    print(f"[QUOTA] Rows after HAC filter [{hac_min},{hac_max}]: {len(pop)}")

    total = pop["count"].sum()

    print(f"[QUOTA] Total molecules in filtered population: {total}")

    pop["combined_weight"] = pop["count"] / total

    if uniform:
        print(f"[QUOTA] Applying uniform smoothing (alpha={u_lambda})")
        pop["combined_weight"] = pop["combined_weight"] ** u_lambda
        pop["combined_weight"] /= pop["combined_weight"].sum()

    pop["expected"] = pop["combined_weight"] * size
    pop["target_n"] = np.floor(pop["expected"]).astype(int)

    print("[QUOTA] Initial quota assigned")

    if frac is not None:
        print(f"[QUOTA] Enforcing minimum fraction per DB: {frac}")
        min_frac_samples = int(size * frac)
        pop["target_n"] = np.maximum(pop["target_n"], min_frac_samples)

    if min_sample is not None:
        print(f"[QUOTA] Enforcing minimum sample per DB: {min_sample}")
        pop["target_n"] = np.maximum(pop["target_n"], min_sample)

    # Target sample number can't be larger than the DB size
    pop["target_n"] = np.minimum(pop["target_n"], pop["count"])

    remaining = size - pop["target_n"].sum()

    print(f"[QUOTA] Remaining molecules after rounding: {remaining}")

    if remaining > 0:

        print("[QUOTA] Distributing remaining samples")

        remainders = (pop["expected"] - np.floor(pop["expected"]))
        order = remainders.sort_values(ascending=False).index

        for i in order[:remaining]:
            if pop.loc[i, "target_n"] < pop.loc[i, "count"]:
                pop.loc[i, "target_n"] += 1

    elif remaining < 0:

        print("[QUOTA] Reducing excess samples")

        order = (pop["expected"] - np.floor(pop["expected"])).sort_values().index

        for i in order:

            if remaining == 0:
                break

            if min_sample is None or pop.loc[i, "target_n"] > min_sample:
                pop.loc[i, "target_n"] -= 1
                remaining += 1

    print("[QUOTA] Final quota assigned")

    pop.to_csv(f"population_weight_{u_lambda}.csv", index=False)
    print("[QUOTA] Population weights saved")

    quota = {}

    for row in pop.itertuples():

        if row.target_n <= 0:
            continue

        hac = int(row.hac_avg)
        db = int(row.db_id)

        quota.setdefault(hac, {})
        quota[hac][db] = int(row.target_n)

    with open(f'sample_quota_{size}.json', 'w') as f:
        json.dump(quota, f)
    print(f"[QUOTA] Quota saved to sample_quota_{size}.json")

    total_planned = sum(sum(db.values()) for db in quota.values())
    print(f"[QUOTA] TOTAL PLANNED SAMPLE SIZE: {total_planned}")

    return quota


def iter_sample_hac(mol_dir, hac, quota, chunk_size=100000, seed=1234):
    """Sample molecules for a specific HAC according to DB quotas.

    Args:
        mol_dir: root directory containing HAC folders
        hac: HAC value being sampled
        quota: dictionary {db_id: target_samples}
        chunk_size: number of rows per batch read from parquet
        seed: RNG seed

    Yields:
        Sampled molecule DataFrames, one per input batch.  This avoids holding
        all selected molecules for a HAC in memory.
    """
    print(f"\n========== SAMPLING HAC {hac} ==========")

    rng = np.random.default_rng(seed)

    files = list_hac_files(mol_dir, hac)

    if len(files) == 0:
        print(f"[HAC {hac}] No files found")
        return

    rng.shuffle(files)

    print(f"[HAC {hac}] Files shuffled")

    remaining = quota.copy()

    print(f"[HAC {hac}] Initial quotas: {remaining}")

    sampled_count = 0

    file_counter = 0
    batch_counter = 0

    for f in files:

        if all(v <= 0 for v in remaining.values()):
            print(f"[HAC {hac}] Quota satisfied early")
            break

        file_counter += 1
        print(f"[HAC {hac}] Processing file {file_counter}/{len(files)}: {f}")

        # --- FILE OPEN PROTECTION ---
        try:
            pf = pq.ParquetFile(f)
        except Exception as e:
            print(f"[ERROR] Cannot open file {f}: {e}")
            continue

        # --- ITERATOR PROTECTION (CRITICAL FIX) ---
        try:
            for batch_idx, batch in enumerate(pf.iter_batches(batch_size=chunk_size)):

                batch_counter += 1
                print(f"[HAC {hac}] Batch {batch_counter}")

                # --- BATCH CONVERSION PROTECTION ---
                try:
                    df = batch.to_pandas()
                    print(f"[HAC {hac}] Batch size: {len(df)}")
                except Exception as e:
                    print(f"[ERROR] Failed converting batch {batch_idx} in {f}: {e}")
                    break

                if all(v <= 0 for v in remaining.values()):
                    print(f"[HAC {hac}] Quota satisfied before reading next batch")
                    break

                df["db_id"] = df["db_id"].astype(int)

                groups = df.groupby("db_id")
                chosen_idx = []

                for db_id, group in groups:

                    if db_id not in remaining:
                        continue

                    needed = remaining[db_id]

                    if needed <= 0:
                        continue

                    indices = group.index.to_numpy()
                    take = min(len(indices), needed)
                    sampled = rng.choice(indices, take, replace=False)
                    chosen_idx.append(sampled)
                    remaining[db_id] -= take

                    print(f"[HAC {hac}] DB {db_id} sampled {take}, remaining {remaining[db_id]}")

                    if all(v <= 0 for v in remaining.values()):
                        print(f"[HAC {hac}] Quota satisfied inside batch")
                        break

                if chosen_idx:
                    selected = df.loc[np.concatenate(chosen_idx)].copy()
                    selected["hac_avg"] = hac
                    print(f"[HAC {hac}] Selected {len(selected)} molecules")
                    sampled_count += len(selected)
                    # The caller writes this before the next batch is decoded.
                    yield selected

                if all(v <= 0 for v in remaining.values()):
                    print(f"[HAC {hac}] Quota satisfied")
                    break

        except OSError as e:
            print(f"[ERROR] I/O error while reading {f}: {e}")
            continue

    if sampled_count == 0:
        print(f"[HAC {hac}] No molecules sampled")
        return

    print(f"[HAC {hac}] Total sampled molecules: {sampled_count}")


MAX_ROWS_PER_PARQUET = 20_000_000


def execute_sampling(mol_dir, quota, outname, seed=1234,
                     max_rows_per_file=MAX_ROWS_PER_PARQUET,
                     chunk_size=100000):
    """Execute the sampling plan and produce the final dataset.

    Samples over ``max_rows_per_file`` are written to numbered parquet files in
    a directory named after ``outname`` without its suffix.
    """
    print("\n========== EXECUTING SAMPLING ==========")

    outpath = Path(outname)
    planned_size = sum(sum(db.values()) for db in quota.values())
    partitioned_output = planned_size > max_rows_per_file
    output_dir = outpath.with_suffix("") if partitioned_output else None

    if partitioned_output:
        if output_dir.exists() and any(output_dir.glob("part-*.parquet")):
            raise FileExistsError(
                f"Refusing to mix a new sample with existing parquet parts in: {output_dir}"
            )
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"[EXECUTE] Writing parquet parts of up to {max_rows_per_file:,} rows to: {output_dir}")
    else:
        outpath.parent.mkdir(parents=True, exist_ok=True)

    writer = None
    rows_in_current_file = 0
    part_number = 0
    total_written = 0

    def open_writer(table):
        """Open the next output parquet writer using the table schema."""
        nonlocal writer, rows_in_current_file, part_number

        if partitioned_output:
            part_number += 1
            part_path = output_dir / f"part-{part_number:05d}.parquet"
        else:
            part_path = outpath

        print(f"[EXECUTE] Creating Parquet writer: {part_path}")
        writer = pq.ParquetWriter(part_path, table.schema)
        rows_in_current_file = 0

    def write_dataframe(df):
        """Write one selected batch, rolling over at a part boundary."""
        nonlocal writer, rows_in_current_file, total_written

        start = 0
        while start < len(df):
            if writer is None:
                schema_table = pa.Table.from_pandas(
                    df.iloc[start:start + 1], preserve_index=False
                )
                open_writer(schema_table)

            available = max_rows_per_file - rows_in_current_file
            stop = min(start + available, len(df))
            table = pa.Table.from_pandas(df.iloc[start:stop], preserve_index=False)
            writer.write_table(table)

            written = stop - start
            rows_in_current_file += written
            total_written += written
            start = stop
            print(f"[EXECUTE] Written {written} rows (TOTAL {total_written})")

            if partitioned_output and rows_in_current_file == max_rows_per_file:
                writer.close()
                writer = None

    for hac, db_quota in quota.items():
        hac = int(hac)
        print(f"\n[EXECUTE] Starting HAC {hac}")

        hac_written = 0
        for df in iter_sample_hac(
            mol_dir, hac, db_quota, chunk_size=chunk_size, seed=seed + hac
        ):
            hac_written += len(df)
            write_dataframe(df)

        if hac_written == 0:
            print(f"[EXECUTE] HAC {hac} returned empty sample")
        else:
            print(f"[EXECUTE] HAC {hac} sampled {hac_written} molecules")

    if total_written == 0:
        raise RuntimeError("No molecules sampled")

    if writer is not None:
        writer.close()

    print(f"\n[EXECUTE] FINAL TOTAL WRITTEN: {total_written}")
    print(f"[EXECUTE] Sample saved to: {output_dir if partitioned_output else outpath}")

    return total_written


def db_display_name(db_id):
    """Return a readable DB label from numeric or string DB identifiers."""
    try:
        key = f"{int(db_id):03d}"
    except (TypeError, ValueError):
        key = str(db_id)

    return BASE_DB.get(key, f"DB_{key}")


def plot_hac_db_histogram(sample_file, plot_out, hac_min=None, hac_max=None):
    """Plot the sampled HAC distribution as a stacked histogram by database.

    Args:
        sample_file: Parquet file produced by execute_sampling
        plot_out: Output image path
        hac_min: Optional minimum HAC shown on the x-axis
        hac_max: Optional maximum HAC shown on the x-axis
    """
    print("\n========== PLOTTING HAC HISTOGRAM ==========")
    sample_path = Path(sample_file)
    files = [sample_path] if sample_path.is_file() else sorted(sample_path.glob("*.parquet"))
    if not files:
        print("[PLOT] No parquet files found; skipping plot")
        return

    # Reading all rows with pandas would OOM for multi-billion-row samples.
    # Keep only the small HAC x DB aggregate while scanning parquet batches.
    aggregate = {}
    observed_hacs = []
    for file_number, parquet_path in enumerate(files, start=1):
        print(f"[PLOT] Scanning {file_number}/{len(files)}: {parquet_path}")
        parquet_file = pq.ParquetFile(parquet_path)
        for batch in parquet_file.iter_batches(
            batch_size=100000, columns=["db_id", "hac_avg"]
        ):
            batch_df = batch.to_pandas()
            batch_df["db_id"] = batch_df["db_id"].astype(int)
            batch_df["hac_avg"] = batch_df["hac_avg"].astype(int)
            observed_hacs.extend((batch_df["hac_avg"].min(), batch_df["hac_avg"].max()))
            for (batch_hac, db_id), count in batch_df.groupby(["hac_avg", "db_id"]).size().items():
                key = (int(batch_hac), int(db_id))
                aggregate[key] = aggregate.get(key, 0) + int(count)

    if not aggregate:
        print("[PLOT] Sample is empty; skipping plot")
        return

    if hac_min is None:
        hac_min = min(observed_hacs)
    if hac_max is None:
        hac_max = max(observed_hacs)

    records = [
        {"hac_avg": batch_hac, "db_id": db_id, "count": count}
        for (batch_hac, db_id), count in aggregate.items()
        if hac_min <= batch_hac <= hac_max
    ]
    counts = (
        pd.DataFrame(records)
        .pivot(index="hac_avg", columns="db_id", values="count")
        .fillna(0)
        .astype(np.int64)
        .sort_index()
    )

    all_hacs = np.arange(hac_min, hac_max + 1)
    counts = counts.reindex(all_hacs, fill_value=0)
    counts = counts.reindex(sorted(counts.columns), axis=1)
    counts.columns = [db_display_name(db_id) for db_id in counts.columns]

    print(f"[PLOT] HAC range: {hac_min}-{hac_max}")
    print(f"[PLOT] Databases in sample: {list(counts.columns)}")

    fig_width = max(10, min(28, 0.35 * len(all_hacs)))
    fig, ax = plt.subplots(figsize=(fig_width, 7))

    colors = plt.get_cmap("tab20").colors
    bottom = np.zeros(len(counts), dtype=int)

    for i, column in enumerate(counts.columns):
        values = counts[column].to_numpy()
        ax.bar(
            counts.index,
            values,
            bottom=bottom,
            width=0.9,
            label=column,
            color=colors[i % len(colors)],
            edgecolor="white",
            linewidth=0.2,
        )
        bottom += values

    ax.set_xlabel("HAC")
    ax.set_ylabel("Number of sampled molecules")
    ax.set_title("Sampled HAC distribution by database")
    ax.set_xlim(hac_min - 0.5, hac_max + 0.5)
    ax.set_xticks(all_hacs)
    ax.tick_params(axis="x", rotation=90)
    ax.legend(title="Database", bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.grid(axis="y", alpha=0.25)

    fig.tight_layout()
    fig.savefig(plot_out, dpi=300)
    plt.close(fig)

    print(f"[PLOT] Histogram saved to: {plot_out}")


if __name__ == "__main__":
    import argparse

    print("\n========== PARSING ARGUMENTS ==========")

    parser = argparse.ArgumentParser(description="Script to perform stratified sampling on large DBs")

    parser.add_argument("--met_file",
                        help="Path to file containing total DB analysis",
                        required=True)
    parser.add_argument("--mol_dir",
                        help="Path to directory containing the molecules",
                        required=True)
    parser.add_argument("--size",
                        help="Number of molecules to sample from the DB",
                        required=True,
                        type=int)
    parser.add_argument("--hac_min",
                        help="Define HAC minimum to sample from",
                        type=int,
                        default=4)
    parser.add_argument("--hac_max",
                        help="Define HAC maximum to sample from",
                        type=int,
                        default=80)
    parser.add_argument("--frac",
                        help="Minimum percentage of molecules to sample per DB",
                        type=float)
    parser.add_argument("--min_sample",
                        help="Minimum number of molecules to sample per DB",
                        type=int)
    parser.add_argument("--outname",
                        help="Output name for the parquet file containing the sampled set",
                        required=True)
    parser.add_argument("--uniform",
                        help="Apply uniform smoothing to DB-HAC weights",
                        action="store_true")
    parser.add_argument("--u_lambda",
                        help="Smoothing strength (0=uniform, 1=original)",
                        type=float)
    parser.add_argument("--seed",
                        help="Seed for reproducible results",
                        type=int,
                        default=1234)
    parser.add_argument("--chunk_size",
                        help="Input rows decoded at once; lower values reduce peak RAM",
                        type=int,
                        default=100000)
    parser.add_argument("--plot_hac_hist",
                        help="Plot a stacked HAC histogram coloured by source database",
                        action="store_true")
    parser.add_argument("--hac_hist_out",
                        help="Output image path for --plot_hac_hist",
                        default=None)

    args = parser.parse_args()

    print("\n========== PARAMETERS ==========")
    print(args)
    met_file = args.met_file
    mol_dir = args.mol_dir
    size = args.size
    hac_min = args.hac_min
    hac_max = args.hac_max
    frac = args.frac
    min_sample = args.min_sample
    outname = args.outname
    uniform = args.uniform
    u_lambda = args.u_lambda
    seed = args.seed
    chunk_size = args.chunk_size
    plot_hac_hist = args.plot_hac_hist
    hac_hist_out = args.hac_hist_out

    print("\n========== STEP 1: COMPUTE QUOTA ==========")
    # Sampling HACs
    quota = get_sample_quota(met_file, size, hac_min, hac_max, frac, min_sample, uniform, u_lambda)

    print("\n========== STEP 2: SAMPLE DATABASE ==========")
    total_written = execute_sampling(mol_dir, quota, outname, seed,
                                     chunk_size=chunk_size)

    if plot_hac_hist:
        if hac_hist_out is None:
            hac_hist_out = str(Path(outname).with_suffix("")) + "_hac_hist_by_db.png"
        print("\n========== STEP 3: PLOT HAC HISTOGRAM ==========")
        plot_hac_db_histogram(outname, hac_hist_out, hac_min, hac_max)

    print("\n========== SAMPLING FINISHED ==========")
    print(f"TOTAL MOLECULES WRITTEN: {total_written}")
