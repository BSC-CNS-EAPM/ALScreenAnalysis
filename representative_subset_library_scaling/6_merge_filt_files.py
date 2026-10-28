import os
import math
import random
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def get_parquet_files(input_dir):

    print("\n========== STEP 1: DISCOVERING BINS ==========")

    parquet_files = []

    for subdir, _, files in os.walk(input_dir):
        for f in files:
            if f.endswith("_filtered.parquet"):

                filepath = os.path.join(subdir, f)
                parquet_files.append(filepath)
                print(f"[FOUND BIN FILE] {filepath}")

    if len(parquet_files) == 0:
        raise RuntimeError("No parquet files found")

    print(f"\n[TOTAL BINS FOUND] {len(parquet_files)}")

    return sorted(parquet_files)


def stratified_sample_by_database(table, quota, database_column, rng, bin_name, seed=1234):

    print(f"\n----- DATABASE STRATIFIED SAMPLING ({bin_name}) -----")

    df = table.to_pandas()
    total_rows = len(df)

    print(f"[BIN SIZE] {total_rows}")

    if quota >= total_rows:
        print("[QUOTA >= BIN SIZE] Taking entire bin")
        return pa.Table.from_pandas(df, preserve_index=False)

    grouped = df.groupby(database_column)
    db_counts = grouped.size()

    print("[DB DISTRIBUTION BEFORE SAMPLING]")

    for db, count in db_counts.items():
        print(f"  DB {db}: {count}")

    print("\n[ALLOCATING QUOTA PER DATABASE]")

    allocation_info = []

    #
    # First pass:
    # compute exact expectations
    #
    for db_name, group in grouped:

        proportion = len(group) / total_rows
        expected = proportion * quota

        base = math.floor(expected)
        remainder = expected - base

        allocation_info.append({
            "db": db_name,
            "group": group,
            "expected": expected,
            "base": base,
            "remainder": remainder
        })

        print(
            f"  DB {db_name} | "
            f"proportion={proportion:.4f} | "
            f"expected={expected:.2f} | "
            f"base={base} | "
            f"remainder={remainder:.4f}"
        )

    #
    # Largest remainder method
    #
    assigned_total = sum(x["base"] for x in allocation_info)

    remaining = quota - assigned_total

    print(f"\n[FLOOR ASSIGNED TOTAL] {assigned_total}")
    print(f"[REMAINING TO DISTRIBUTE] {remaining}")

    allocation_info.sort(
        key=lambda x: x["remainder"],
        reverse=True
    )

    for i in range(remaining):
        allocation_info[i]["base"] += 1

    #
    # Respect group sizes
    #
    allocation = {}

    for info in allocation_info:

        db = info["db"]
        group_size = len(info["group"])

        allocation[db] = min(info["base"], group_size)

    current_total = sum(allocation.values())

    print(f"[TOTAL AFTER ROUNDING] {current_total}")

    #
    # If any group saturated, redistribute leftovers
    #
    while current_total < quota:

        added = False

        for info in allocation_info:

            db = info["db"]
            group_size = len(info["group"])

            if allocation[db] < group_size:

                allocation[db] += 1
                current_total += 1
                added = True

                if current_total == quota:
                    break

        if not added:
            break

    print(f"[FINAL ASSIGNED TOTAL] {current_total}")

    sampled_parts = []

    print("\n[FINAL DATABASE ALLOCATION]")

    for info in allocation_info:

        db = info["db"]
        group = info["group"]

        n_samples = allocation[db]

        print(
            f"  DB {db} | "
            f"available={len(group)} | "
            f"assigned={n_samples}"
        )

        if n_samples > 0:
            sampled_parts.append(
                group.sample(
                    n=n_samples,
                    random_state=seed
                )
            )

    sampled_df = pd.concat(sampled_parts, ignore_index=False)

    print(f"[TOTAL ASSIGNED IN BIN] {current_total}")
    print(f"[ACTUAL SAMPLED FROM BIN] {len(sampled_df)}")

    #
    # Safety check
    #
    if len(sampled_df) != quota:

        print(
            f"[WARNING] Expected {quota} rows "
            f"but obtained {len(sampled_df)}"
        )

        if len(sampled_df) < quota:

            deficit = quota - len(sampled_df)

            print(f"[FILLING DEFICIT] {deficit}")

            remaining_df = df.drop(sampled_df.index)

            extra = remaining_df.sample(
                n=deficit,
                random_state=seed
            )

            sampled_df = pd.concat(
                [sampled_df, extra],
                ignore_index=False
            )

            print(
                f"[SIZE AFTER DEFICIT FILL] "
                f"{len(sampled_df)}"
            )

    print(f"[FINAL BIN SAMPLE SIZE] {len(sampled_df)}")

    return pa.Table.from_pandas(
        sampled_df,
        preserve_index=False
    )


def sample_from_bin(parquet_file, quota, database_column, rng, seed):

    print(f"\n========== PROCESSING BIN ==========")
    print(f"[BIN FILE] {parquet_file}")

    table = pq.read_table(parquet_file)

    sampled = stratified_sample_by_database(table, quota, database_column, rng, parquet_file, seed)

    return sampled


def scan_bin_sizes(parquet_files):

    print("\n========== STEP 2: SCANNING BIN SIZES ==========")

    bin_sizes = {}
    total_population = 0

    for f in parquet_files:

        pf = pq.ParquetFile(f)
        size = pf.metadata.num_rows

        bin_sizes[f] = size
        total_population += size

        print(f"[BIN SIZE] {f}  {size}")

    print(f"\n[TOTAL POPULATION ACROSS BINS] {total_population}")

    return bin_sizes, total_population


def compute_bin_quota(bin_sizes, target_total):

    print("\n========== STEP 3: COMPUTING ADAPTIVE BIN QUOTA ==========")

    remaining_bins = dict(bin_sizes)
    remaining_target = target_total

    quota = {}
    iteration = 0

    while remaining_bins:

        iteration += 1
        print(f"\n[QUOTA ITERATION {iteration}]")

        avg_quota = remaining_target / len(remaining_bins)
        print(f"[CURRENT AVG QUOTA] {avg_quota:.2f}")

        finished = []

        for bin_file, size in remaining_bins.items():

            if size <= avg_quota:

                quota[bin_file] = size
                remaining_target -= size
                finished.append(bin_file)

                print(f"[TAKE ALL] {bin_file} {size}")

        if not finished:
            
            n_bins = len(remaining_bins)
            base_quota = remaining_target // n_bins
            remainder = remaining_target % n_bins

            for i, bin_file in enumerate(remaining_bins):
                quota[bin_file] = base_quota

                if i < remainder:
                    quota[bin_file] += 1

                print(f"[ASSIGN PROPORTIONAL] {bin_file} {quota[bin_file]}")
            break

        for f in finished:
            del remaining_bins[f]

    total_assigned = sum(quota.values())
    print(f"\n[TOTAL ASSIGNED ACROSS BINS] {total_assigned}")
    print(f"[TARGET TOTAL] {target_total}")

    return quota


def merge_diverse_sample_set(input_dir, output_file, total_samples, database_column, seed):

    print("\n========== STARTING SAMPLING PIPELINE ==========")

    rng = random.Random(seed)

    parquet_files = get_parquet_files(input_dir)

    bin_sizes, total_population = scan_bin_sizes(parquet_files)
    
    if total_samples > total_population: 
        print("\n[WARNING] Requested sample larger than dataset")
        total_samples = total_population
    
    bin_quota = compute_bin_quota(bin_sizes, total_samples)

    sampled_tables = []
    total_collected = 0

    print("\n========== STEP 2: BIN SAMPLING ==========")
    
    for file in parquet_files:
        quota = bin_quota[file]
        sampled = sample_from_bin(file, quota, database_column, rng, seed) 
        sampled_tables.append(sampled)

        total_collected += len(sampled)
        print(f"[RUNNING TOTAL COLLECTED] {total_collected}")
    
    print("\n========== STEP 3: MERGING TABLES ==========")

    merged_table = pa.concat_tables(sampled_tables, promote=True)
    print(f"[MERGED SIZE BEFORE TRIMMING] {len(merged_table)}")

    if len(merged_table) > total_samples:

        print("\n========== STEP 4: TRIMMING ==========")
        print(f"[TRIMMING FROM {len(merged_table)} TO {total_samples}]")

        indices = rng.sample(range(len(merged_table)), total_samples)
        merged_table = merged_table.take(pa.array(indices))

    print("\n========== FINAL DATABASE DISTRIBUTION ==========")

    df_final = merged_table.to_pandas()

    final_counts = df_final.groupby(database_column).size()

    for db, count in final_counts.items():

        proportion = count / len(df_final)
        print(f"  DB {db}: {count} ({proportion:.4f})")

    print("\n========== WRITING OUTPUT ==========")

    pq.write_table(merged_table, output_file)

    print(f"[FINAL DATASET SIZE] {len(merged_table)}")
    print(f"[OUTPUT FILE] {output_file}")
    print("\n========== SAMPLING COMPLETE ==========")


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Physicochemical + DB balanced sampler")

    parser.add_argument("--input_dir",
                        help="",
                        required=True)
    parser.add_argument("--outname",
                        help="",
                        required=True)
    parser.add_argument("--sample_size",
                        help="",
                        type=int,
                        required=True)
    parser.add_argument("--db_col",
                        help="",
                        default="db_id")
    parser.add_argument("--seed",
                        help="",
                        type=int,
                        default=1234)

    args = parser.parse_args()

    input_dir = args.input_dir
    outname = args.outname
    sample_size = args.sample_size
    db_col = args.db_col
    seed = args.seed

    merge_diverse_sample_set(input_dir, outname, sample_size, db_col, seed)
