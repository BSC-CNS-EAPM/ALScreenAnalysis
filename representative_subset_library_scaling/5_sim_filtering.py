import os
import pyarrow as pa
import pyarrow.parquet as pq


def dir_traversal(input_dir):
    """
    Finds parquet + similarity files inside bin directories
    """

    file_pairs = []

    for subdir, _, files in os.walk(input_dir):

        parquet_file = None
        sim_file = None

        for f in files:

            if f.endswith(".parquet") and not f.endswith("_filtered.parquet"):
                parquet_file = os.path.join(subdir, f)

            elif f.endswith("_sim.txt"):
                sim_file = os.path.join(subdir, f)

        if parquet_file and sim_file:
            print(f"[PAIR FOUND]")
            print(f"  parquet : {parquet_file}")
            print(f"  sim     : {sim_file}\n")
            file_pairs.append((parquet_file, sim_file))

        elif parquet_file:
            print(f"[SKIPPED] No sim file for {parquet_file}")

    file_pairs.sort()
    return file_pairs


def load_similarity_pairs(sim_file):
    """
    Reads TAPNN similarity output and returns 1-based molecule indices
    to remove.
    """

    remove_set = set()
    total_pairs = 0

    with open(sim_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue

            try:
                i = int(parts[0])
                j = int(parts[1])
            except ValueError:
                continue

            if i==j:
                continue

            remove_set.add(max(i, j))
            total_pairs += 1
    
    print(f"[INFO] Similarity adges read: {total_pairs}")
    print(f"[INFO] Molecules marked for removal: {len(remove_set)}")

    return remove_set


def filter_parquet(input_parquet, remove_indices, batch_size):
    """
    Stream the parquet file in record batches and write the 
    filtered rows to a new parquet file without loading the
    full table into pandas.
    """

    print(f"[READING] {input_parquet}")
    output_file = input_parquet.replace(".parquet", "_filtered.parquet")

    pf = pq.ParquetFile(input_parquet)
    total_rows = pf.metadata.num_rows
    remove_zero = {idx - 1 for idx in remove_indices}
    
    print(f"[SIZE BEFORE] {total_rows}")
    print(f"[WRITING] {output_file}")
    
    writer = None
    rows_written = 0
    global_row = 0

    try:
        for batch_id, batch in enumerate(pf.iter_batches(batch_size=batch_size), start=1):
            table = pa.Table.from_batches([batch])
            batch_rows = table.num_rows

            keep_positions = [pos for pos in range(batch_rows) if (global_row + pos) not in remove_zero]

            if keep_positions:
                filtered = table.take(pa.array(keep_positions,
                                               type=pa.int64()))
                if writer is None:
                    writer = pq.ParquetWriter(output_file, filtered.schema)

                writer.write_table(filtered)
                rows_written += filtered.num_rows

            global_row += batch_rows
            print(f"[BATCH {batch_id}] processed={global_row} kept={rows_written}")

        if writer is None:
            empty_table = pa.Table.from_arrays(
                [pa.array([], type=field.type) for field in pf.schema_arrow],
                names=pf.schema_arrow.names,)

            writer = pq.ParquetWriter(output_file, empty_table.schema)
            writer.write_table(empty_table)
    
    finally:
        if writer is not None:
            writer.close()

    print(f"[SIZE AFTER] {rows_written}\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Remove redundant molecules from TAPNN similarity results")

    parser.add_argument("--input_dir",
                        help="Input directory with bins containing TAPNN results",
                        required=True)
    parser.add_argument("--batch_size",
                        help="Batch size to stream the input file",
                        default=250000)

    args = parser.parse_args()
    input_dir = args.input_dir
    batch_size = args.batch_size

    file_pairs = dir_traversal(input_dir)
    print(f"\nTotal bin pairs detected: {len(file_pairs)}\n")
    
    total_mol_rem = 0

    for parquet_file, sim_file in file_pairs:

        remove_indices = load_similarity_pairs(sim_file)
        total_mol_rem += len(remove_indices)

        if len(remove_indices) == 0:
            print("[SKIP] No redundant molecules detected\n")
            continue

        filter_parquet(parquet_file, remove_indices, batch_size)

    print(f"There is a total of: {total_mol_rem} redundant molecules to remove")
