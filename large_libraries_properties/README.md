# Large libraries properties

These scripts are meant for the analysis of **ultra-large libraries of compounds**.

In our study, we applied them to **Enamine REAL
5B and 65B** to describe the chemical space the ALScreen was searching through,
and the resulting figures are in the Supplementary Information of the paper.

Three properties are computed over the whole library:

- **Molecular weight**, with `calculate_mw_database.py` and `plot_mw_database.py`. The weights
  are computed once and binned later, so you can choose the bin width when you plot.
- **Linearity ratio**, with `calculate_linearity_ratios_database.py` and
  `plot_linearity_ratios_database.py`. It is the graph diameter divided by the number of heavy
  atoms minus one, so a fully linear chain gives 1 and a compact fused system gives much less.
- **Unique Bemis-Murcko scaffolds**, with `calculate_unique_scaffolds_database.py`, which counts
  how many distinct scaffolds the library actually contains.

Nothing here is specific to Enamine REAL. Any library works as long as it is prepared as
described below, so you can point these at your own collection.

## Preparing the library

Before running anything, make sure you have your libraries **downloaded and prepared in Parquet
format**, split into partition files inside one directory.

Parquet is a much better choice than
csv or SMILES text files at this scale: it is columnar, so the scripts read only the SMILES
column and skip everything else, and it is compressed, which saves a lot of space and a lot of
reading time when you are going over billions of molecules.

The input directory should look like this, with any number of partitions:

```
Enamine_REAL_65B_parquet/
├── part_0000.parquet
├── part_0001.parquet
└── ...
```

Each file needs a column with the SMILES, called `SMILES` by default, and the linearity script
also uses an `ID` column. If your columns are named differently, pass `--smiles_column` and
`--id_column`.

## Running on your own cluster

These jobs go over billions of molecules, so they are written for an HPC cluster.

 The three `run_*.sh` files are SLURM submission scripts, and in each one the lines you
need to modify are marked with a `# --- Cluster specific ---` flag: the account and
QoS your cluster expects, and conda environment.


The cores, wall time and memory constraints at the top are the values we used, so adjust them to
your nodes too.

If a job dies halfway it is not a problem. Finished partitions are skipped, so you can resubmit
exactly the same command and it carries on where it stopped.
