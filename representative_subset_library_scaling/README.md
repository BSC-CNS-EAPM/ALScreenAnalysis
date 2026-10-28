# Representative set using stratified sampling

This directory contains a seven-step workflow for selecting a diverse molecular subset from the **Enamine REAL 5B library**. Initial sampling is stratified by heavy atom count (HAC), followed by scaffold filtering, physicochemical binning and fingerprint similarity filtering.

The workflow combines coverage of molecular size and physicochemical properties with reduced structural redundancy.

## Sampling strategy

```text
Enamine REAL 5B, organized by HAC
    ↓
1. Allocate HAC quotas and extract an initial sample
    ↓
2. Retain one molecule per Bemis–Murcko scaffold
    ↓
3. Partition by MW, LogP and TPSA quartiles: 4 × 4 × 4 bins
    ↓
4. Identify similar molecules using Morgan fingerprints and TAPNN
    ↓
5. Remove redundant molecules within each bin
    ↓
6. Sample across filtered bins and merge the final subset
    ↓
7. Calculate iSIM and visualize the subset using UMAP
```

The scripts were originally written to support multiple source databases. For this application, **all input records and population counts describe Enamine REAL 5B only**. A single consistent `db_id` is retained for compatibility with the code. With one source, source/HAC stratification reduces to HAC stratification, and source balancing within each property bin has no additional effect.

Stage 1 can allocate quotas proportionally to the library's HAC distribution or apply weight smoothing to increase the contribution of less populated HAC classes. Subsequent scaffold filtering, similarity filtering and bin balancing change the sampled distribution. The final subset should therefore be described as a diversity-oriented stratified sample.

## Input data

### Molecular records

Store Parquet files directly inside directories named for their HAC:

```text
enamine_real_5b/
    HAC_20/
        part_001.parquet
        part_002.parquet
    HAC_21/
        part_001.parquet
    ...
```

The complete workflow requires:

| Column | Description |
| --- | --- |
| `SMILES` | Molecular structure in SMILES notation. |
| `db_id` | One consistent integer-compatible identifier for the Enamine REAL 5B dataset. |
| Molecular identifier, if available | Retain the original library identifier to connect sampled molecules with docking results and subsequent retrieval. The scripts do not require a particular identifier column name. |

### HAC population counts

Provide a CSV containing `db_id`, `hac_avg` and `count`. The counts must match the supplied Enamine REAL 5B files. Repeated rows for a source/HAC combination are summed. The same source identifier must be used in the CSV and molecular records.

## Script descriptions

### 1. HAC-stratified extraction — `1_database_sampler.py`

This script reads the population CSV, restricts it to the selected HAC range and calculates sampling quotas. For a single source, the population weight of HAC class `h` is:

```text
p(h) = count(h) / SUM(counts in the selected HAC range)
```

Without smoothing, the expected quota is `sample_size × p(h)`. With `--uniform --u_lambda alpha`, the weights become:

```text
w(h) = p(h)^alpha / SUM(p(h)^alpha)
expected_quota(h) = sample_size × w(h)
```

For positive-population strata, `alpha = 1` preserves proportional weights, `alpha = 0` gives equal weights, and intermediate values smooth the distribution toward greater representation of smaller strata. For the case of Enamine REAL 5B a `alpha = 0.5` was used. The `--uniform` option enables this weight transformation; it does not establish uniform random sampling across every molecular record in a HAC class.

After quota allocation, the script shuffles file order within each HAC directory, reads batches and samples rows without replacement from each encountered batch until its quota is filled. Later batches may never be visited. Accordingly, the implementation should not be described as uniform random sampling over the entire HAC stratum.

| Main option | Purpose |
| --- | --- |
| `--met_file` | HAC population CSV. |
| `--mol_dir` | Root of the HAC-organized Enamine REAL 5B files. |
| `--size` | Requested initial sample size, before downstream filtering. |
| `--outname` | Sampled Parquet output path. |
| `--hac_min`, `--hac_max` | Inclusive HAC interval. |
| `--uniform`, `--u_lambda` | Optional smoothing; supply the exponent when enabling it. |
| `--seed` | Sampling seed; default 1234. |
| `--chunk_size` | Read batch size; default 100,000. |
| `--plot_hac_hist`, `--hac_hist_out` | Optional sampled HAC histogram and output image path. |

Outputs include the sampled data, `population_weight_<u_lambda>.csv` and `sample_quota_<size>.json`. Planned samples above 20 million rows are split into numbered Parquet files in a directory derived from `--outname`.

### 2. Unique-scaffold selection — `2_sample_scaffold.py`

RDKit parses each structure and generates its Bemis–Murcko scaffold SMILES. The script appends a `scaffolds` column and retains the first encountered molecule for each distinct scaffold string. Invalid SMILES are excluded during selection.

- Input: `--sample` from stage 1.
- Controls: `--processors`, `--batch_size`.
- Outputs: `<sample_stem>_refine.parquet` and `<sample_stem>_refine_unique.parquet`.

The retained representative depends on input order. Acyclic molecules generally share an empty Murcko scaffold, so only one such molecule is retained. This stage intentionally removes repeated scaffold representatives; it is not the later approach that allows multiple analogues per scaffold.

### 3. Physicochemical stratification — `3_sample_prop_division.py`

The script calculates molecular weight (`MolWt`), calculated LogP (`MolLogP`) and TPSA using RDKit. It determines the 25th, 50th and 75th percentiles for each property in the scaffold-filtered sample, then assigns molecules to the corresponding combinations of quartiles.

- Input: `--sample` from stage 2.
- Controls: `--processors`, `--batch_size`, `--bins_dir`.
- Output: `<bins_dir>/MWqX_LogPqY_TPSAqZ/data.parquet`, with `MW`, `LogP` and `TPSA` columns added.

There are 64 possible bins, although correlations and tied property values can produce unequal or empty bins. Quartiles describe the processed sample, not the full 5B library. The script makes two passes, computes properties twice and holds the first-pass property values in memory.

### 4. Similarity calculation — `4_get_fps_tapnn.py`

Each populated bin is processed independently. The script calculates binary Morgan fingerprints, exports active bits to a sparse text file and invokes TAPNN to identify molecular pairs meeting its Tanimoto threshold criterion.

- Input directory: `--samples`.
- Fingerprint defaults: radius 4 and 1,024 bits, configurable with `--radius` and `--fpsize`.
- Threshold: `--t_sim`; the historical launcher uses 0.5.
- Other controls: `--processors`, `--batch_size`, `--root_tapnn`.
- Outputs per bin: `data.txt` and `data_sim.txt`.

The CLI also requires `--outname`, although actual result filenames are derived from the input file. TAPNN must be available at `<root_tapnn>/build/apss`.

**Only molecules in the same bin are compared.** This step does not impose a global Tanimoto threshold on the merged sample.

### 5. Similarity filtering — `5_sim_filtering.py`

For each reported non-self pair `(i, j)`, the script marks the larger row index for removal. It streams the molecular file and writes the retained records to `data_filtered.parquet`.

- Input: `--input_dir` containing bins and TAPNN results.
- Default batch size: 250,000 rows.
- Output: `<input_stem>_filtered.parquet` for bins with marked removals.

### 6. Balanced selection and merging — `6_merge_filt_files.py`

The script scans `_filtered.parquet` files and allocates the requested final size approximately equally across available bins. Bins with insufficient molecules contribute all their records; the remaining allocation is redistributed across the other bins.

Within each bin, molecules are sampled without replacement using the specified seed. Because this application contains only Enamine REAL 5B, the original within-bin source allocation reduces to sampling from one source.

- Required: `--input_dir`, `--outname`, `--sample_size`.
- Optional: `--seed` and `--db_col` (default `db_id`).
- Output: the final sampled Parquet file.

If too few molecules remain, the similarity filtering is reduced in order to obtain the required sample size.

### 7. Final-set analysis — `7_final_analysis_isim_umap.py`

The script generates Morgan fingerprints, stores packed fingerprints on disk and accumulates fingerprint bit counts. It uses `calculate_isim(..., n_ary="JT")` to report global and per-source iSIM, with diversity reported as `1 − iSIM`. For a valid single-source dataset with at least two molecules, the global and single-source iSIM values should agree apart from presentation precision.

It also fits a two-dimensional UMAP on a reservoir reference sample using Jaccard distance, transforms the remaining fingerprints in batches and renders the embedding with Datashader. Source coloring has only one category in this application.

- Required: `--infile`.
- Main controls: `--workdir`, `--outname`, `--processors`, `--batch_size`, `--transform_batch_size`.
- Defaults: radius 4, 1,024 bits, UMAP reference size 500,000, 15 neighbors, minimum distance 0.1 and UMAP random state 42.
- Outputs: similarity statistics in the log, fingerprint/intermediate files, embedding parts and a UMAP image.

Use a separate work directory for every analysis run. iSIM and UMAP describe the sampled set; neither verifies a maximum pairwise similarity or demonstrates coverage of every region of the 5B library.

## Environment and configuration

The scripts import NumPy, pandas, PyArrow, RDKit and Matplotlib. Analysis additionally uses `umap-learn`, Dask DataFrame, Datashader, SciPy, seaborn and scikit-learn. A pinned dependency environment is not supplied.

The recovered directory also requires external components that are not included:

- An importable `utils.comp.calculate_isim` implementation. Stages 2 and 3 import it even though their main workflows do not use it; stage 7 uses it for iSIM calculations.
- A compatible compiled TAPNN executable at `<root_tapnn>/build/apss`.

