# Analogue expansion analysis

Scripts for preparing analogue libraries, examining molecular properties and chemical space, comparing Glide and PELE scores, and selecting analogues relative to experimental parent hits.

The full analogue libraries, PELE results and visually accepted poses are omitted from the supplied data because of their size. The examples below require researcher-supplied inputs; run them from `analogue_analysis/`. Lower Glide and PELE scores are better, and selection requires scores strictly below the parent or specified cutoff.

## Installation

Use the repository's analysis environment. To run the three helpers (`get_sdf_from_csv.py`, `get_moldb_from_sdf.py` and `intersect_moldbs.py`), you must also clone [MolecularAnalysis](https://github.com/IFilella/MolecularAnalysis.git) and install its dependencies in the environment used by the helpers.

From `analogue_analysis/`, clone it beside ALScreenAnalysis and make it importable:

```bash
git clone https://github.com/IFilella/MolecularAnalysis.git ../../MolecularAnalysis
export PYTHONPATH="$(cd ../.. && pwd)${PYTHONPATH:+:$PYTHONPATH}"
```

## Prepare analogue libraries

### `sdf_creator.py`

Converts the CSV files in each input folder to SDF using the CSV-to-SDF helper.

```bash
python sdf_creator.py --directories data/analogue_hunter/IL1b \
  --script ./get_sdf_from_csv.py --id-column name --smiles-column smiles \
  --properties space smiles --delimiter ','
```

Flags: `--directories`: input folders; `--script`: helper path; `--id-column`: compound ID column; `--smiles-column`: SMILES column; `--properties`: columns copied to SDF properties; `--delimiter`: CSV separator.

### `moldb_creator.py`

Converts the SDF files in each input folder to serialized MolDB libraries.

```bash
python moldb_creator.py --directories data/analogue_hunter/IL1b \
  --script ./get_moldb_from_sdf.py
```

Flags: `--directories`: input folders; `--script`: SDF-to-MolDB helper path.

### `intersect_analogues.py`

Intersects matching MolDB files from two source folder trees and exports the shared compounds.

```bash
python intersect_analogues.py --directories data/analogue_hunter/IL1b \
  --other-base data/motif_matcher --script ./intersect_moldbs.py \
  --similarity 1 --output-dir results/intersections
```

Flags: `--directories`: first-source folders; `--other-base`: second-source parent folder with matching subfolder names; `--script`: intersection helper; `--similarity`: similarity cutoff; `--output-dir`: output folder; `--no-sdf`: disable SDF export.

All three batch wrappers also accept `--python` to select the helper's interpreter and `--dry-run` to print commands without executing them.

## Run individual helpers

### `get_sdf_from_csv.py`

Converts one compound CSV to SDF, retaining the requested properties.

```bash
python get_sdf_from_csv.py --csv data/analogues.csv --id name --smi smiles \
  --params space smiles -d ',' -o results/analogues
```

Flags: `--csv`: input CSV; `--id`: ID column; `--smi`: SMILES column; `--params`: exported property columns; `-d`: separator; `-o`: output prefix; `--sim_filt`: supply `True` to enable similarity filtering, or omit to disable it.

### `get_moldb_from_sdf.py`

Converts one or more SDF files to MolDB libraries.

```bash
python get_moldb_from_sdf.py --sdfs results/analogues.sdf -o results/analogues
```

Flags: `--sdfs`: input SDF files; `-o` / `--outname`: corresponding output prefixes, one per input.

### `intersect_moldbs.py`

Finds the intersection of two MolDB libraries at the chosen similarity cutoff.

```bash
python intersect_moldbs.py --moldb1 results/source_A.moldb \
  --moldb2 results/source_B.moldb --sim 1 --sdf -o results/analogues_intersect
```

Flags: `--moldb1`, `--moldb2`: input libraries; `--sim`: similarity cutoff; `--sdf`: also export SDF; `-o` / `--outname`: output prefix.

## Molecular properties and chemical space

### `get_metrics.py`

Summarises supplied molecular-property columns or calculates molecular weight, LogP and TPSA from SMILES.

```bash
python get_metrics.py --directory data/analogues --compute-from-smiles \
  --smiles-column smiles -o results/property_summary.csv
```

Flags: `--directory`: input folder; `--pattern`: filename pattern (default `*.csv`); `--compute-from-smiles`: calculate properties; `--smiles-column`: SMILES column; `--columns`: existing property columns when not calculating; `-o` / `--output`: summary CSV.

### `tsne_plot.py`

Plots the chemical space of a CSV or SDF using Morgan fingerprints and t-SNE, optionally highlighting a parent hit.

```bash
python tsne_plot.py --input data/union/HIT.csv \
  --smiles-column smiles --group-column space \
  --experimental-csv data/experimental_hits.csv --hit HIT \
  --experimental-delimiter ';' --radius 4 --n-bits 2048 \
  --perplexity 50 --iterations 5000 --seed 42 -o results/HIT_tsne
```

Flags: `--input`: CSV or SDF; `--smiles-column`: CSV SMILES column; `--group-column`: group column or SDF property; `--delimiter`: input CSV separator; `--experimental-csv`: optional parent-hit table; `--hit`: parent ID (required with that table); `--experimental-id-column`, `--experimental-smiles-column`: parent-table columns; `--experimental-delimiter`: parent-table separator; `--radius`, `--n-bits`: Morgan fingerprint settings; `--perplexity`: t-SNE neighbourhood setting, smaller than the sample count; `--iterations`: optimisation steps; `--seed`: random seed; `--metric`: `euclidean` or `jaccard`; `-o` / `--output-prefix`: PDF and embedding-CSV prefix.

## Docking and PELE inputs

`pele_plots.py` and `final_selection.py` require a CSV manifest with one row per target, parent hit and constraint. File paths are relative to the manifest location.

| Field | Description |
| --- | --- |
| `target`, `hit`, `constraint` | Target label, exact parent ID and constraint label |
| `analogue_csv`, `hit_csv`, `pele_csv` | Analogue docking, parent docking and PELE tables |
| `visual_dir` | Accepted PDB poses; required for final selection |
| `pose_dir` | Source PDB poses; required with `--copy-poses` |
| `gscore_threshold`, `pele_threshold` | Optional additional score cutoffs; blank means none |
| `pele_suffix`, `visual_suffix`, `pose_suffix` | Optional naming suffixes; default to `constraint` |

Docking tables use `ID`, `SMILES` and `glide_gscore`; plots require only ID and score. PELE tables use `molecule` and the requested metric columns. PELE identifiers follow `<ID>_<constraint>` and pose filenames follow `<ID>_<constraint>.pdb`, unless suffix overrides are supplied. Provide best-pose tables with unique IDs per constraint. `original_thresholds.csv` lists the original Glide cutoffs for optional use in the manifest.

### `pele_analysis.py`

Counts candidates with PELE metrics below the chosen parent-hit condition and saves comparison and summary CSVs.

```bash
python pele_analysis.py --pele_df data/pelemetrics.csv --hit HIT \
  --const HP_HB --constraints unbias HP HB HP_HB \
  --metrics all_BFE lowBFE_minBE --output-prefix results/HIT_pele
```

Flags: `--pele_df` / `--pele-csv`: PELE table; `--hit`: parent ID; `--const` / `--reference-constraint`: parent condition used for comparison; `--constraints`: all constraint suffixes present; `--metrics`: PELE columns to compare; `--molecule-column`: identifier column; `--output-prefix`: output prefix.

### `pele_plots.py`

Plots Glide scores against PELE metrics by constraint and reports correlations and improvements relative to the parent hit.

```bash
python pele_plots.py --manifest manifest.csv --metrics all_BFE lowBFE_minBE \
  --reference-constraint HP_HB --output-dir results/pele_plots
```

Flags: `--manifest`: input manifest; `--metrics`: PELE columns; `--reference-constraint`: optional additional comparison against one parent condition; `--output-dir`: output folder; `--id-column`, `--score-column`: docking-table columns; `--pele-id-column`: PELE identifier column.

### `final_selection.py`

Selects visually accepted analogues with Glide and PELE scores below the parent hit and any additional manifest cutoffs, reporting their similarity to the parent.

```bash
python final_selection.py --manifest manifest.csv --pele-metric lowBFE_minBE \
  --radius 4 --n-bits 2048 --output-dir results/final_selection --copy-poses
```

Flags: `--manifest`: input manifest; `--pele-metric`: PELE selection column; `--radius`, `--n-bits`: Morgan fingerprint settings; `--output-dir`: output folder; `--copy-poses`: copy selected PDBs; `--id-column`, `--smiles-column`, `--score-column`: docking-table columns; `--pele-id-column`: PELE identifier column.
