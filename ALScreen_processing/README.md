# ALScreen processing

Scripts for processing Glide docking outputs and analysing ALScreen predictions. They select the best pose per ligand, filter structures, compare score distributions, evaluate prediction errors, and generate the tables and figures for Round 1&4 and the Supplementary Information.

Lower Glide gscores indicate more favourable docking scores.

## Requirements

The CSV processing and analysis scripts use Python with the following packages:

```bash
python -m pip install numpy pandas matplotlib seaborn scipy statsmodels
```

`filter_maegz.py` requires a Schrödinger installation and its Python environment. Run it with `$SCHRODINGER/run`.

## Prepare docking outputs

### Select the best pose from a Glide CSV

`process_glide_csv_output.py` selects the lowest `r_i_glide_gscore` for each ligand title, after dropping rows without that score. Malformed CSV lines are skipped.

```bash
python process_glide_csv_output.py --csv docking.csv
```

The input must contain `title`, `SMILES`, and `r_i_glide_gscore`. It writes two files alongside the input:

- `docking_best.csv`: selected poses with all original columns.
- `docking_best_simple.csv`: selected poses with columns `id`, `smiles`, and `glide_gscore`.

### Filter a Maestro pose file

`filter_maegz.py` retains the lowest-scoring qualifying pose for each ligand title and writes `<output>.maegz`.

```bash
"$SCHRODINGER/run" filter_maegz.py \
  --maegz docking.maegz -o docking_filtered \
  --metric gscore --metric_threshold -7 \
  --neutral --mw 400 minor
```

Only `--maegz` and `-o` are required. Optional filters are:

- `--metric gscore` (default) or `--metric docking_score`: use `r_i_glide_gscore` or `r_i_docking_score` to select poses.
- `--metric_threshold VALUE`: discard poses with scores above the threshold; equality is retained.
- `--neutral`: discard structures containing any atom with a nonzero formal charge, including structures with net charge zero but charged atoms.
- `--mw VALUE minor`: retain molecular weights at or below `VALUE`.
- `--mw VALUE major`: retain molecular weights at or above `VALUE`.

Use an output stem without `.maegz`; the extension is appended automatically. The script treats structures lacking the selected metric as receptor structures, writes the last such structure first, and expects at least one to be present.

## Analyse individual datasets

### Count compounds below score thresholds

`glide_score_threshold.py` prints the number of rows with `glide_gscore` strictly below each threshold, for each input CSV in the supplied order. It does not write output files.

```bash
python glide_score_threshold.py \
  --db_csvs initial.csv pred1.csv pred2.csv \
  --gscore_thresh -6 -7 -8
```

### Plot score distributions

`plot_glide_score_AL_hist.py` overlays histograms, kernel density curves, and mean-score lines for datasets containing `glide_gscore`. It prints each dataset's mean and standard deviation and can add reference scores as dashed vertical lines.

```bash
python plot_glide_score_AL_hist.py \
  --db_csvs initial.csv pred1.csv pred2.csv \
  --db_labels "Initial 50K" "Prediction 1" "Prediction 2" \
  --db_colors gray royalblue indianred \
  --ref_csvs reference.csv --ref_labels "Reference compounds" \
  --ref_colors firebrick \
  -b 100 --xlim -10 -6 --ylim 0 500 -o glide_histogram
```

It writes `<output>.pdf` and, when `--xlim` or `--ylim` is provided, `<output>_zoom.pdf`. The default output stem is `glide_histogram` and the default bin count is 100. Rows with `glide_gscore >= 50` are excluded.

Use `--ref_gscores -8.5 -7.2` to supply reference values directly. Colours and labels correspond to inputs by position. Prefer one reference input method per run: CSV references and explicit scores share the same colour and label lists.

Implementation details to bear in mind:

- `-f` normalises histogram bars to density; KDE curves always use density scaling, including when histogram bars show counts.
- The main PDF has a fixed y-axis range of 0–5000. Use `--ylim` to set the zoomed PDF's range, particularly with `-f`.

### Rank reference compounds

`rank_ref_cmpds.py` ranks reference compounds within the combined database and reference rows, sorted by ascending `glide_gscore`. Both files must contain `ID` and `glide_gscore`.

```bash
python rank_ref_cmpds.py --ref_csv reference.csv --db_csv database.csv
```

It excludes scores at or above 50, then prints the database's minimum and maximum scores and each reference compound's one-based rank over the combined total. No files are written. IDs should be unique across the combined inputs: duplicate IDs are reported using their first matching rank.

### Extract compounds with the smallest and largest prediction errors

`training_score_diff.py` prints MAE and RMSE, computes the absolute difference between `glide_gscore` and `dockscore`, and extracts the `n` smallest and `n` largest differences.

```bash
python training_score_diff.py --csv merged_predictions.csv -n 1000 -o training
```

The input must contain `ID`, `SMILES`, `glide_gscore`, and `dockscore`. Rows with `glide_gscore >= 50` are excluded. Outputs are:

- `<output>_top_<n>k.csv`: smallest errors, retaining input columns and adding `score_diff`.
- `<output>_bottom_<n>k.csv`: largest errors, retaining input columns and adding `score_diff`.
- `<output>_diff.csv`: concatenated selections with only `ID`, `SMILES`, and `glide_gscore`.

`-n` is a row count, despite the literal `k` in the output filenames. The two selections can overlap when fewer than `2 × n` rows remain.

## Generate the target-by-bias analysis

The following three scripts take no arguments and can be run in any order:

```bash
python plot_grid_gscore_histograms.py
python prediction_errors.py
python statistical_analysis.py
```

They resolve input paths relative to their own location, rather than the current working directory. Place the scripts in a directory whose parent contains `data/ALScreen/`, for example:

```text
project/
├── data/ALScreen/<target>/<bias>/
│   ├── initial.csv
│   ├── pred1.csv
│   ├── pred2.csv
│   ├── reference.csv
│   ├── pred1_alscreen.csv
│   └── pred2_alscreen.csv
└── scripts/
    ├── plot_grid_gscore_histograms.py
    ├── prediction_errors.py
    ├── statistical_analysis.py
    └── results/                         # Created automatically
```

The configured target/bias combinations are:

| Target | Display label | Bias folders |
| --- | --- | --- |
| `IL1b` | IL-1β | `unbias`, `bias1` |
| `IL2` | IL-2 | `unbias`, `bias1`, `bias2` |
| `IL36g` | IL-36γ | `unbias`, `bias1`, `bias2` |

`initial.csv`, `pred1.csv`, `pred2.csv`, and `reference.csv` require `glide_gscore`. The `pred1_alscreen.csv` and `pred2_alscreen.csv` files require paired `dockscore` predictions and observed `glide_gscore` values on the same rows. The scripts do not merge or rename these input files automatically.

- `plot_grid_gscore_histograms.py` overlays the initial set and both AL predictions in one panel per target and bias. Reference compounds appear as dashed lines, with an inset covering scores from −10 to −6. It writes `score_histograms_grid.png` and `.pdf`.
- `prediction_errors.py` computes MAE, RMSE, Pearson correlation, and Spearman correlation for each prediction. Error fold improvement is `prediction 1 error / prediction 2 error`; values above 1 indicate reduced error. It writes `prediction_errors.csv`, `.png`, and `.pdf`; the figure displays MAE, RMSE, and fold improvement.
- `statistical_analysis.py` compares `initial` vs `pred1`, `pred1` vs `pred2`, and `initial` vs `pred2` using two-sided Mann–Whitney U and Kolmogorov–Smirnov tests, plus Cliff's delta. Benjamini–Hochberg correction is applied across all comparisons separately for each test. It writes `statistical_summary.csv` and `cliffs_delta_heatmap.png`/`.pdf`.

All these outputs are saved in `results/` beside the scripts. Re-running a script overwrites its existing outputs.
