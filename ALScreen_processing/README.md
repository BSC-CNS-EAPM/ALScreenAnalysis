# ALScreen processing

Scripts for processing Glide docking outputs and analysing ALScreen predictions. They select the best pose per ligand, filter structures, compare score distributions, evaluate prediction errors, and generate the tables and figures for Round 1&4 and the Supplementary Information.

Examples:

```bash
python process_glide_csv_output.py \
  --csv ../data/ALScreen/IL1b/unbias/pred2.csv \
  -o results/IL1b_unbias_pred2
python rank_ref_cmpds.py \
  --db_csv ../data/ALScreen/IL1b/unbias/pred2.csv \
  --ref_csv ../data/ALScreen/IL1b/unbias/reference.csv \
  -o results/IL1b_unbias_pred2
python glide_score_threshold.py \
  --db_csvs \
    ../data/ALScreen/IL1b/unbias/initial.csv \
    ../data/ALScreen/IL1b/unbias/pred1.csv \
    ../data/ALScreen/IL1b/unbias/pred2.csv \
  --gscore_thresh -7 -7.5 -8 -8.5 -9 -9.5 -10 -10.5 \
  -o results/IL1b_unbias
python plot_glide_score_AL_hist.py \
  --db_csvs \
    ../data/ALScreen/IL1b/unbias/initial.csv \
    ../data/ALScreen/IL1b/unbias/pred1.csv \
    ../data/ALScreen/IL1b/unbias/pred2.csv \
  --ref_csvs ../data/ALScreen/IL1b/unbias/reference.csv \
  --db_colors grey blue red \
  --ref_colors red \
  --db_labels 'Initial 50k' 'Prediction 1' 'Prediction 2' \
  --ref_labels 'Reference Compounds' \
  -o results/IL1b_init_pred1_pred2_hist_unbias
python training_score_diff.py \
  --csv ../data/ALScreen/IL1b/unbias/pred2_alscreen.csv \
  -n 7500 \
  -o results/IL1b_unbias_train_3_data
python plot_grid_gscore_histograms.py
python prediction_errors.py
python statistical_analysis.py
```

Default output files go into `ALScreen_processing/results/`.

`filter_maegz.py` is the exception: no `.maegz` pose viewer is supplied, so it still
requires an additional original Glide file and the Schrödinger environment.


## Script description

- `process_glide_csv_output.py` selects the lowest Glide score for each ligand, using either raw Glide exports or simplified docking tables. It writes `<output>_best.csv` with all original columns and `<output>_best_simple.csv` with `ID`, `SMILES`, and `glide_gscore`.
- `rank_ref_cmpds.py` ranks reference compounds within the combined database and reference dataset by ascending `glide_gscore`, with lower scores ranked first. It prints the database's minimum and maximum scores and writes each reference compound's one-based rank over the combined total into `<output_ref_rank.txt>`.
- `glide_score_threshold.py` counts compounds with `glide_gscore` strictly below each specified threshold, separately for each input dataset. It writes `<output>_gscore_threshold.txt`.
- `plot_glide_score_AL_hist.py` overlays Glide score histograms and kernel density curves for the supplied datasets. Dataset means appear as solid vertical lines, and reference scores can be added as dashed lines. It prints each dataset's mean and standard deviation and writes `<output>.pdf`, plus `<output>_zoom.pdf` when axis limits are specified.
- `training_score_diff.py` computes MAE and RMSE between ALScreen predictions (`dockscore`) and Glide docking scores (`glide_gscore`), then selects the `n` smallest and `n` largest absolute errors. It prints both error metrics and writes `<output>_top_<n>k.csv`, `<output>_bottom_<n>k.csv`, and `<output>_diff.csv`; the combined selection contains `ID`, `SMILES`, and `glide_gscore`.
- `plot_grid_gscore_histograms.py` overlays the initial set and both AL predictions in one panel per target and bias. Reference compounds appear as dashed lines, with an inset covering scores from −10 to −6. It writes `score_histograms_grid.png` and `.pdf`.
- `prediction_errors.py` computes MAE, RMSE, Pearson correlation, and Spearman correlation for each prediction. Error fold improvement is `prediction 1 error / prediction 2 error`; values above 1 indicate reduced error. It writes `prediction_errors.csv`, `.png`, and `.pdf`; the figure displays MAE, RMSE, and fold improvement.
- `statistical_analysis.py` compares `initial` vs `pred1`, `pred1` vs `pred2`, and `initial` vs `pred2` using two-sided Mann–Whitney U and Kolmogorov–Smirnov tests, plus Cliff's delta. Benjamini–Hochberg correction is applied across all comparisons separately for each test. It writes `statistical_summary.csv` and `cliffs_delta_heatmap.png`/`.pdf`.

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

