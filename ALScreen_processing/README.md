# ALScreen processing

Scripts that process the output of the ALScreen runs, the docking results stored in
`data/ALScreen/<target>/<bias>/`, and turn them into the tables and figures of Round 1 and the
Supplementary Information.

- `plot_grid_gscore_histograms.py` draws the Glide gscore distributions of the initial random
  set and the two AL predictions, one panel per target and bias, with the docked reference
  compounds marked and an inset zooming into the best-scoring tail.
- `prediction_errors.py` compares the ALScreen predicted score with the Glide gscore of the same
  compounds (MAE, RMSE, Pearson and Spearman) and the fold improvement from prediction 1 to
  prediction 2.
- `statistical_analysis.py` tests whether the score distributions really improve between rounds,
  with Mann-Whitney U and Kolmogorov-Smirnov tests, the Cliff's delta effect size and a
  Benjamini-Hochberg correction.

They take no arguments and read the data folder directly, so they can be run in any order. Each
one writes its output to `results/`.

```bash
python plot_grid_gscore_histograms.py
python prediction_errors.py
python statistical_analysis.py
```
