# ALScreenAnalysis

This repository holds the analysis code behind the paper:

> **Stress-testing Billion-Compound Libraries Beyond Druggable Space: Active Learning for Interleukin PPI Disruptors**
>
> *DOI: to be added once the paper is published.*

The study introduces **ALScreen**, an active-learning framework based on the ScoreFormer deep surrogate docking model. It screens ultra-large chemical libraries against three interleukin protein-protein interaction targets: **IL-1β**, **IL-2** and **IL-36γ**. ALScreen was initially applied to exhaustively evaluate 5B compounds from the Enamine REAL database, and furthered scaled to evaluate up to 65 billion compounds from Enamine REAL Space.

The scripts here prepare and analyse the data for the figures and tables reported in the study. The links related to the data required for such analysis are found within the [data/README.md](data/README.md).

![ALScreen workflow and selection criteria](figures/figure1.png)


## Contents in each subdirectory

| Subdirectory | Paper section | What it does |
|---|---|---|
| [ALScreen_processing/](ALScreen_processing/) | **Rounds 1 and 4**, and Supplementary | Glide pose processing, score distributions, prediction errors and statistical comparisons |
| [interaction_fingerprints/](interaction_fingerprints/) | **Round 2**, interaction fingerprint-driven selection | Protein-ligand interaction fingerprints, clustering and interaction frequencies |
| [analogue_analysis/](analogue_analysis/) | **Round 3**, analogue expansion | Analogue preparation, Glide/PELE comparisons and final selection |
| [representative_subset_library_scaling/](representative_subset_library_scaling/) | **Round 4**, representative subset and library scaling | Construction and analysis of a diverse Enamine REAL 5B subset for active-learning seeding |
| [large_libraries_properties/](large_libraries_properties/) | Supplementary | Molecular weight, linearity ratio and unique Bemis-Murcko scaffolds of the 5B and 65B libraries |
| [reference_compounds/](reference_compounds/) | Supplementary | Molecular and affinity-derived properties of known binders and experimentally tested compounds |
| [chemical_diversity/](chemical_diversity/) | Supplementary | Chemotype and scaffold diversity, Sankey diagrams and chemical-space UMAPs |
| [data/](data/) | | Input data and format descriptions |
| [figures/](figures/) | | Workflow figure |

Most analysis figures and tables are written to `results/` folders. Some preparation scripts write beside their inputs or accept an output path. See each folder’s README for details.

## Installation

Create and activate the base environment:

```bash
conda env create -f environment.yml
conda activate alscreen-analysis
```

Additional components are required for particular scripts:

- **Schrödinger installation and license** for `ALScreen_processing/filter_maegz.py`.

Check the corresponding subdirectory README for specific dependencies and configuration.

## Data

The analysis data is too large for git and is deposited in Zenodo:

> **Zenodo DOI: to be added.**

Download all ZIP files and exrtact them into [data/](data/). 

Library-wide analyses require the Enamine REAL 5B and 65B libraries as partitioned Parquet files. These are distributed by [Enamine](https://enamine.net/compound-collections/real-compounds) and are not redistributed here. 
Similarly, data generated from the analogue search analysis are not included due to their large volume. These data are available upon reasonable request to the corresponding author of the paper.


## How to run the analyses

Execute the following scripts from the repository root. Follow the execution order from the blocks so that intermediate data required for the plots or other analysis are generated. Subdirectory READMEs describe inputs, outputs and parameters; scripts accepting arguments also provide `--help`.

### Round 1 to 4 and supplementary analysis

Process raw Glide output if necessary, then generate score distributions, prediction errors and statistical comparisons:

```bash
cd ALScreen_processing

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

See [ALScreen_processing/README.md](ALScreen_processing/README.md) for additional filtering, ranking and plotting tools.

### Round 2, interaction fingerprints

Compute fingerprints once from docked poses and a prepared receptor, then reuse the saved file for clustering and analysis:

```bash
cd interaction_fingerprints

python compute_fingerprints.py --target IL1b \
    --pdb ../data/interaction_fingerprints/IL1b/receptor.pdb \
    --sdf ../data/interaction_fingerprints/IL1b/poses.sdf

python analyse_fingerprints.py --target IL1b --ifp results/IL1b_ifp.pkl
```

See [interaction_fingerprints/README.md](interaction_fingerprints/README.md) for input details.

### Round 3, analogue expansion

The workflow prepares and intersects analogue libraries, examines molecular properties
and chemical space, compares Glide docking scores with PELE metrics, and selects
candidates using docking scores, visual inspection and PELE improvement relative to
the parent hit.

Because of their large volume, the full analogue libraries, PELE outputs and
visual-inspection pose datasets are not included in the supplied data. The complete
Round 3 workflow therefore cannot be executed using the supplied data alone. These data are available upon reasonable request to the corresponding author of the paper. 

See [analogue_analysis/README.md](analogue_analysis/README.md) for the manifest format,
library preparation, helper dependencies and selection criteria.

### Round 4, representative subset and library scaling

The representative-subset workflow starts from Enamine REAL 5B and applies HAC-stratified extraction, unique-scaffold selection, physicochemical binning, within-bin fingerprint similarity filtering and balanced merging. The final subset is assessed using iSIM and UMAP.

You must download the Enamine REAL 5B database and prepare it in the input format
described in the workflow README.

See [representative_subset_library_scaling/README.md](representative_subset_library_scaling/README.md) for the complete workflow and external requirements.

### Supplementary, reference compounds

Calculate properties for each reference or experimental table before generating the comparison plot:

```bash
cd reference_compounds

python compute_properties.py --input ../data/reference_compounds/IL1b_reference.csv

python compute_properties.py --input ../data/experimental_hits/IL1b_round1.csv \
    --affinity_column KD_uM --affinity_units uM

python compute_properties.py --input ../data/experimental_hits/IL1b_round2.csv \
    --affinity_column KD_uM --affinity_units uM

python compute_properties.py --input ../data/experimental_hits/IL1b_round3.csv \
    --affinity_column KD_uM --affinity_units uM

python compute_properties.py --input ../data/reference_compounds/IL2_reference.csv

python compute_properties.py --input ../data/experimental_hits/IL2_round1.csv \
    --affinity_column KD_uM --affinity_units uM

python compute_properties.py --input ../data/experimental_hits/IL2_round2.csv \
    --affinity_column KD_uM --affinity_units uM

python compute_properties.py --input ../data/experimental_hits/IL2_round3.csv \
    --affinity_column KD_uM --affinity_units uM

python compute_properties.py --input ../data/reference_compounds/IL36g_reference.csv

python compute_properties.py --input ../data/experimental_hits/IL36g_round1.csv \
    --affinity_column KD_uM --affinity_units uM

python compute_properties.py --input ../data/experimental_hits/IL36g_round2.csv \
    --affinity_column KD_uM --affinity_units uM

python compute_properties.py --input ../data/experimental_hits/IL36g_round3.csv \
    --affinity_column KD_uM --affinity_units uM

python plot_strip_plot.py
```

Repeat property calculation for all datasets needed by the plot. See [reference_compounds/README.md](reference_compounds/README.md).

### Supplementary, chemical diversity

Compute chemotypes first, then generate the plots independently:

```bash
cd chemical_diversity

python compute_chemotypes.py

python plot_chemical_diversity.py

python plot_sankey_grid.py

python plot_umap_grid.py
```

### Supplementary, large libraries properties

These analyses process full libraries and are intended for an HPC cluster. Prepare partitioned Parquet inputs and adapt the SLURM account, QoS, environment and resource settings to your cluster.

For example, calculate molecular weights and plot their distributions:

```bash
cd large_libraries_properties

mkdir -p logs

sbatch --array=0-47 run_calculate_mw_database.sh \
    /path/to/Enamine_REAL_65B_parquet results/mw/enamine_65B

sbatch --array=0-6 run_calculate_mw_database.sh \
    /path/to/Enamine_REAL_5B_parquet results/mw/enamine_5B

python plot_mw_database.py results/mw/enamine_5B results/mw/enamine_65B
```

Linearity-ratio and unique-scaffold calculations follow the same submission pattern. Finished partitions are skipped when jobs are resubmitted. Prepare reference-compound properties before plots that overlay them.

See [large_libraries_properties/README.md](large_libraries_properties/README.md) for input preparation and cluster configuration.

## Citation

If you use this code, please cite the paper. The reference and DOI will be added once it is published.

## Contact

Questions about the code or data are welcome:

- **Isaac Soul Garcia**, [isaac.soulgarcia@bsc.es](mailto:isaac.soulgarcia@bsc.es)
- **Júlia Vilalta Mor**, [julia.vilalta@bsc.es](mailto:julia.vilalta@bsc.es)
- **Dr. Isaac Filella-Merce**, [ifilella@nostrumbiodiscovery.com](mailto:ifilella@nostrumbiodiscovery.com)

## Acknowledgments

This project, with reference number CPP2022-009737, has been funded by the State Research Agency (AEI) through the State Program to Promote Scientific-Technical Research and its Transfer, part of the State Plan for Scientific, Technical, and Innovation Research 20212023, within the framework of the Recovery, Transformation, and Resilience Plan (Next Generation EU funds).
