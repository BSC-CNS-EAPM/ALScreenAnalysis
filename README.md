# ALScreenAnalysis

This repository holds the analysis code behind the paper:

> **Stress-testing Billion-Compound Libraries Beyond Druggable Space: Active Learning for Interleukin PPI Disruptors**
> *DOI: to be added. The paper will be published soon and the link will be added here.*

The study introduces **ALScreen**, an active-learning (AL) framework based on the ScoreFormer deep surrogate docking. The workflow enables scalable virtual screening of ultra-large chemical libraries against three interleukin protein-protein interaction targets: **IL-1β**
(PDB 8C3U), **IL-2** (1PY2) and **IL-36γ** (6P9E). It was used to screen 5 billion compounds from the Enamine REAL library, with scalability evaluated up to 65 billion compounds. Each active-learning cycle trains on real Glide docking scores, predicts over the whole library, docks the top-ranked compounds, and feeds them back into the training pool.
Compounds were subsequently selected for experimental validation by surface plasmon resonance (SPR) across four rounds, each using a different selection criterion.

The scripts collected here are the ones that produce the figures and
tables reported for those rounds.

![ALScreen workflow and selection criteria](figures/figure1.png)

Every target was screened under 2/3 docking setups, which we call bias conditions: `unbias`
(no constraints), `bias1` and `bias2` (hydrophobic and/or H-bond constraints). They are described in [data/README.md](data/README.md).

## Contents in each folder

| Folder | Paper section | What it does |
|---|---|---|
| [ALScreen/](ALScreen/) | **Round 1**, medicinal chemistry-driven selection, and Supplementary | Glide pose handling, score distributions across the active-learning rounds, DSD prediction error, and the statistical comparison between rounds |
| [interaction_fingerprints/](interaction_fingerprints/) | **Round 2**, interaction fingerprint-driven selection | ProLIF protein-ligand interaction fingerprints of the docked poses, their clustering, and the per-cluster interaction frequencies |
| [analogue_analysis/](analogue_analysis/) | **Round 3**, analogue expansion | Analysis of the analogues selected around the confirmed hits |
| [representative_subset_library_scaling/](representative_subset_library_scaling/) | **Round 4**, representative subset and library scaling | 50K vs 500K seeding of the active-learning cycle, and the scaling from Enamine REAL 5B to 65B |
| [large_libraries_properties/](large_libraries_properties/) | Supplementary | Molecular weight, linearity ratio and unique Bemis-Murcko scaffolds of the full Enamine REAL 5B and 65B libraries |
| [reference_compounds/](reference_compounds/) | Supplementary | Properties (KD, ligand efficiency, MW, linearity ratio) of the known interleukin binders and of the compounds we tested experimentally |
| [chemical_diversity/](chemical_diversity/) | Supplementary | Chemotype and scaffold diversity along the rounds: diversity metrics, Sankey diagrams, and a UMAP of the chemical space explored |
| [data/](data/) | | Input data, downloaded from Zenodo (see below) |
| [figures/](figures/) | | Workflow figure |

Each script writes its output to a `results/` folder next to it. Those folders are not tracked
by git, except for the figures and tables that made it into the paper.

## Installation

```bash
conda env create -f environment.yml
conda activate alscreen-analysis
```

That environment covers everything except [ALScreen/filter_maegz.py](ALScreen/filter_maegz.py),
which reads Glide pose viewer files through the Schrödinger Python API and therefore needs a
**Schrödinger license and installation**. Run it with the Schrödinger interpreter instead:

```bash
$SCHRODINGER/run ALScreen/filter_maegz.py --maegz glide_pv.maegz -o filtered --max_score -8 --neutral
```

## Data

The input data is too large for git, so it lives in Zenodo:

> **Zenodo DOI: 10.5281/zenodo.XXXXXXX** *(placeholder, the final record will be added here)*

Download the archive and extract it into [data/](data/). The folder structure it should end up
with, the columns of every file, and the definition of the bias conditions are all described in
[data/README.md](data/README.md).

The library-wide analyses in `large_libraries_properties/` also need the Enamine REAL libraries
themselves (5B and 65B) as Parquet files. Those are distributed by
[Enamine](https://enamine.net/compound-collections/real-compounds) and are not redistributed here.

## Examples of how to run the analyses

A few examples to give you an idea of how the scripts are used. They are all run from inside
their own folder, and the ones that compute something have to go before the ones that plot it.
Each script explains what it does and what it needs in the docstring at the top of the file,
and the ones that take arguments also respond to `python <script>.py --help`, so those are the
best places to look when you want something other than what is shown here.

### Round 1 and the supplementary docking analysis

Say you are starting from a raw Glide output and want the score distributions and the
statistics that go with them:

```bash
cd ALScreen
python get_best_glide_pose.py --csv glide_output.csv   # only if you start from raw Glide output
python plot_score_histograms.py                        # gscore distributions per target and bias
python prediction_errors.py                            # MAE, RMSE and correlation of the DSD predictions
python statistical_analysis.py                         # Mann-Whitney, KS and Cliff's delta between rounds
```

The three analysis scripts read the data straight from `data/DSD/`, so if you downloaded the
Zenodo archive you can skip the first line and run them as they are.

### Round 2, interaction fingerprints

Computing the fingerprints is the slow part, so you may want to try a single target first:

```bash
cd interaction_fingerprints
python compute_fingerprints.py IL1b     # one target; leave it out to do all three
python analyse_fingerprints.py          # DBSCAN clustering and interaction-frequency heatmaps
```

### Round 3, analogue expansion

See [analogue_analysis/README.md](analogue_analysis/README.md).

### Round 4, representative subset and library scaling

See [representative_subset_library_scaling/README.md](representative_subset_library_scaling/README.md).

### Supplementary, reference compounds

Properties first, plot second. This one is worth running early, because the library-wide plots
further down reuse its output:

```bash
cd reference_compounds
python compute_properties.py            # KD, MW, heavy atoms, ΔG, ligand efficiency, linearity ratio
python plot_strip_plot.py               # hits of each round against the reference compounds
```

### Supplementary, chemical diversity

The chemotypes are computed once and then reused by the three plots, which you can run in any
order or on their own:

```bash
cd chemical_diversity
python compute_chemotypes.py            # chemotypes and scaffold chemotypes per round
python plot_chemical_diversity.py
python plot_sankey_grid.py
python plot_umap_grid.py
```

### Supplementary, large libraries properties

These go over the full Enamine REAL libraries, so they are meant for an HPC cluster rather than
a laptop. The `run_*.sh` scripts are SLURM jobs, and you will need to edit the account, QoS and
module lines for your own cluster. Finished partitions are skipped, so if a job is killed you
can resubmit the same command and it picks up where it left off.
[large_libraries_properties/README.md](large_libraries_properties/README.md) goes through all
of this in detail, including how the library has to be prepared.

As an example, here is the molecular weight distribution of both libraries, with the 65B one
split over 48 nodes and the 5B one over 7:

```bash
cd large_libraries_properties
mkdir -p logs
sbatch --array=0-47 run_calculate_mw_database.sh /path/to/Enamine_REAL_65B_parquet results/mw/enamine_65B
sbatch --array=0-6  run_calculate_mw_database.sh /path/to/Enamine_REAL_5B_parquet  results/mw/enamine_5B
python plot_mw_database.py results/mw/enamine_5B results/mw/enamine_65B
```

The linearity ratios and the unique scaffolds follow the same pattern:

```bash
sbatch run_calculate_linearity_ratios_database.sh /path/to/Enamine_REAL_65B_parquet results/linearity/enamine_65B
sbatch run_calculate_linearity_ratios_database.sh /path/to/Enamine_REAL_5B_parquet  results/linearity/enamine_5B
python plot_linearity_ratios_database.py results/linearity/enamine_5B results/linearity/enamine_65B

sbatch run_calculate_unique_scaffolds_database.sh /path/to/Enamine_REAL_65B_parquet results/scaffolds/enamine_65B
sbatch run_calculate_unique_scaffolds_database.sh /path/to/Enamine_REAL_5B_parquet  results/scaffolds/enamine_5B
```

Both plotting scripts draw the reference compounds on top of the library distributions, so
`reference_compounds/compute_properties.py` has to have been run first.

## Citation

If you use this code, please cite the paper (reference and DOI to be added once it is published).

## Contact

Questions about the code or the data are welcome:

**Júlia Vilalta Mor**, [julia.vilalta@bsc.es](mailto:julia.vilalta@bsc.es)

## Acknowledgments

This project, with reference number CPP2022-009737, has been funded by the State Research Agency
(AEI) through the State Program to Promote Scientific-Technical Research and its Transfer, part
of the State Plan for Scientific, Technical, and Innovation Research 2021-2023, within the
framework of the Recovery, Transformation, and Resilience Plan (Next Generation EU funds).
