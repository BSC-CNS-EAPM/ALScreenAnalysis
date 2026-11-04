# Reference compounds

`compute_properties.py` works on **any csv with SMILES data**, one file at a time.

It adds the molecular properties (MW, heavy atoms, ΔG, ligand efficiency and linearity ratio) and writes the result to a new csv, keeping every column of the input.

In the paper we computed them for both the reference compounds and the compounds that reported a KD affinity in our selection rounds.

`plot_strip_plot.py` is the comparison itself, the strip plot shown in the Supplementary
Information: one row per target and one column per property, with the hits of each round placed against the reference compounds.

Example:
```bash
python compute_properties.py --input ../data/reference_compounds/IL1b_reference.csv

python compute_properties.py --input ../data/experimental_hits/IL1b_round1.csv \
    --affinity_column KD_uM --affinity_units uM

python compute_properties.py --input ../data/reference_compounds/IL2_reference.csv

python compute_properties.py --input ../data/experimental_hits/IL2_round1.csv \
    --affinity_column KD_uM --affinity_units uM

python compute_properties.py --input ../data/reference_compounds/IL36g_reference.csv

python compute_properties.py --input ../data/experimental_hits/IL36g_round1.csv \
    --affinity_column KD_uM --affinity_units uM

python plot_strip_plot.py
```
