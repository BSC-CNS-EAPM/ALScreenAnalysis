# Data

The data needed to run the git is too large and is deposited in Zenodo. **[DOI: to be added]**

Download and extract the ZIP files to follow the data distribution seen below:


```
data/
├── ALScreen/                            # docking results of the AL screen
│   ├── IL1b/
│   │   ├── unbias/
│   │   │   ├── initial.csv              # random initial 50K set (Glide)
│   │   │   ├── pred1.csv                # compounds selected by AL prediction 1 (Glide)
│   │   │   ├── pred1_alscreen.csv       # ALScreen predicted score vs Glide gscore, prediction 1
│   │   │   ├── pred2.csv                # compounds selected by AL prediction 2 (Glide)
│   │   │   ├── pred2_alscreen.csv       # ALScreen predicted score vs Glide gscore, prediction 2
│   │   │   └── reference.csv            # reference compounds docked with the same setup
│   │   └── bias1/ ...
│   ├── IL2/{unbias,bias1,bias2}/ ...
│   └── IL36g/{unbias,bias1,bias2}/ ...
├── interaction_fingerprints/
│   └── {IL1b,IL2,IL36g}/
│       ├── receptor.pdb                 # prepared receptor structure
│       └── poses.sdf                    # best Glide pose of the selected compounds
├── reference_compounds/
│   └── {IL1b,IL2,IL36g}_reference.csv   # known binders from the literature
└── experimental_hits/
    └── <target>_round<N>.csv            # compounds tested by SPR in each experimental round
```

## File formats

| Files | Columns |
|---|---|
| `initial.csv`, `pred1.csv`, `pred2.csv`, `reference.csv` | `ID`, `SMILES`, `glide_gscore` (best pose per compound) |
| `pred1_alscreen.csv`, `pred2_alscreen.csv` | `ID`, `SMILES`, `dockscore` (ALScreen prediction), `glide_gscore` |
| `<target>_reference.csv` | `ID`, `DB` (source), `SMILES`, `Affinity` (nM), `Affinity_Type` (IC50/KD) |
| `<target>_round<N>.csv` | `ID`, `SMILES`, `KD_uM` (SPR KD in µM, empty if not measured) |

A `glide_gscore` of 10000 marks a compound that Glide failed to dock; these are discarded by the scripts.

## Bias conditions

| Target (PDB) | `unbias` | `bias1` | `bias2` |
|---|---|---|---|
| IL-1β (8C3U) | no constraints | hydrophobic constraint | – |
| IL-2 (1PY2) | no constraints | H-bond with Glu62 | hydrophobic + H-bond with Glu62 |
| IL-36γ (6P9E) | no constraints | hydrophobic constraint | hydrophobic + H-bond with residue 163 |
