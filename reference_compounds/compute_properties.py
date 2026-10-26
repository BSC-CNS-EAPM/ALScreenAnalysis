"""Molecular properties of a set of compounds, from a csv with their SMILES and affinity.

Takes one csv and writes the same table with these columns added:
    KD_uM             binding affinity in µM, converted from the affinity column
    MW                molecular weight (Da)
    HeavyAtoms        number of heavy atoms
    DeltaG            binding free energy, RT·ln(KD) (kcal/mol)
    LigandEfficiency  |DeltaG| / HeavyAtoms
    linearity_ratio   graph diameter / (HeavyAtoms - 1)

Example:
    python compute_properties.py --input ../data/reference_compounds/IL1b_reference.csv
    python compute_properties.py --input ../data/experimental_hits/IL1b_round1.csv \
        --affinity_column KD_uM --affinity_units uM
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

OUT = Path(__file__).resolve().parent / "results"

RT = 0.001987 * 298.15  # kcal/mol
TO_uM = {"M": 1e6, "mM": 1e3, "uM": 1.0, "nM": 1e-3, "pM": 1e-6}


def properties(smiles, kd_uM):
    """Structural properties of one molecule, plus the affinity-derived ones when there is a KD."""
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return pd.Series({"MW": np.nan, "HeavyAtoms": np.nan, "DeltaG": np.nan,
                          "LigandEfficiency": np.nan, "linearity_ratio": np.nan})
        
    heavy_atoms = mol.GetNumHeavyAtoms()
    diameter = Chem.GetDistanceMatrix(Chem.RemoveHs(mol)).max()
    delta_g = RT * np.log(kd_uM * 1e-6) if pd.notna(kd_uM) and kd_uM > 0 else np.nan
    
    return pd.Series({"MW": Descriptors.MolWt(mol), "HeavyAtoms": heavy_atoms, "DeltaG": delta_g,
                      "LigandEfficiency": abs(delta_g / heavy_atoms) if pd.notna(delta_g) else np.nan,
                      "linearity_ratio": round(diameter / (heavy_atoms - 1), 3) if heavy_atoms > 1 else 1.0})


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--input", required=True, type=Path, help="csv with the compounds")
parser.add_argument("--output", type=Path, help="output csv (default: results/<input name>_properties.csv)")
parser.add_argument("--smiles_column", default="SMILES")
parser.add_argument("--affinity_column", default="Affinity", help="column holding the measured affinity")
parser.add_argument("--affinity_units", default="nM", choices=list(TO_uM), help="units of that column")
args = parser.parse_args()

if not args.input.is_file():
    raise FileNotFoundError(f"{args.input} does not exist")

df = pd.read_csv(args.input)
for column in (args.smiles_column, args.affinity_column):
    if column not in df:
        raise ValueError(f"{args.input} has no column '{column}'; its columns are {list(df)}")

df["KD_uM"] = df[args.affinity_column] * TO_uM[args.affinity_units]

df = df.join(df.apply(lambda row: properties(row[args.smiles_column], row["KD_uM"]), axis=1))

n_invalid = int(df["MW"].isna().sum())
if n_invalid:
    print(f"warning: {n_invalid} compound(s) could not be read from their SMILES")

output = args.output or OUT / f"{args.input.stem}_properties.csv"
output.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(output, index=False)
print(f"{args.input.stem}: {len(df)} compounds saved to {output}")
