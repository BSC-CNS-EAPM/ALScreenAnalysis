"""Compute the protein-ligand interaction fingerprints (ProLIF) of the docked poses of one target.

Reads the prepared receptor (pdb) and the docked poses (sdf) and saves the fingerprint
table as results/<target>_ifp.pkl.

    python compute_fingerprints.py --target IL1b \
        --pdb ../data/interaction_fingerprints/IL1b/receptor.pdb \
        --sdf ../data/interaction_fingerprints/IL1b/poses.sdf
"""
import argparse
from pathlib import Path
import prolif as plf
from rdkit import Chem

OUT = Path(__file__).resolve().parent / "results"

INTERACTIONS = ["Hydrophobic", "HBDonor", "HBAcceptor", "PiStacking", "Anionic", "Cationic", "CationPi",
                "PiCation", "XBAcceptor", "XBDonor", "MetalAcceptor", "MetalDonor"]

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--target", required=True, help="target name, used to name the output file")
parser.add_argument("--pdb", required=True, type=Path, help="prepared receptor structure (pdb)")
parser.add_argument("--sdf", required=True, type=Path, help="docked poses of the compounds (sdf)")
parser.add_argument("--output", type=Path, help="output pickle (default: results/<target>_ifp.pkl)")
args = parser.parse_args()

for path in (args.pdb, args.sdf):
    if not path.is_file():
        raise FileNotFoundError(f"{path} does not exist")

receptor = Chem.MolFromPDBFile(str(args.pdb), removeHs=False)
if receptor is None:
    raise ValueError(f"RDKit could not read {args.pdb}; check that it is a prepared receptor structure")
protein = plf.Molecule(receptor)
poses = plf.sdf_supplier(str(args.sdf))

fp = plf.Fingerprint(interactions=INTERACTIONS)
fp.run_from_iterable(poses, protein)

df = fp.to_dataframe()
df.index = [mol.GetProp("_Name") for mol in poses]

output = args.output or OUT / f"{args.target}_ifp.pkl"
output.parent.mkdir(parents=True, exist_ok=True)
df.to_pickle(output)
print(f"{args.target}: {df.shape[0]} compounds, {df.shape[1]} residue-interaction pairs saved to {output}")
