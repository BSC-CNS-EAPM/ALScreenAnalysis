"""Chemotype and scaffold diversity of the compounds selected in each AL round.

Molecules are grouped into chemotypes: connected components of the graph that
links two molecules when the Tanimoto similarity of their Morgan fingerprints
(radius 2, 2048 bits) is >= 0.8. The same grouping is applied to the
Bemis-Murcko scaffolds of the molecules (scaffold chemotypes).

Outputs (results/):
    summary_chemotypes.csv           diversity metrics per target, bias and round
    scaffold_chemotype_evolution.csv scaffold chemotypes retained, lost and new between
                                     rounds (initial -> pred1 -> pred2)
"""
from pathlib import Path
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

DATA = Path(__file__).resolve().parents[1] / "data" / "ALScreen"
OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)

TARGETS = {"IL1b": ["unbias", "bias1"], "IL2": ["unbias", "bias1", "bias2"], "IL36g": ["unbias", "bias1", "bias2"]}
ROUNDS = ["initial", "pred1", "pred2"]
THRESHOLD = 0.8

def build_chemotypes(fps):
    """Connected components of the similarity graph (Tanimoto >= THRESHOLD)."""
    visited = np.zeros(len(fps), dtype=bool)
    chemotypes = []
    for i in range(len(fps)):
        if visited[i]:
            continue
        visited[i] = True
        stack, cluster = [i], []
        while stack:
            idx = stack.pop()
            cluster.append(idx)
            sims = np.array(DataStructs.BulkTanimotoSimilarity(fps[idx], fps))
            new = np.where((sims >= THRESHOLD) & ~visited)[0]
            visited[new] = True
            stack.extend(new)
        chemotypes.append(cluster)
    return chemotypes

def analyse_round(target, bias, rnd):
    df = pd.read_csv(DATA / target / bias / f"{rnd}.csv")
    df = df[df["glide_gscore"] != 10000].dropna(subset=["ID", "SMILES", "glide_gscore"])
    df["mol"] = df["SMILES"].apply(Chem.MolFromSmiles)
    df = df.dropna(subset=["mol"])
    mols = list(df["mol"])
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    scaffolds = [MurckoScaffold.MurckoScaffoldSmiles(mol=m) for m in mols]

    mol_fps = [fpgen.GetFingerprint(m) for m in mols]
    scaf_fps = [fpgen.GetFingerprint(Chem.MolFromSmiles(s)) for s in scaffolds if s]
    mol_chemotypes = build_chemotypes(mol_fps)
    scaf_chemotypes = build_chemotypes(scaf_fps)

    n_top = max(1, int(len(df) * 0.01))
    metrics = {"Target": target, "Bias": bias, "Iteration": rnd,
               "n_molecules": len(mols),
               "n_chemotypes": len(mol_chemotypes),
               "chemotype_diversity": len(mol_chemotypes) / len(mols),
               "redundancy": 1 - len(mol_chemotypes) / len(mols),
               "avg_cluster_size": np.mean([len(c) for c in mol_chemotypes]),
               "max_cluster_size": max(len(c) for c in mol_chemotypes),
               "n_scaffold_chemotypes": len(scaf_chemotypes),
               "scaffold_chemotype_ratio": len(scaf_chemotypes) / len(scaf_fps),
               "unique_scaffolds_exact": len(set(scaffolds)),
               "top1_glide": df["glide_gscore"].nsmallest(n_top).mean()}
    # The first member of each scaffold chemotype is used as its representative
    representatives = [scaf_fps[c[0]] for c in scaf_chemotypes]
    return metrics, representatives

def match_chemotypes(reps_a, reps_b):
    """Greedy one-to-one matching of chemotype representatives (Tanimoto >= THRESHOLD)."""
    matched, retained = set(), 0
    for fp in reps_a:
        sims = DataStructs.BulkTanimotoSimilarity(fp, reps_b)
        candidates = [(s, j) for j, s in enumerate(sims) if s >= THRESHOLD and j not in matched]
        if candidates:
            matched.add(max(candidates)[1])
            retained += 1
    return retained, len(reps_a) - retained, len(reps_b) - retained

jobs = [(t, b, r) for t, biases in TARGETS.items() for b in biases for r in ROUNDS]
results = Parallel(n_jobs=-1, verbose=10)(delayed(analyse_round)(*job) for job in jobs)
pd.DataFrame([m for m, _ in results]).to_csv(OUT / "summary_chemotypes.csv", index=False)

reps = {job: r for job, (_, r) in zip(jobs, results)}
rows = []
for target, biases in TARGETS.items():
    for bias in biases:
        init, pred1, pred2 = (reps[(target, bias, r)] for r in ROUNDS)
        ret1, lost1, new1 = match_chemotypes(init, pred1)
        ret2, lost2, new2 = match_chemotypes(pred1, pred2)
        rows.append({"Target": target, "Bias": bias,
                     "Initial_Scaffolds": len(init), "Pred1_Scaffolds": len(pred1), "Pred2_Scaffolds": len(pred2),
                     "Initial_to_Pred1_Retained": ret1, "Initial_to_Pred1_Retained_%": 100 * ret1 / len(init),
                     "Initial_to_Pred1_Lost": lost1, "Initial_to_Pred1_Lost_%": 100 * lost1 / len(init),
                     "Pred1_New_Scaffolds": new1, "Pred1_New_Scaffolds_%": 100 * new1 / len(pred1),
                     "Pred1_to_Pred2_Retained": ret2, "Pred1_to_Pred2_Retained_%": 100 * ret2 / len(pred1),
                     "Pred1_to_Pred2_Lost": lost2, "Pred1_to_Pred2_Lost_%": 100 * lost2 / len(pred1),
                     "Pred2_New_Scaffolds": new2, "Pred2_New_Scaffolds_%": 100 * new2 / len(pred2)})

pd.DataFrame(rows).round(2).to_csv(OUT / "scaffold_chemotype_evolution.csv", index=False)
