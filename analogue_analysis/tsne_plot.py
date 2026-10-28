"""Plot chemical space for a CSV or SDF, optionally highlighting a reference hit."""
import argparse
import inspect
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.manifold import TSNE
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analogue_common import read_table, single_hit, output_stem


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', required=True, help='CSV or SDF')
    p.add_argument('--smiles-column', default='smiles')
    p.add_argument('--group-column', default='space')
    p.add_argument('--delimiter', default=',')
    p.add_argument('--experimental-csv')
    p.add_argument('--hit')
    p.add_argument('--experimental-id-column', default='ID')
    p.add_argument('--experimental-smiles-column', default='SMILES')
    p.add_argument('--experimental-delimiter', default=';')
    p.add_argument('--radius', type=int, default=4)
    p.add_argument('--n-bits', type=int, default=2048)
    p.add_argument('--perplexity', type=float, default=50)
    p.add_argument('--iterations', type=int, default=5000)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--metric', choices=['euclidean', 'jaccard'], default='euclidean')
    p.add_argument('-o', '--output-prefix', required=True)
    a = p.parse_args()
    if bool(a.experimental_csv) != bool(a.hit):
        p.error('Supply both --experimental-csv and --hit, or neither')
    if a.radius < 0 or a.n_bits < 1 or a.iterations < 250:
        p.error('radius >= 0, n-bits > 0 and iterations >= 250 are required')
    if Path(a.input).suffix.lower() == '.sdf':
        records = []
        for i, mol in enumerate(Chem.SDMolSupplier(a.input)):
            if mol is None:
                raise ValueError(f'Invalid SDF record {i + 1}')
            if not mol.HasProp(a.group_column):
                raise ValueError(f'SDF record {i + 1} lacks property {a.group_column!r}')
            records.append({a.smiles_column: Chem.MolToSmiles(mol), a.group_column: mol.GetProp(a.group_column)})
        df = pd.DataFrame(records, columns=[a.smiles_column, a.group_column])
    else:
        df = pd.read_csv(a.input, sep=a.delimiter, dtype=str, keep_default_na=False)
    for col in [a.smiles_column, a.group_column]:
        if col not in df or df[col].str.strip().eq('').any():
            raise ValueError(f'Missing/blank required column {col!r}')
    if df.empty:
        raise ValueError('No analogue molecules')
    smiles = df[a.smiles_column].tolist()
    if a.experimental_csv:
        exp = pd.read_csv(a.experimental_csv, sep=a.experimental_delimiter, dtype=str, keep_default_na=False)
        for col in [a.experimental_id_column, a.experimental_smiles_column]:
            if col not in exp:
                raise ValueError(f'Experimental CSV lacks {col!r}')
        ref = single_hit(exp, a.hit, a.experimental_id_column, a.experimental_csv)
        smiles.append(ref[a.experimental_smiles_column])
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=a.radius, fpSize=a.n_bits)
    vectors = []
    for i, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f'Invalid SMILES at combined row {i + 1}: {smi!r}')
        vectors.append(gen.GetFingerprintAsNumPy(mol))
    if not 0 < a.perplexity < len(vectors):
        p.error(f'Perplexity must be positive and less than total sample count ({len(vectors)})')
    kwargs = dict(n_components=2, perplexity=a.perplexity, random_state=a.seed,
                  metric=a.metric, init='random', learning_rate=200.0)
    # Support both scikit-learn parameter spellings without fixing one release.
    key = 'max_iter' if 'max_iter' in inspect.signature(TSNE).parameters else 'n_iter'
    kwargs[key] = a.iterations
    embeddings = TSNE(**kwargs).fit_transform(np.asarray(vectors, dtype=bool if a.metric == 'jaccard' else np.float32))
    df = df.copy()
    df['tsne_x'], df['tsne_y'], df['is_reference'] = embeddings[:len(df), 0], embeddings[:len(df), 1], False
    fig, ax = plt.subplots(figsize=(10, 8))
    for label, group in df.groupby(a.group_column, sort=False):
        ax.scatter(group['tsne_x'], group['tsne_y'], alpha=0.6, marker='.',
                   label=f'{label}: {100 * len(group) / len(df):.3f}%')
    if a.experimental_csv:
        ax.scatter(*embeddings[-1], color='red', marker='X', s=100, label=f'Experimental hit: {a.hit}')
        refrow = {a.smiles_column: smiles[-1], a.group_column: '__reference__',
                  'tsne_x': embeddings[-1, 0], 'tsne_y': embeddings[-1, 1], 'is_reference': True}
        df = pd.concat([df, pd.DataFrame([refrow])], ignore_index=True)
    ax.set(title='t-SNE of molecules grouped by database', xlabel='t-SNE 1', ylabel='t-SNE 2')
    ax.legend(title='Databases (analogue percentages)')
    fig.tight_layout()
    stem = output_stem(a.output_prefix)
    fig.savefig(str(stem) + '.pdf')
    plt.close(fig)
    df.to_csv(str(stem) + '_embeddings.csv', index=False)


if __name__ == '__main__':
    main()
