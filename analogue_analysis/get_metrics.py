"""Summarise supplied molecular properties, or calculate them from SMILES."""
import argparse
from pathlib import Path
import pandas as pd
from analogue_common import read_table, output_stem


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', required=True)
    p.add_argument('--pattern', default='*.csv')
    p.add_argument('--columns', nargs='+', default=['mw', 'logp', 'tpsa'])
    p.add_argument('--compute-from-smiles', action='store_true')
    p.add_argument('--smiles-column', default='SMILES')
    p.add_argument('-o', '--output', required=True, help='Summary CSV path')
    a = p.parse_args()
    files = sorted(Path(a.directory).glob(a.pattern))
    files = [f for f in files if f.is_file() and f.resolve() != Path(a.output).resolve()]
    if not files:
        p.error('No matching input files')
    rows = []
    for path in files:
        if a.compute_from_smiles:
            from rdkit import Chem
            from rdkit.Chem import Descriptors
            df = read_table(path, [a.smiles_column])
            properties = []
            for i, smi in enumerate(df[a.smiles_column]):
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    raise ValueError(f'{path}: invalid SMILES on data row {i + 1}')
                properties.append([Descriptors.MolWt(mol), Descriptors.MolLogP(mol), Descriptors.TPSA(mol)])
            df = pd.DataFrame(properties, columns=['mw', 'logp', 'tpsa'])
            columns = ['mw', 'logp', 'tpsa']
        else:
            columns = a.columns
            df = read_table(path, columns, numeric=columns)
        for col in columns:
            values = df[col]
            rows.append(dict(file=str(path), property=col, n=len(values), mean=values.mean(),
                             max=values.max(), min=values.min(), std=values.std()))
    summary = pd.DataFrame(rows)
    summary.to_csv(output_stem(a.output), index=False)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
