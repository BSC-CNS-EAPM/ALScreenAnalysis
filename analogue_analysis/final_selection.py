"""Select analogues passing Glide, visual-inspection and PELE criteria."""
import argparse
import shutil
from pathlib import Path
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from analogue_common import (manifest, read_table, single_hit, suffix,
                             molecule_name, optional_number)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--id-column', default='ID')
    p.add_argument('--smiles-column', default='SMILES')
    p.add_argument('--score-column', default='glide_gscore')
    p.add_argument('--pele-id-column', default='molecule')
    p.add_argument('--pele-metric', default='lowBFE_minBE')
    p.add_argument('--radius', type=int, default=4)
    p.add_argument('--n-bits', type=int, default=2048)
    p.add_argument('--copy-poses', action='store_true')
    a = p.parse_args()
    if a.radius < 0 or a.n_bits < 1:
        p.error('radius must be nonnegative and n-bits positive')
    config = manifest(a.manifest, extra=['visual_dir'])
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=a.radius, fpSize=a.n_bits)
    results, reports, copies = [], [], []
    columns = ['ID', 'SMILES', 'similarity_to_hit', 'glide_gscore',
               'glide_constraint', 'PELE_metric', 'PELE_value']
    for _, row in config.iterrows():
        hit, constraint = row['hit'], row['constraint']
        analogues = read_table(row['analogue_csv'], [a.id_column, a.smiles_column, a.score_column],
                               [a.score_column], [a.id_column])
        hits = read_table(row['hit_csv'], [a.id_column, a.smiles_column, a.score_column],
                          [a.score_column], [a.id_column])
        parent = single_hit(hits, hit, a.id_column, row['hit_csv'])
        pele = read_table(row['pele_csv'], [a.pele_id_column, a.pele_metric],
                          [a.pele_metric], [a.pele_id_column])
        parent_pele = single_hit(pele, molecule_name(hit, suffix(row)), a.pele_id_column, row['pele_csv'])
        mol = Chem.MolFromSmiles(parent[a.smiles_column])
        if mol is None:
            raise ValueError(f'Invalid hit SMILES: {hit}')
        parent_fp = gen.GetFingerprint(mol)
        visual_dir = Path(row['visual_dir'])
        if not visual_dir.is_dir():
            raise FileNotFoundError(visual_dir)
        # Match complete expected names, so HP never matches HP_HB.
        visual_names = {f.stem for f in visual_dir.glob('*.pdb')}
        selected = analogues.loc[analogues[a.id_column].ne(hit) &
                                 analogues[a.score_column].lt(parent[a.score_column])].copy()
        threshold = optional_number(row, 'gscore_threshold')
        if threshold is not None:
            selected = selected.loc[selected[a.score_column].lt(threshold)].copy()
        n_glide = len(selected)
        selected = selected.loc[selected[a.id_column].map(
            lambda v: molecule_name(v, suffix(row, 'visual')) in visual_names)].copy()
        n_visual = len(selected)
        selected['_pele_key'] = selected[a.id_column].map(lambda v: molecule_name(v, suffix(row)))
        selected = selected.merge(pele[[a.pele_id_column, a.pele_metric]].rename(
            columns={a.pele_id_column: '_pele_key', a.pele_metric: '_pele_value'}),
            on='_pele_key', how='left', validate='one_to_one')
        if selected['_pele_value'].isna().any():
            missing = selected.loc[selected['_pele_value'].isna(), '_pele_key'].tolist()
            raise ValueError(f'Missing PELE metrics for visually accepted Glide candidates: {missing}')
        selected = selected.loc[selected['_pele_value'].lt(parent_pele[a.pele_metric])].copy()
        pele_threshold = optional_number(row, 'pele_threshold')
        if pele_threshold is not None:
            selected = selected.loc[selected['_pele_value'].lt(pele_threshold)]
        outrows = []
        for _, candidate in selected.iterrows():
            identifier = candidate[a.id_column]
            mol = Chem.MolFromSmiles(candidate[a.smiles_column])
            if mol is None:
                raise ValueError(f'Invalid analogue SMILES: {identifier}')
            outrows.append([identifier, candidate[a.smiles_column],
                            DataStructs.TanimotoSimilarity(parent_fp, gen.GetFingerprint(mol)),
                            candidate[a.score_column], constraint, a.pele_metric, candidate['_pele_value']])
            if a.copy_poses:
                if not row.get('pose_dir', ''):
                    raise ValueError('--copy-poses requires pose_dir in every manifest row')
                source = Path(row['pose_dir']) / (molecule_name(identifier, suffix(row, 'pose')) + '.pdb')
                # Candidate IDs become filenames only when copying poses.
                if Path(identifier).name != identifier or '\\' in identifier or identifier in ('.', '..'):
                    raise ValueError(f'Unsafe pose filename: {identifier}')
                if not source.is_file():
                    raise FileNotFoundError(source)
                copies.append((source, row['target'], hit, constraint, identifier))
        results.append((row['target'], hit, constraint, pd.DataFrame(outrows, columns=columns)))
        reports.append({'target': row['target'], 'hit': hit, 'constraint': constraint,
                        'input_rows': len(analogues), 'after_glide': n_glide,
                        'after_visual': n_visual, 'after_pele': len(selected)})
    # Validate the whole run before writing results.
    root = Path(a.output_dir)
    for (target, hit), group in config.groupby(['target', 'hit'], sort=False):
        folder = root / target / hit
        folder.mkdir(parents=True, exist_ok=True)
        frames = []
        for t, h, constraint, df in results:
            if (t, h) == (target, hit):
                df.to_csv(folder / f'{target}_{hit}_{constraint}_analogues.csv', index=False)
                frames.append(df)
        pd.concat(frames, ignore_index=True).to_csv(folder / f'{target}_{hit}_analogues.csv', index=False)
    for source, target, hit, constraint, identifier in copies:
        dest = root / target / hit / 'PELE_poses' / constraint
        dest.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest / f'{identifier}.pdb')
    pd.DataFrame(reports).to_csv(root / 'selection_summary.csv', index=False)
    print(pd.DataFrame(reports).to_string(index=False))


if __name__ == '__main__':
    main()
