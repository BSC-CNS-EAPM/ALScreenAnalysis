"""Compare PELE metrics against an explicitly chosen reference hit."""
import argparse
from pathlib import Path
import pandas as pd
from analogue_common import read_table, single_hit, output_stem


def parse_names(df, constraints, molecule_column='molecule'):
    """Use longest exact suffix; preserve IDs and compound constraints."""
    labels = sorted(set(constraints), key=lambda c: (-len(c), c))
    records = []
    for name in df[molecule_column]:
        matches = [c for c in labels if name.endswith('_' + c)]
        if not matches:
            raise ValueError(f'No declared constraint suffix matches {name!r}')
        label = matches[0]
        identifier = name[:-(len(label) + 1)]
        if not identifier:
            raise ValueError(f'Missing molecule ID in {name!r}')
        records.append((identifier, label))
    result = df.copy()
    result['ID'] = [r[0] for r in records]
    result['constraint'] = [r[1] for r in records]
    if result.duplicated(['ID', 'constraint']).any():
        raise ValueError('Duplicate molecule/constraint combinations')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pele_df', '--pele-csv', dest='pele_csv', required=True)
    p.add_argument('--hit', required=True)
    p.add_argument('--const', '--reference-constraint', dest='reference', required=True)
    p.add_argument('--constraints', nargs='+', required=True, help='All suffix labels in the PELE file')
    p.add_argument('--metrics', nargs='+', default=['all_BFE', 'lowBFE_minBE'])
    p.add_argument('--molecule-column', default='molecule')
    p.add_argument('--output-prefix', required=True)
    a = p.parse_args()
    if a.reference not in a.constraints:
        p.error('Reference constraint must be in --constraints')
    df = read_table(a.pele_csv, [a.molecule_column, *a.metrics], a.metrics, [a.molecule_column])
    df = parse_names(df, a.constraints, a.molecule_column)
    ref = single_hit(df.loc[df['constraint'].eq(a.reference)], a.hit, 'ID', a.pele_csv)
    candidates = df.loc[df['ID'].ne(a.hit)].copy()
    rows = []
    stem = output_stem(a.output_prefix)
    for metric in a.metrics:
        candidates[f'below_{metric}'] = candidates[metric].lt(ref[metric])
        for constraint, group in candidates.groupby('constraint'):
            better = group.loc[group[f'below_{metric}']]
            rows.append(dict(metric=metric, constraint=constraint, threshold=ref[metric],
                             n_rows=len(better), n_unique=better['ID'].nunique()))
        better = candidates.loc[candidates[f'below_{metric}']]
        rows.append(dict(metric=metric, constraint='__ALL__', threshold=ref[metric],
                         n_rows=len(better), n_unique=better['ID'].nunique()))
    candidates.to_csv(str(stem) + '_comparisons.csv', index=False)
    summary = pd.DataFrame(rows)
    summary.to_csv(str(stem) + '_summary.csv', index=False)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
