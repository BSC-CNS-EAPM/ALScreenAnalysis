"""Join Glide and PELE by exact molecule and constraint; plot and count improvements."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from analogue_common import manifest, read_table, single_hit, suffix, molecule_name


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--metrics', nargs='+', default=['all_BFE', 'lowBFE_minBE'])
    p.add_argument('--id-column', default='ID')
    p.add_argument('--score-column', default='glide_gscore')
    p.add_argument('--pele-id-column', default='molecule')
    p.add_argument('--reference-constraint', help='Also count unique IDs against this constraint across all data')
    a = p.parse_args()
    config = manifest(a.manifest)
    out = Path(a.output_dir)
    joined, refs, correlations, summaries = [], [], [], []
    for _, row in config.iterrows():
        analogues = read_table(row['analogue_csv'], [a.id_column, a.score_column],
                               [a.score_column], [a.id_column])
        hits = read_table(row['hit_csv'], [a.id_column, a.score_column], [a.score_column], [a.id_column])
        parent = single_hit(hits, row['hit'], a.id_column, row['hit_csv'])
        pele = read_table(row['pele_csv'], [a.pele_id_column, *a.metrics], a.metrics, [a.pele_id_column])
        parent_pele = single_hit(pele, molecule_name(row['hit'], suffix(row)), a.pele_id_column, row['pele_csv'])
        # Exclude the designated parent and every reference compound listed in hit_csv.
        candidates = analogues.loc[~analogues[a.id_column].isin(hits[a.id_column])].copy()
        candidates = candidates[[a.id_column, a.score_column]].rename(
            columns={a.id_column: 'ID', a.score_column: 'glide_gscore'})
        candidates['_key'] = candidates['ID'].map(lambda v: molecule_name(v, suffix(row)))
        combined = candidates.merge(pele[[a.pele_id_column, *a.metrics]].rename(
            columns={a.pele_id_column: '_key'}), on='_key', how='left', validate='one_to_one', indicator=True)
        missing = combined['_merge'].ne('both')
        if missing.any():
            print(f"{row['target']}/{row['hit']}/{row['constraint']}: {missing.sum()} docking rows lack PELE metrics; omitted")
        combined = combined.loc[~missing].drop(columns=['_merge', '_key'])
        for key in ['target', 'hit', 'constraint']:
            combined[key] = row[key]
        joined.append(combined)
        refs.append({'target': row['target'], 'hit': row['hit'], 'constraint': row['constraint'],
                     'glide_gscore': parent[a.score_column], **{m: parent_pele[m] for m in a.metrics}})
    data, references = pd.concat(joined, ignore_index=True), pd.DataFrame(refs)
    for (target, hit), group in references.groupby(['target', 'hit'], sort=False):
        folder = out / target / hit
        folder.mkdir(parents=True, exist_ok=True)
        candidates = data.loc[data['target'].eq(target) & data['hit'].eq(hit)].copy()
        reference_constraint = a.reference_constraint
        if reference_constraint and reference_constraint not in set(group['constraint']):
            raise ValueError(f'{target}/{hit}: reference constraint {reference_constraint!r} absent')
        for metric in a.metrics:
            candidates[f'below_{metric}'] = False
            fig, ax = plt.subplots(figsize=(8, 6))
            for i, (_, ref) in enumerate(group.iterrows()):
                constraint = ref['constraint']
                points = candidates.loc[candidates['constraint'].eq(constraint)]
                color = plt.get_cmap('tab20')(i % 20)
                ax.scatter(points['glide_gscore'], points[metric], color=color, alpha=0.65, label=constraint)
                ax.scatter([ref['glide_gscore']], [ref[metric]], color=color, marker='*', s=170,
                           edgecolors='black', label=f'Hit ({constraint})')
                ax.axvline(ref['glide_gscore'], color=color, ls='--', lw=0.8)
                ax.axhline(ref[metric], color=color, ls='--', lw=0.8)
                valid = len(points) >= 2 and points['glide_gscore'].nunique() > 1 and points[metric].nunique() > 1
                r, pv = pearsonr(points['glide_gscore'], points[metric]) if valid else (np.nan, np.nan)
                correlations.append(dict(target=target, hit=hit, constraint=constraint, metric=metric,
                                         n=len(points), Pearson_r=r, Pearson_pvalue=pv))
                flag = f'below_{metric}'
                candidates.loc[points.index, flag] = points[metric].lt(ref[metric]).to_numpy()
                better = points.loc[points[metric].lt(ref[metric])]
                summaries.append(dict(target=target, hit=hit, metric=metric, constraint=constraint,
                                      reference=constraint, threshold=ref[metric], n_rows=len(better),
                                      n_unique=better['ID'].nunique()))
            passed = candidates.loc[candidates[f'below_{metric}'].eq(True)]
            summaries.append(dict(target=target, hit=hit, metric=metric, constraint='__ALL__',
                                  reference='matched constraint', threshold=np.nan,
                                  n_rows=len(passed), n_unique=passed['ID'].nunique()))
            if reference_constraint:
                ref = group.loc[group['constraint'].eq(reference_constraint)].iloc[0]
                better = candidates.loc[candidates[metric].lt(ref[metric])]
                summaries.append(dict(target=target, hit=hit, metric=metric, constraint='__ALL__',
                                      reference=reference_constraint, threshold=ref[metric],
                                      n_rows=len(better), n_unique=better['ID'].nunique()))
            ax.set(xlabel=a.score_column, ylabel=metric, title=f'{target} · {hit}')
            ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(folder / f'{metric}_vs_glide.pdf')
            plt.close(fig)
        candidates.to_csv(folder / 'joined_metrics.csv', index=False)
    pd.DataFrame(correlations).to_csv(out / 'correlations.csv', index=False)
    pd.DataFrame(summaries).to_csv(out / 'improvement_counts.csv', index=False)
    references.to_csv(out / 'reference_metrics.csv', index=False)


if __name__ == '__main__':
    main()
