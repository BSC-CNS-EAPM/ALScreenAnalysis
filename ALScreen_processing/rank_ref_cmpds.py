import pandas as pd
from pathlib import Path


def main(ref_csv, db_csv, output_prefix=None):

    df_ref = pd.read_csv(ref_csv)
    df_db = pd.read_csv(db_csv)
    df_ref = df_ref[df_ref['glide_gscore'] < 50]
    df_db = df_db[df_db['glide_gscore'] < 50]

    print('Min DB gscore: %.4f' % min(df_db['glide_gscore'].to_numpy()))
    print('Max DB gscore: %.4f' % max(df_db['glide_gscore'].to_numpy()))

    print('\n')

    df_total = pd.concat([df_db, df_ref])
    df_sort = df_total.sort_values(by='glide_gscore', ignore_index=True)
    df_ref_sort = df_ref.sort_values(by='glide_gscore', ignore_index=True)

    lines = []
    for compound in df_ref_sort['ID'].tolist():
        print(compound)
        index = df_sort.index[df_sort['ID'] == compound].tolist()[0]+1
        print('%d/%d' % (index, len(df_total)))
        lines.append(f'{compound}\t{index}/{len(df_total)}')

    prefix = Path(output_prefix) if output_prefix is not None else Path(__file__).resolve().parent / 'results' / 'IL1b_unbias_pred2'
    output = Path(str(prefix) + '_ref_rank.txt')
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text('Reference compound\tRank/total\n' + '\n'.join(lines) + '\n', encoding='utf-8')
    print(f'Reference ranks saved to {output}')


if __name__ == '__main__':
    # Run from ALScreen_processing/ (shell commands):
    # python rank_ref_cmpds.py --db_csv ../data/ALScreen/IL1b/unbias/pred2.csv --ref_csv ../data/ALScreen/IL1b/unbias/reference.csv

    import argparse
    import sys
    from pathlib import Path

    data = Path(__file__).resolve().parents[1] / 'data' / 'ALScreen' / 'IL1b' / 'unbias'
    # Executed automatically by python rank_ref_cmpds.py; edit these arguments to change the run.
    default_run = ['--db_csv', str(data / 'pred2.csv'), '--ref_csv', str(data / 'reference.csv')]
    parser = argparse.ArgumentParser(description='rank reference compounds\
                                     in the context of db compounds')
    requiredArguments = parser.add_argument_group('required arguments')
    requiredArguments.add_argument('--ref_csv',
                                   help='glide simplified output csv for\
                                   reference compounds', required=True)
    requiredArguments.add_argument('--db_csv',
                                   help='glide simplified output csv for\
                                   DB compounds', required=True)
    parser.add_argument('-o', '--output-prefix',
                        default=str(Path(__file__).resolve().parent / 'results' / 'IL1b_unbias_pred2'),
                        help='Output prefix; writes <output>_ref_rank.txt')
    args = parser.parse_args(default_run if len(sys.argv) == 1 else None)

    ref_csv = args.ref_csv
    db_csv = args.db_csv

    main(ref_csv, db_csv, args.output_prefix)
