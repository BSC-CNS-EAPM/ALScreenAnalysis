import time
import pandas as pd
from pathlib import Path


def get_best_glide_docking_pose(path_to_csv, to_csv=True, output_prefix=None):
    """ It returns a df with the best poses for each ligand,
        based on the glide gscore """

    df = pd.read_csv(path_to_csv, on_bad_lines='skip')
    # Accept both raw Glide exports and the supplied best-pose tables.
    score = 'r_i_glide_gscore' if 'r_i_glide_gscore' in df else 'glide_gscore'
    identifier = next((c for c in ('title', 'ID', 'id') if c in df), None)
    if identifier is None or score not in df:
        raise ValueError('Expected title/ID/id and r_i_glide_gscore/glide_gscore columns')
    _df = df.dropna(subset=[score])
    idx = _df.groupby(identifier)[score].idxmin()
    df_best = _df.loc[idx]
    prefix = Path(output_prefix) if output_prefix is not None else Path(path_to_csv).with_suffix('')
    if to_csv:
        prefix.parent.mkdir(parents=True, exist_ok=True)
        df_best.to_csv(str(prefix) + '_best.csv', index=False)

    df_best_simple = df_best.rename(columns={'title': 'ID', 'id': 'ID',
                                             'r_i_glide_gscore': 'glide_gscore',
                                             'smiles': 'SMILES', 'SMILE': 'SMILES'})
    df_best_simple = df_best_simple[['ID', 'SMILES', 'glide_gscore']]

    if to_csv:
        df_best_simple.to_csv(str(prefix) + '_best_simple.csv',
                              index=False)

    return df_best


if __name__ == '__main__':
    # Run from ALScreen_processing/ (shell commands):
    # python process_glide_csv_output.py --csv ../data/ALScreen/IL1b/unbias/pred2.csv
    # The supplied tables already contain the best pose; raw Glide CSVs are also accepted.

    import argparse
    import sys

    here = Path(__file__).resolve().parent
    data = here.parent / 'data' / 'ALScreen' / 'IL1b' / 'unbias'
    # Executed automatically by python process_glide_csv_output.py; edit these arguments to change the run.
    default_run = [
        '--csv', str(data / 'pred2.csv'),
        '-o', str(here / 'results' / 'IL1b_unbias_pred2'),
    ]
    parser = argparse.ArgumentParser(description='Load a glide csv output\
                                     file. Filter by best pose (gscore)\
                                     and generate two csv, one with all\
                                     glide parameters and another\
                                     simplified')
    requiredArguments = parser.add_argument_group('required arguments')
    requiredArguments.add_argument('--csv', help='glide output csv',
                                   required=True)
    parser.add_argument('-o', '--output-prefix',
                        help='Output prefix (default with --csv: alongside the input)')
    args = parser.parse_args(default_run if len(sys.argv) == 1 else None)

    glide_csv = args.csv

    df = get_best_glide_docking_pose(glide_csv, output_prefix=args.output_prefix)
