import time
import pandas as pd


def get_best_glide_docking_pose(path_to_csv, to_csv=True):
    """ It returns a df with the best poses for each ligand,
        based on the glide gscore """

    df = pd.read_csv(path_to_csv, on_bad_lines='skip')
    _df = df.dropna(subset=['r_i_glide_gscore'])
    idx = _df.groupby('title')['r_i_glide_gscore'].idxmin()
    df_best = _df.loc[idx]
    if to_csv:
        df_best.to_csv(path_to_csv.replace('.csv', '_best.csv'), index=False)

    df_best_simple = df_best.rename(columns={'title': 'id',
                                             'r_i_glide_gscore': 'glide_gscore',
                                             'SMILES': 'smiles'})
    df_best_simple = df_best_simple[['id', 'smiles', 'glide_gscore']]

    if to_csv:
        df_best_simple.to_csv(path_to_csv.replace('.csv', '_best_simple.csv'),
                              index=False)

    return df_best


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Load a glide csv output\
                                     file. Filter by best pose (gscore)\
                                     and generate two csv, one with all\
                                     glide parameters and another\
                                     simplified')
    requiredArguments = parser.add_argument_group('required arguments')
    requiredArguments.add_argument('--csv', help='glide output csv',
                                   required=True)
    args = parser.parse_args()

    glide_csv = args.csv

    df = get_best_glide_docking_pose(glide_csv)
