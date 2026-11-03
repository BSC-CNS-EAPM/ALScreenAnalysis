import pandas as pd


def gscore_threshold(dfs_db, gscore_thresh):
    for db_csv in dfs_db:
        for gscore in gscore_thresh:
            gscore = float(gscore)
            counts = db_csv[db_csv['glide_gscore'] < gscore].shape[0]
            print(f"Number of compounds below {gscore}: {counts}")
        print('----------------------------')


if __name__ == '__main__':
    # Run from ALScreen_processing/ (shell commands):
    # python glide_score_threshold.py --db_csvs ../data/ALScreen/IL1b/unbias/initial.csv ../data/ALScreen/IL1b/unbias/pred1.csv ../data/ALScreen/IL1b/unbias/pred2.csv --gscore_thresh -7 -7.5 -8 -8.5 -9 -9.5 -10 -10.5

    import argparse
    import sys
    from pathlib import Path

    data = Path(__file__).resolve().parents[1] / 'data' / 'ALScreen' / 'IL1b' / 'unbias'
    # Executed automatically by python glide_score_threshold.py; edit these arguments to change the run.
    default_run = [
        '--db_csvs', str(data / 'initial.csv'), str(data / 'pred1.csv'), str(data / 'pred2.csv'),
        '--gscore_thresh', '-7', '-7.5', '-8', '-8.5', '-9', '-9.5', '-10', '-10.5',
    ]
    parser = argparse.ArgumentParser(description='Obtain the number of\
    compounds below a gscore value threshold')
    requiredArguments = parser.add_argument_group('required arguments')
    requiredArguments.add_argument('--db_csvs',
                                   nargs='+',
                                   help='List of final docking CSVs',
                                   required=True)
    requiredArguments.add_argument('--gscore_thresh',
                                   nargs='+',
                                   help='List of glide gscore thresholds',
                                   required=True)
    args = parser.parse_args(default_run if len(sys.argv) == 1 else None)

    # Load arguments
    dfs_db = []
    for db_csv in args.db_csvs:
        dfs_db.append(pd.read_csv(db_csv))

    gscore_thresh = args.gscore_thresh

    gscore_threshold(dfs_db, gscore_thresh)
