import pandas as pd


def gscore_threshold(dfs_db, gscore_thresh):
    for db_csv in dfs_db:
        for gscore in gscore_thresh:
            gscore = float(gscore)
            counts = db_csv[db_csv['glide_gscore'] < gscore].shape[0]
            print(f"Number of compounds below {gscore}: {counts}")
        print('----------------------------')


if __name__ == '__main__':
    import argparse
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
    args = parser.parse_args()

    # Load arguments
    dfs_db = []
    for db_csv in args.db_csvs:
        dfs_db.append(pd.read_csv(db_csv))

    gscore_thresh = args.gscore_thresh

    gscore_threshold(dfs_db, gscore_thresh)
