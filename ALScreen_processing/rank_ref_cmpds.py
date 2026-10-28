import pandas as pd


def main(ref_csv, db_csv):

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

    for compound in df_ref_sort['ID'].tolist():
        print(compound)
        index = df_sort.index[df_sort['ID'] == compound].tolist()[0]+1
        print('%d/%d' % (index, len(df_total)))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='rank reference compounds\
                                     in the context of db compounds')
    requiredArguments = parser.add_argument_group('required arguments')
    requiredArguments.add_argument('--ref_csv',
                                   help='glide simplified output csv for\
                                   reference compounds', required=True)
    requiredArguments.add_argument('--db_csv',
                                   help='glide simplified output csv for\
                                   DB compounds', required=True)
    args = parser.parse_args()

    ref_csv = args.ref_csv
    db_csv = args.db_csv

    main(ref_csv, db_csv)
