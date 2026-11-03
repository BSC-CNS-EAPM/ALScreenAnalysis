import pandas as pd
from pathlib import Path


def score_diff(csv, n, outname):
    csv = csv.rename(columns={'SMILE': 'SMILES'}) if 'SMILES' not in csv else csv
    csv = csv[csv['glide_gscore'] < 50].copy()
    Path(outname).parent.mkdir(parents=True, exist_ok=True)
    
    #Get MAE and RMSE between glide_gscore and dockscore
    mae = (abs(csv['glide_gscore'] - csv['dockscore'])).mean()
    rmse = ((csv['glide_gscore'] - csv['dockscore'])**2).mean()**0.5
    
    print(f"MAE: {mae}")
    print(f"RMSE: {rmse}")
    
    csv['score_diff'] = abs(csv['glide_gscore'] - csv['dockscore'])

    df_sorted = csv.sort_values(by='score_diff')

    lowest_diff = df_sorted.head(n)
    highest_diff = df_sorted.tail(n)

    lowest_diff.to_csv(f"{outname}_top_{n}k.csv", index=False)
    highest_diff.to_csv(f"{outname}_bottom_{n}k.csv", index=False)

    lowest_highest_diff = pd.concat([lowest_diff, highest_diff],
                                    ignore_index=True)
    print(lowest_highest_diff)

    lowest_highest_diff = lowest_highest_diff[['ID', 'SMILES', 'glide_gscore']]

    lowest_highest_diff.to_csv(outname + '_diff.csv', index=False)


if __name__ == '__main__':
    # Run from ALScreen_processing/ (shell commands):
    # python training_score_diff.py --csv ../data/ALScreen/IL1b/unbias/pred2_alscreen.csv -n 7500 -o ../data/ALScreen/IL1b/unbias/train_3_data

    import argparse
    import sys

    here = Path(__file__).resolve().parent
    data = here.parent / 'data' / 'ALScreen' / 'IL1b' / 'unbias'
    # Executed automatically by python training_score_diff.py; edit these arguments to change the run.
    default_run = [
        '--csv', str(data / 'pred2_alscreen.csv'), '-n', '7500',
        '-o', str(here / 'results' / 'IL1b_unbias_train_3_data'),
    ]
    parser = argparse.ArgumentParser(description='This script takes a csv file\
    with the predicted and real gscores and computes the absolute difference\
    between them. Furthermore, it stores the output in a new csv file as well\
    as outputting two more csv files corresponding to the top and bottom best\
    and worse predicted gscores.')
    requiredArguments = parser.add_argument_group('required arguments')
    requiredArguments.add_argument('--csv',
                                   help='Merged csv file',
                                   required=True)
    requiredArguments.add_argument('-n',
                                   type=int,
                                   help='Number of best and worst predicted molecules to extract',
                                   required=True)
    requiredArguments.add_argument('-o',
                                   '--outname',
                                   help='Output name for the files',
                                   required=True)
    args = parser.parse_args(default_run if len(sys.argv) == 1 else None)

    # Load arguments
    csv = pd.read_csv(args.csv)
    n = args.n
    outname = args.outname

    score_diff(csv, n, outname)
