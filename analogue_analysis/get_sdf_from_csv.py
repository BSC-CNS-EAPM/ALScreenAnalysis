from MolecularAnalysis.moldb import Mol, MolDB
import pandas as pd


def generated_sdf_from_csv(df_csv, idt, smi, params, sim_filt, output):
    mols = []
    for ind, row in df_csv.iterrows():
        try:
            ident = row[idt]
            smile = row[smi]
            molobj = Mol(smile=smile, name=str(ident))
            for param in params:
                value = row[param]
                molobj.molrdkit.SetProp(str(param), str(value))
            mols.append(molobj)
        except Exception as e:
            print(f"Error processing row {ind}: {row}")
            print(f"Exception: {e}")
            exit()

    db = MolDB(molList=mols, verbose=True)

    if sim_filt is not False:
        db.filterSimilarity(verbose=True)

    db.saveToSdf(output)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='This script outputs an sdf\
    from an input csv file.')
    requiredArguments = parser.add_argument_group('required arguments')
    requiredArguments.add_argument('--csv', help='input csv file',
                                   required=True)
    requiredArguments.add_argument('--id', help='Add ID to the sdf,\
                        insert ID col name', type=str, required=True)
    requiredArguments.add_argument('--smi', help='Add SMILE to the sdf,\
                         insert SMILE col name', type=str, required=True)
    parser.add_argument('--sim_filt', help='Apply similarity filter (bool))',
                        default=False)
    parser.add_argument('--params', nargs='+',
                        help='List of parameters to add to the sdf',
                        default=None)
    parser.add_argument('-o', help='Output file name')
    parser.add_argument('-d',
                        help='Separator of csv columns',
                        default=',')
    args = parser.parse_args()

    delimiter = args.d
    # Check if the columns exist in the CSV file:
    csv_columns = pd.read_csv(args.csv, sep=delimiter).columns.tolist()
    required_columns = [args.id] + [args.smi] + (args.params if args.params else [])
    for column in filter(None, required_columns):
        if column not in csv_columns:
            raise ValueError('Column "{}" must exist in the CSV file'.format(column))

    # Load CSV as a dataframe
    df_csv = pd.read_csv(args.csv, sep=delimiter)

    # Get args
    idt = str(args.id)
    smi = args.smi
    sim_filt = args.sim_filt
    params = args.params if args.params else []
    output = args.o if args.o else args.csv.split('.')[0:-1]

    generated_sdf_from_csv(df_csv, idt, smi, params, sim_filt, output)
