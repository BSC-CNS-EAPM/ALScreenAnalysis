from MolecularAnalysis.moldb import MolDB
from MolecularAnalysis.moldb import intersectMolDBs
import os


def intersect_moldbs(moldb1, moldb2, outname, sdf, sim):
    intersect = intersectMolDBs(moldb1, moldb2, simt=sim)
    print(intersect.size)

    if sim:
        intersect.filterSimilarity(simt=sim, verbose=False)

    if sdf:
        intersect.saveToSdf(outname)
    if intersect.size > 0:
        intersect.saveToPickle(outname)
    for f in os.listdir('.'):
        if f.endswith('.pickle'):
            os.rename(f, f.replace('.pickle', '.moldb'))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='This script finds the\
    intersection between two moldb objects.')
    requiredArguments = parser.add_argument_group('required arguments')
    requiredArguments.add_argument('--moldb1',
                                   help='First moldb file',
                                   required=True)
    requiredArguments.add_argument('--moldb2',
                                   help='Second moldb file',
                                   required=True)
    requiredArguments.add_argument('-o',
                                   '--outname',
                                   help='Output file name',
                                   required=True)
    parser.add_argument('--sim',
                        help='Apply similarity filter to the intersect',
                        default=1)
    parser.add_argument('--sdf',
                        help='Save the intersect to an sdf file',
                        action='store_true')

    args = parser.parse_args()

    # Load arguments
    moldb1 = MolDB(molDB=args.moldb1)
    moldb2 = MolDB(molDB=args.moldb2)
    outname = args.outname
    sim = float(args.sim)
    sdf = args.sdf

    intersect_moldbs(moldb1, moldb2, outname, sdf, sim)
