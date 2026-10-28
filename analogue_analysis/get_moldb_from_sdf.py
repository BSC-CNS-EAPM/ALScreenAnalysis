#!/usr/bin/env python3

from MolecularAnalysis.moldb import MolDB
import os
import sys


def get_moldb_from_sdf(sdfs, outname):
    for i, sdf in enumerate(sdfs):
        moldbobj = MolDB(sdfDB=sdf, paramaters=True)
        moldbobj.saveToPickle(outname[i])

    for f in os.listdir('./'):
        if f.endswith('.pickle'):
            os.rename(f, f.replace('.pickle', '.moldb'))
    """
    outpaths = []
    for out in outname:
        outpath = os.path.dirname(out) + '/'
        outpaths.append(outpath)
    for outpath in outpaths:
        for f in os.listdir(outpath):
            if f.endswith('.pickle'):
                os.rename(outpath + f,
                          (outpath + f).replace('.pickle', '.moldb'))
    """


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='This script outputs a\
    moldb object from an input sdf file.')
    requiredArguments = parser.add_argument_group('required arguments')
    requiredArguments.add_argument('--sdfs',
                                   help='List of sdf files',
                                   nargs='+',
                                   required=True)
    requiredArguments.add_argument('-o',
                                   '--outname',
                                   nargs='+',
                                   help='Output name for the files',
                                   required=True)
    args = parser.parse_args()

    # Checks
    if len(args.outname) != len(args.sdfs):
        raise ValueError('Length of output names must be the same as\
                          the length of input sdf files')

    # Load Arguments
    sdfs = args.sdfs
    outname = args.outname

    get_moldb_from_sdf(sdfs, outname)
