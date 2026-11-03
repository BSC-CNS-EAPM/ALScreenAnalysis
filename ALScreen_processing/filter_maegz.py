# $SCHRODINGER/run filt_maegz.py
from schrodinger import structure
import numpy as np
# from rdkit import Chem
from schrodinger.analysis.transformations import structure_to_rdkit_smiles

if __name__ == '__main__':
    
    # Run from ALScreen_processing/ (shell commands):
    # Missing input: data/ contains no Glide pose-viewer .maegz files.
    # Requires the original pose viewer and a licensed Schrodinger installation.
    # With that additional input (replace /path/to/glide_pv.maegz):
    # $SCHRODINGER/run filter_maegz.py --maegz /path/to/glide_pv.maegz -o filtered --metric_threshold -8 --neutral
    # The supplied SDF/PDB files are not pose-viewer .maegz files.

    import argparse
    parser = argparse.ArgumentParser(description='Filter maegz poses\
                                     by the best pose for each ligand\
                                     and if needed by gscore')
    requiredArguments = parser.add_argument_group('required arguments')
    requiredArguments.add_argument('--maegz',
                                   help='glide maegz to be filter',
                                   required=True)
    requiredArguments.add_argument('-o',
                                   help='output name for the filtrated maegz',
                                   required=True)
    requiredArguments.add_argument('--metric',
                                   help='metric to filter the maegz\
                                   (default: gscore). Options: gscore, docking_score',
                                   default='gscore')
    parser.add_argument('--metric_threshold',
                        help='metric threshold\
                        (molecules above the trheshold will be discarded)',
                        default=None)
    parser.add_argument('--neutral',
                        help='neutral threshold\
                        (molecules ionized will be discarded)',
                        default=False,
                        action='store_true')
    parser.add_argument('--mw',
                        help='molecular weight filter e.g 400 minor.\
                        (molecules above 400 will be discarded \
                        would be discarded). Format [mw_threshold, \'minor\'/\'mayor\']',
                        nargs='+',
                        default=None),

    args = parser.parse_args()
    maegz = args.maegz
    out = args.o
    metric = args.metric
    metric_threshold = args.metric_threshold
    if metric == 'gscore':
        metric = 'r_i_glide_gscore'
    elif metric == 'docking_score':
        metric = 'r_i_docking_score'
    if metric_threshold:
        metric_threshold = float(metric_threshold)
    else:
        metric_threshold = np.inf
    neutral = args.neutral
    mw_threshold = args.mw
    print(mw_threshold)
    if mw_threshold:
        if len(mw_threshold) != 2:
            raise ValueError('Error: --mw must have two arguments: first a real number\
            , second a string with only two options: minor or major')
        try:
            mw_threshold[0] = float(mw_threshold[0])
        except Exception as typ:
            print('Error: --mw first argument must be a real_number')
            print(f'Error: {typ}')
            exit()
        if mw_threshold[1] != 'minor' and mw_threshold[1] != 'major':
            raise ValueError('Error: --mw second argument must be either \'minor\'\
            or \'major\'')

    maegz = structure.StructureReader(maegz)
    ligands_dic = {}
    count_str = 0
    for i, struc in enumerate(maegz):
        print(i)
        title = struc.property['s_m_title']
        props = struc.property.keys()
        if metric not in props:
            # prot structure
            target = struc
            count_str += 1
            continue
        else:
            # small molecule
            try:
                smile = struc.property['s_sd_SMILE']
            except:
                smile = structure_to_rdkit_smiles(struc)
            metric_value = float(struc.property[metric])
            if metric_value > metric_threshold:
                continue
            if neutral:
                aux = False
                for atom in struc.atom:
                    charge = atom.formal_charge
                    if charge != 0:
                        aux = True
                        break
                if aux:
                    continue
            if mw_threshold:
                mw = struc.total_weight
                if mw_threshold[1] == 'major':
                    if mw < mw_threshold[0]:
                        continue
                elif mw_threshold[1] == 'minor':
                    if mw > mw_threshold[0]:
                        continue
            if title in ligands_dic.keys():
                if metric_value < ligands_dic[title][0]:
                    ligands_dic[title] = (metric_value, struc)
            else:
                ligands_dic[title] = (metric_value, struc)

    print(count_str)
    print(len(list(ligands_dic.keys())))
    # Save the final filtrated maegz
    with structure.StructureWriter(out+'.maegz') as writer:
        writer.append(target)
        for title in ligands_dic.keys():
            writer.append(ligands_dic[title][1])
