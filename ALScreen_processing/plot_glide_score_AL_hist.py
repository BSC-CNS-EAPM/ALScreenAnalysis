import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def generate_colors_from_palette(palette_name, n):
    cmap = plt.get_cmap(palette_name)
    norm = mcolors.Normalize(vmin=0, vmax=n-1)
    scalar_map = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    colors = [scalar_map.to_rgba(i) for i in range(n)]
    return colors


def plot_hists(dfs_dbs, dfs_refs, colors_dbs, colors_refs,
               labels_refs, ref_gscores, labels_dbs, transparency,
               bins, xlim, ylim, frequency, output):

    plt.figure(figsize=(8, 6), dpi=2200)

    for i, df_db in enumerate(dfs_dbs):
        df_db = df_db[df_db['glide_gscore'] < 50]

        values = df_db['glide_gscore'].to_numpy()
        mean_gscore = df_db['glide_gscore'].mean()
        std_gscore = df_db['glide_gscore'].std()

        print(labels_dbs[i], mean_gscore, std_gscore)

        counts, bin_edges, _ = plt.hist(values,
                                        bins=bins,
                                        density=frequency,
                                        color=colors_dbs[i],
                                        alpha=0.45,
                                        edgecolor=colors_dbs[i],
                                        linewidth=0.5,
                                        label=f"{labels_dbs[i]} (mean={mean_gscore:.2f})")

        if len(values) > 1:
            sns.kdeplot(values,
                        color=colors_dbs[i],
                        linewidth=2)

            sns.kdeplot(values,
                        color=colors_dbs[i],
                        fill=True,
                        alpha=0.15,
                        linewidth=0,
                        zorder=2)

        plt.axvline(mean_gscore,
                    color=colors_dbs[i],
                    linewidth=1)

    if ref_gscores:
        for i, gscore in enumerate(map(float, ref_gscores)):
            plt.axvline(gscore,
                        color=colors_refs[i],
                        ls='--',
                        linewidth=1.1,
                        alpha=0.9,
                        zorder=5,
                        label=labels_refs[i])

    if dfs_refs:
        for i, df_ref in enumerate(dfs_refs):
            df_ref = df_ref[df_ref['glide_gscore'] < 50]
            for gscore in df_ref['glide_gscore'].to_numpy():
                plt.axvline(gscore,
                            color=colors_refs[i],
                            ls='--',
                            linewidth=1.1,
                            alpha=0.7,
                            zorder=2,
                            label=labels_refs[i])

    plt.xlabel("Glide GScore (kcal/mol)")
    plt.ylabel("Density" if frequency else "Count")
    plt.title("Distribution of Glide GScore")
                                                                                         
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    plt.ylim(0, 5000)

    sns.despine()

    plt.tight_layout()
    plt.savefig(output + '.pdf', dpi=2200)

    if xlim:
        plt.xlim(xlim[0], xlim[1])
    if ylim:
        plt.ylim(ylim[0], ylim[1])
    if xlim or ylim:
        plt.tight_layout()
        plt.savefig(output + '_zoom.pdf', dpi=2200)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='This script takes various csv\
    files and outputs their histogram plots. Furthermore, you can also\
    modify the plot parameters (colors, labels, etc.)')
    requiredArguments = parser.add_argument_group('required arguments')
    requiredArguments.add_argument('--db_csvs', nargs='+',
                                   help='List of glide simplified \
                                   output csv for DBs compounds',
                                   required=True)
    parser.add_argument('--ref_csvs', nargs='+',
                        help='List of glide simplified output csv\
                        for reference compounds', default=None)
    parser.add_argument('--ref_colors', nargs='+', help='List of\
                        colors for reference compounds vlines',
                        default=None)
    parser.add_argument('--ref_gscores', nargs='+',
                        help='List of gscore values to print as vertical lines.',)
    parser.add_argument('--db_colors', nargs='+', help='List of colors\
                        for DB compound hists')
    parser.add_argument('--db_labels', nargs='+', help='List of labels\
                        for DB compounds hist', default=None)
    parser.add_argument('--ref_labels', nargs='+', help='List of labels\
                        for ref compounds hist', default=None)
    parser.add_argument('-t', help='Transparency value of the hists',
                        default=0.6)
    parser.add_argument('-b', help='Bins value for hists', default=100)
    parser.add_argument('--xlim', nargs='+',
                        help='Hist xlims [left_value, rigth_value]',
                        default=None)
    parser.add_argument('--ylim', nargs='+',
                        help='Hist ylims [left_value, right_value]',
                        default=None)
    parser.add_argument('-o', help='Output file name',
                        default='glide_histogram')
    parser.add_argument('-f', help='Frequency instead of counts in hists',
                        default=False, action='store_true')
    args = parser.parse_args()

    """
    # Check args
    if args.ref_colors:
        if len(args.ref_colors) != len(args.ref_csvs) or len(args.ref_gscores):
            raise ValueError('If provided, ref colors should be of the same\
            length as ref csvs')
    if args.db_colors:
        if len(args.db_colors) != len(args.db_csvs):
            raise ValueError('If provided, db colors should be of the same\
            length as db csvs')
    if args.ref_labels:
        if len(args.ref_labels) != len(args.ref_csvs) or len(args.ref_gscores):
            raise ValueError('If provided, ref labels should be of the same\
            length as ref csvs')
    if args.db_labels:
        if len(args.db_labels) != len(args.db_csvs):
            raise ValueError('If provided, db labels should be of the same\
            length as db csvs')
    if args.xlim:
        if len(args.xlim) != 2:
            raise ValueError('xlim must have left and right value')
    #    else:
    #        if args.xlim[0] > args.xlim[1]:
    #            raise ValueError('left xlim value should be smaller than right\
    #                             xlim value')

    if args.ylim:
        if len(args.ylim) != 2:
            raise ValueError('ylim must have left and right value')
    #    else:
    #        if args.ylim[0] > args.ylim[1]:
    #            raise ValueError('left ylim value should be smaller than right\
    #                             ylim value')
    """
    transparency = float(args.t)
    bins = int(args.b)
    ref_gscores = args.ref_gscores
    output = args.o
    frequency = args.f

    # Load csv as dataframes
    dfs_dbs = []
    for csv in args.db_csvs:
        dfs_dbs.append(pd.read_csv(csv))

    if args.ref_csvs:
        dfs_refs = []
        for csv in args.ref_csvs:
            dfs_refs.append(pd.read_csv(csv))
        # Get colors
        if not args.ref_colors:
            colors_refs = generate_colors_from_palette('tab20', len(dfs_refs))
        else:
            colors_refs = args.ref_colors
        # Get labels
        if not args.ref_labels:
            labels_refs = [None]*len(dfs_refs)
        else:
            labels_refs = args.ref_labels
    else:
        dfs_refs = None
        labels_refs = None
        colors_refs = None

    if ref_gscores:
        if not args.ref_colors:
            colors_refs = generate_colors_from_palette('tab20', len(ref_gscores))
        else:
            colors_refs = args.ref_colors
        # Get labels
        if not args.ref_labels:
            labels_refs = [None]*len(ref_gscores)
        else:
            labels_refs = args.ref_labels
    else:
        ref_gscores = None
        labels_refs_gscores = None
        colors_refs_gscores = None

    if not args.db_colors:
        colors_dbs = generate_colors_from_palette('tab20c', len(dfs_dbs))
    else:
        colors_dbs = args.db_colors

    if not args.db_labels:
        labels_dbs = [None]*len(dfs_dbs)
    else:
        labels_dbs = args.db_labels

    # Get xlim
    if not args.xlim:
        xlim = args.xlim
    else:
        xlim = [float(args.xlim[0]), float(args.xlim[1])]
    if not args.ylim:
        ylim = args.ylim
    else:
        ylim = [int(args.ylim[0]), int(args.ylim[1])]

    plot_hists(dfs_dbs, dfs_refs, colors_dbs, colors_refs, labels_refs,
               ref_gscores, labels_dbs, transparency, bins, xlim, ylim, frequency, output)
