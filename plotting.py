# ==============================================================================
__author__ = "Nicholas Nuechterlein"
__license__ = "MIT"
__maintainer__ = "Nicholas Nuechterlein"
# ==============================================================================
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def scna_manhattan_plot(scna_df,
                        gene_loc_df,scna_df2 = None,
                        take_difference=False,
                        chroms=['6'],
                        bp_start_col='start',
                        chr_col='chr',
                        exclude_genes=[],
                        tick_size=12,
                        label_size=18,
                        ylim=1.05,
                        title='Manattan Plot',
                        save=False,
                        outfile='temp.png',
                        xlabel='\nChromosome',
                        chrom_label_ypos=1.12,
                        chrom_text_dict={},
                        chrom_label_size=14,
                        label_chroms=True,
                        norm=False,
                        sort=True,
                        title_size=30,
                        figsize=(15,8),
                        gain_color='r',
                        loss_color='b',
                        gain_hatch=None,
                        loss_hatch=None,
                        xlabel_fontsize=20,
                        font='Arial'):

    '''
    Plots manhattan plot
    '''
    ## font parameters
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = font
    # reduce the scna data to just genes for which we have location data and in the chromosome we want
    # amd group -2, -1 and 1, 2 together
    common_genes = [x for x in scna_df.columns if x in gene_loc_df.index]
    gene_loc_df = gene_loc_df.loc[common_genes]
    gene_loc_df = gene_loc_df.loc[gene_loc_df[chr_col].isin(chroms)].sort_values(['numeric_chr', bp_start_col])
    genes = gene_loc_df.index
    scna_df = scna_df[genes].replace({-2:-1, 2:1})

    # drop chosen genes (this isn't necessary anymore (I don't have random gains/losses of misordered genes anymore))
    exclude_genes = [x for x in exclude_genes if x in scna_df.columns]
    scna_df = scna_df.drop(columns=exclude_genes)

    # get starting postion for each gene
    start_dict, end_dict = {str(int(chroms[0])-1):0}, {str(int(chroms[0])-1):0}
    for chrom in chroms:
        start_dict[chrom] = np.min(gene_loc_df[gene_loc_df[chr_col] == str(chrom)][bp_start_col]) + end_dict[str(int(chrom)-1)]
        end_dict[chrom] = np.max(gene_loc_df[gene_loc_df[chr_col] == str(chrom)]['end']) + end_dict[str(int(chrom)-1)]

    # get gene --> x-axis coordinates
    gene_loc_df['chrom_start'] = gene_loc_df[chr_col].copy()
    gene_loc_df['chrom_start'] = gene_loc_df['chrom_start'].replace(start_dict)
    gene_loc_df['global_start'] = gene_loc_df[bp_start_col] + gene_loc_df['chrom_start']

    if take_difference:
        scna_df2 = scna_df2[genes].drop(columns=exclude_genes)
        scna_df2 = scna_df2.replace({-2:-1, 2:1})

        scna_df2_avg = scna_df2.sum(axis=0)/scna_df2.shape[0]
        scna_df_avg = scna_df.sum(axis=0)/scna_df.shape[0]
        scna_diff_df =  scna_df2_avg - scna_df_avg

        ys_pos = scna_diff_df.copy()
        ys_pos[ys_pos < 0] = 0

        ys_neg = scna_diff_df.copy()
        ys_neg[ys_neg > 0] = 0

        plot_df = pd.concat([gene_loc_df['global_start'], ys_neg, ys_pos], axis=1).rename(columns={0:'losses', 1:'gains'}).dropna()


    else:
        #### plot difference in LOSSES ###
        # get losses in df1
        _scna_df = scna_df.copy()
        _scna_df[_scna_df > 0] = 0
        negitives = _scna_df.sum(axis=0)/_scna_df.shape[0]

        #### plot difference in GAINS ###
        # get gains in df1
        _scna_df = scna_df.copy()
        _scna_df[_scna_df < 0] = 0
        positives = _scna_df.sum(axis=0)/_scna_df.shape[0]

        plot_df = pd.concat([gene_loc_df['global_start'], negitives, positives], axis=1).rename(columns={0:'losses', 1:'gains'}).dropna()

    ## figure
    fig, ax = plt.subplots(figsize=figsize)
    plt.title(title, fontsize=title_size)
    plot_df = plot_df.sort_values('global_start')
    plt.plot(plot_df['losses'].values, loss_color)
    plt.fill_between(plot_df['global_start'].values, 0, plot_df['losses'].values, facecolor=loss_color, hatch=loss_hatch)
    plt.plot(plot_df['gains'].values, gain_color)
    plt.fill_between(plot_df['global_start'].values, 0, plot_df['gains'].values, facecolor=gain_color, hatch=gain_hatch)

    for chrom in chroms:
        plt.axvline(x=start_dict[chrom], color='k', linewidth=1, linestyle=':')
        if label_chroms:
            offset = (end_dict[chrom] - start_dict[chrom])/2 - 10000000
            plt.text(x= start_dict[chrom] + offset, y=chrom_label_ypos, s=chrom, color='k', fontsize=chrom_label_size)

    xmin, xmax = np.min(plot_df['global_start'].values), np.max(plot_df['global_start'].values)
    plt.xlim((xmin - xmax*0.0, xmax + xmax*0.0))
    plt.ylim((-ylim, ylim))

    plt.yticks(fontsize=tick_size)
    plt.xticks([])
    plt.xlabel(xlabel, fontsize=xlabel_fontsize)
    plt.ylabel('% with gain or loss', fontsize=20)
    ax.set_yticklabels([str(abs(x)) for x in ax.get_yticks()])
    plt.tight_layout()

    if save:
        print('saving to ', outfile, '. . . ')
        plt.savefig(outfile, transparent=True, dpi=300)
    plt.show()
    return plot_df
