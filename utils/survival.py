# ==============================================================================
__author__ = "Nicholas Nuechterlein"
__license__ = "MIT"
__maintainer__ = "Nicholas Nuechterlein"
# ==============================================================================
import pandas as pd
import numpy as np
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors as mcolors

from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter
from lifelines.utils import median_survival_times
from lifelines.statistics import multivariate_logrank_test

def survival_curves(kmf_df,
                    save=False,
                    outfile='temp.png',
                    title = 'Survival Curve',
                    show_p_value = True,
                    text_xpos=5,
                    text_ypos = 0.5,
                    dpi=300,
                    hline=False,
                    legend=True,
                    legend_outside=False,
                    markerscale=1,
                    colors=[],
                    linestyles=[],
                    alphas=[],
                    linewidths=[],
                    figsize=(15, 8),
                    ci_show=False,
                    exclude_labels=[],
                    style_dict={},
                    pval_dec=4,
                    legend_size=20,
                    show=True,
                    verbose=True,
                    p_val_labels=[],
                    p_longtail=0.2,
                    show_longtail=True,
                    xlim=None,
                    ylim=None,
                    y_label='Survival Time (Years)',
                    x_label='Fraction Surviving',
                    xlabel_size = 24,
                    ylabel_size = 24,
                    xtick_size=18,
                    ytick_size=18,
                    titlesize=32,
                    show_hr=False,
                    hr_replace_dict={},
                    p_size=20,
                    show_os=True,
                    draw_median_line=False,
                    transparent=True,
                    legend_loc='upper right',
                    legend_order=None,
                    text_dict=None,
                    plot_lines_dict=None,
                    font='Arial',
                    show_censors=True):

    ## font parameters
    mpl.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = font

    T = kmf_df[kmf_df.columns[0]]/365.25
    E = kmf_df[kmf_df.columns[1]]
    label_df = kmf_df[kmf_df.columns[2]]
    ax = plt.subplot(111)
    y = label_df.values.flatten()
    survival_stats = {}

    vals, cnts = np.unique(y, return_counts=True)
    labels, counts = list(vals), list(cnts)

    if len(labels) > len(linestyles):
        linestyles = linestyles + ['-']*len(labels)
    if len(labels) > len(colors):
        colors = colors + list(mcolors.TABLEAU_COLORS.values()) + list(mcolors.XKCD_COLORS.values())[:len(labels)]
    if len(labels) > len(alphas):
        alphas = alphas + [0.8]*len(labels)
    if len(labels) > len(linewidths):
        linewidths = linewidths + [2]*len(labels)
    if len(style_dict) == 0:
        for i, label in enumerate(labels):
            if label not in exclude_labels:
                style_dict[label] = {'color':colors[i],  'linestyle':linestyles[i], 'alpha':alphas[i], 'linewidth':linewidths[i]}
    survival_stats['labels'] = {}
    for i, label in enumerate(labels):
        if label not in exclude_labels:

            ix = (label_df == label)
            kmf_temp = KaplanMeierFitter()
            kmf_temp.fit(T[ix], E[ix], label=label)
            median_survival = np.round(kmf_temp.median_survival_time_, 1)
            if draw_median_line:
                plt.axvline(x=kmf_temp.median_survival_time_, ymin=0, ymax=0.50, color='k', linestyle='--')
            longtail_survival = np.round(kmf_temp.percentile(p_longtail), 1)
            survival_stats['labels'][str(label)] = {'median_survival':np.round(kmf_temp.median_survival_time_, 4), 'tail_'+str(p_longtail)+'_percentile':np.round(kmf_temp.percentile(p_longtail), 4)}

            kmf = KaplanMeierFitter()
            if show_os and show_longtail:
                legend_label = str(label) + ' (' + str(median_survival) +' years, '+str(int(p_longtail*100))+'% live '+str(longtail_survival) + ' years)'
            elif not show_os and show_longtail:
                legend_label = str(label) + ' ('+str(int(p_longtail*100))+'% live '+str(longtail_survival) + ' years)'
            elif show_os and not show_longtail:
                legend_label = str(label) + ' (' + str(median_survival) +' years)'
            else:
                legend_label = str(label)
            ax = kmf.fit(T[ix], E[ix], label=label).plot(ax=ax,
                                                      figsize=figsize,
                                                      fontsize=18,
                                                      label=legend_label,
                                                      ci_show=ci_show,
                                                      color=style_dict[label]['color'],
                                                      linestyle=style_dict[label]['linestyle'],
                                                      alpha=style_dict[label]['alpha'],
                                                      linewidth=style_dict[label]['linewidth'],
                                                      show_censors=show_censors)




    p_val_idxs = label_df[~label_df.isin(exclude_labels)].index # calcuate p-value only on curves that are displayed in plot (important!)

    kmf_return_df = pd.concat([T, label_df, E], axis=1, join='inner').loc[p_val_idxs]
    # survival_stats['kmf'] = kmf_return_df
    if show_hr:
        os_col_name, label_col_name, event_col_name = kmf_return_df.columns[0], kmf_return_df.columns[1],kmf_return_df.columns[2]
        # label_col_name = 'cluster'
        hr_kmf_df = copy.deepcopy(kmf_return_df)
        if not hr_replace_dict:
            for i, val in enumerate(np.unique(hr_kmf_df[label_col_name])):
                hr_replace_dict[val] = i

        hr_kmf_df[label_col_name] = hr_kmf_df[label_col_name].replace(hr_replace_dict)
        cph = CoxPHFitter()
        cph.fit(hr_kmf_df, duration_col=os_col_name, event_col=event_col_name)
        # HR = cph.hazard_ratios_[label_col_name]
        HR = cph.summary['coef'][label_col_name]
        survival_stats['HR'] = HR

    results = multivariate_logrank_test(T.loc[p_val_idxs],  label_df.loc[p_val_idxs], E.loc[p_val_idxs])
    p_value = results.p_value
    if p_value >= 1e-3:
        p_value = '= ' + str(np.round(p_value, pval_dec))
    elif p_value < 1e-3 and p_value >= 1e-4:
        p_value = '< 1e-3'
    elif p_value < 1e-4 and p_value >= 1e-5:
        p_value = '< 1e-4'
    elif p_value < 1e-5 and p_value >= 1e-6:
        p_value = '< 1e-5'
    elif p_value < 1e-6 and p_value >= 1e-7:
        p_value = '< 1e-6'
    elif p_value < 1e-7 and p_value >= 1e-8:
        p_value = '< 1e-7'
    elif p_value < 1e-8 and p_value >= 1e-9:
        p_value = '< 1e-8'
    elif p_value < 1e-9 and p_value >= 1e-10:
        p_value = '< 1e-9'
    elif p_value < 1e-10 and p_value >= 1e-11:
        p_value = '< 1e-10'
    elif p_value < 1e-11:
        p_value = '<< 1e-10'
    else:
        p_value = '= '+str(p_value)
    survival_stats['p-value'] = p_value
    if show_p_value and show_hr:
        plt.text(x = text_xpos, y=text_ypos, s=r"HR = "+str(np.round(HR, 2)) + '; p ' + str(p_value), color='black', fontsize=p_size)
    elif show_p_value:
        plt.text(x = text_xpos, y=text_ypos, s=r"p "+str(p_value), color='black', fontsize=p_size)
    elif show_hr:
        plt.text(x = text_xpos, y=text_ypos, s='HR='+str(np.round(HR, 2)), color='black', fontsize=p_size)

    if text_dict:
        plt.text(x = text_dict['x'], y=text_dict['y'], s=text_dict['s'], color=text_dict['color'], fontsize=text_dict['fontsize'])
    if hline:
        plt.axhline(y=0.5, color='k', linestyle=':')
    if plot_lines_dict:
        for key in plot_lines_dict:
            plot_line_dict = plot_lines_dict[key]
            plt.plot(plot_line_dict['xlist'], plot_line_dict['ylist'], color=plot_line_dict['color'], linestyle=plot_line_dict['linestyle'], linewidth=plot_line_dict['linewidth'])
    if legend:
        handles, labels = plt.gca().get_legend_handles_labels()
        if legend_outside:
            if legend_order:
                legend = plt.legend([handles[idx] for idx in legend_order],[labels[idx] for idx in legend_order],
                                    fontsize=legend_size, markerscale=markerscale, bbox_to_anchor=(1.04,1), loc="upper left")
            else:
                plt.legend(fontsize=legend_size, markerscale=markerscale, bbox_to_anchor=(1.04,1), loc="upper left")
        else:
            if legend_order:
                legend = plt.legend([handles[idx] for idx in legend_order],[labels[idx] for idx in legend_order],
                                    fontsize=legend_size, markerscale=markerscale, loc=legend_loc)
            else:
                plt.legend(fontsize=legend_size, markerscale=markerscale, loc=legend_loc)
    else:
        ax.legend().set_visible(False)

    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.title(title, fontsize=titlesize)
    plt.xlabel(y_label, fontsize=xlabel_size)
    plt.ylabel(x_label, fontsize=ylabel_size)
    plt.xticks(fontsize=xtick_size)
    plt.yticks(fontsize=ytick_size)
    plt.tight_layout()
    if save:
        plt.savefig(outfile, transparent=transparent, dpi=dpi)
    if show:
        plt.show()
    survival_stats['kmf'] = kmf_return_df
    survival_stats['ax'] = ax
    return survival_stats
