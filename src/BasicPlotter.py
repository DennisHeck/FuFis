import pandas as pd
from matplotlib import pyplot as plt
import sklearn.metrics as metrics
from matplotlib import cm
import numpy as np
import math
import seaborn as sns
import ColoursAndShapes
import matplotlib_venn
from matplotlib.colors import to_hex
from itertools import chain
import upsetplot
from collections import Counter
from pandas.api.types import is_string_dtype
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches


def basic_bars(plot_df, plot_cols, hue_col=None, title=None, output_path=''):
    """Plots a basic barplot with plot_cols=[x,y], allows to select hue levels."""

    f, ax = plt.subplots(figsize=(8, 6))
    ax.set_axisbelow(True)
    ax.grid(True, axis='y', color='#f2f2f2', linewidth=1, which='major')
    sns.barplot(data=plot_df, x=plot_cols[0], y=plot_cols[1], hue=hue_col, ax=ax, alpha=1, edgecolor='k', linewidth=1,
                color='#2d63ad', palette='tab10' if hue_col else None)
    ax.tick_params(axis='both', labelsize=14)
    ax.set_ylabel(plot_cols[1], fontsize=16)
    ax.set_xlabel(plot_cols[0], fontsize=16)
    plt.title(title, fontsize=18, fontweight='bold')
    f.savefig((output_path + '_'.join(plot_cols) + '_Bars.pdf').replace(' ', ''), bbox_inches='tight')


def basic_hist(plot_df, x_col, hue_col=None, hue_order=None, bin_num=None, title=None, output_path='', stat='count',
               cumulative=False):
    """Plots a basic layered histogram which allows for hue, whose order can be defined as well.
    If x_col is not a column in the df, it will be assumed that hue_col names all the columns which are supposed to be
    plotted."""
    if x_col not in plot_df.columns:  # Reformat so that seaborn can interpret x_col as hue.
        plot_df = pd.DataFrame(list(chain(*[[[c, x] for x in plot_df[c].to_numpy()] for c in hue_order])),
                               columns=[hue_col, x_col])
    f, ax = plt.subplots(figsize=(12, 8))
    ax.set_axisbelow(True)
    ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='major')
    hist = sns.histplot(data=plot_df, x=x_col, hue=hue_col, ax=ax, alpha=0.3, bins=bin_num if bin_num else 'auto',
                        linewidth=1, color='#2d63ad', palette='tab10' if hue_col else None, hue_order=hue_order,
                        multiple='layer', fill=not cumulative, element='step', stat=stat, common_norm=False,
                        cumulative=cumulative)
    for patch in ax.patches:
        clr = patch.get_facecolor()
        patch.set_edgecolor(clr)
    ax.tick_params(axis='both', labelsize=14)
    ax.set_ylabel('Count', fontsize=16)
    ax.set_xlabel(str(x_col), fontsize=16)
    if hue_col:
        plt.setp(hist.get_legend().get_texts(), fontsize=14)
        plt.setp(hist.get_legend().get_title(), fontsize=16)
    plt.title(title, fontsize=18, fontweight='bold')
    f.savefig((output_path + str(x_col) + '_' + str(hue_col) + '_Hist.pdf').replace(' ', ''), bbox_inches='tight')


def basic_violin(plot_df, y_col, x_col, x_order=None, hue_col=None, hue_order=None, title=None, output_path='',
                 numerate=False, ylim=None, palette='tab10', xsize=12, ysize=8, boxplot=False, rotation=None):
    """Plots a basic violin plot which allows for hue, whose order can be defined as well."""
    if palette == 'glasbey':
        palette = ColoursAndShapes.categ_colours
    f, ax = plt.subplots(figsize=(xsize, ysize))
    ax.set_axisbelow(True)
    ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='major')
    if not boxplot:
        vio = sns.violinplot(data=plot_df, y=y_col, x=x_col, order=x_order, hue=hue_col, ax=ax, color='#2d63ad',
                             palette=palette if not hue_col else None, hue_order=hue_order)
    else:
        vio = sns.boxplot(data=plot_df, y=y_col, x=x_col, order=x_order, hue=hue_col, ax=ax, color='#2d63ad',
                          palette=palette if not hue_col else None, hue_order=hue_order)
    ax.tick_params(axis='both', labelsize=18)
    ax.set_ylabel(y_col, fontsize=22)
    ax.set_xlabel(x_col, fontsize=22)
    if ylim:
        ax.set_ylim(ylim)
    if numerate:
        x_counts = Counter(plot_df[x_col].values)
        ax.set_xticklabels([x._text+'\n(#'+str(x_counts[x._text])+')' for x in ax.get_xmajorticklabels()])
    if hue_col:
        plt.setp(vio.get_legend().get_texts(), fontsize=14)
        plt.setp(vio.get_legend().get_title(), fontsize=16)
    if rotation:
        plt.xticks(rotation=rotation, ha='center')
    plt.title(title, fontsize=22, fontweight='bold')
    f.savefig((output_path + x_col + '_' + y_col + '_' + str(hue_col) + '_Violin.pdf').replace(' ', ''), bbox_inches='tight')


def multi_mod_plot(plot_df, score_cols, colour_col=None, marker_col=None, output_path='', diagonal=False, title=None,
                   colour_order=None, marker_order=None, line_plot=False, alpha=0.7, xsize=8, ysize=6,
                   xlim=None, ylim=None, msize=30):
    """Compares two scores. For each entry in plot_df plot one dot with [x,y] based on score_col and
    allows to colour all dots based on colour_col, and if marker_col is selected assigns each class a different
    marker.
    line_plot: 2D list of dots which will be connected to a lineplot."""

    main_list = plot_df[[x for x in score_cols+[colour_col, marker_col] if x is not None]].values.tolist()
    main_idx = {x: i for i, x in enumerate([y for y in score_cols+[colour_col, marker_col] if y is not None])}
    if colour_col == marker_col:
        if colour_order:
            marker_order = colour_order
        elif marker_order:
            colour_order = marker_order

    if colour_col:
        if is_string_dtype(plot_df[colour_col]):
            if len(set(plot_df[colour_col])) > len(ColoursAndShapes.tol_vibrant):
                palette = ColoursAndShapes.categ_colours
            else:
                palette = ColoursAndShapes.tol_vibrant
            if colour_order:
                colour_dict = {c: palette[i] for i, c in
                               enumerate(colour_order)}
            else:
                colour_dict = {c: palette[i] for i, c in
                               enumerate(list(set(plot_df[colour_col].values)))}
        else:
            if plot_df[colour_col].min() < 0 and plot_df[colour_col].max() > 0:
                cmap = cm.get_cmap('coolwarm')
                bound = max([abs(plot_df[colour_col].min()), abs(plot_df[colour_col].max())])
                norm = plt.Normalize(-bound, bound)
            else:
                cmap = cm.get_cmap('viridis')
                norm = plt.Normalize(plot_df[colour_col].min(), plot_df[colour_col].max())

    def give_colour(entry):
        """Takes the entry that is supposed to be plotted, and checks whether we have categorical colours, or numerical
        ones."""
        if colour_col:
            if is_string_dtype(plot_df[colour_col]):
                return colour_dict[entry[main_idx[colour_col]]]
            else:
                return to_hex(cmap(norm(entry[2])))
        else:
            return '#032e99'

    f, ax = plt.subplots(figsize=(xsize, ysize))
    ax.set_axisbelow(True)
    ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='major')
    if marker_col:
        if len(set(plot_df[marker_col].values)) == 2:
            marker_selection = ['o', '*']
        else:
            marker_selection = ColoursAndShapes.marker_shapes * max(1, math.ceil(len(set(plot_df[marker_col])) / len(ColoursAndShapes.marker_shapes)))
        if marker_order:
            marker_dict = {c: marker_selection[i] for i, c in enumerate(marker_order)}
        else:
            marker_dict = {c: marker_selection[i] for i, c in enumerate(list(set(plot_df[marker_col].values)))}

        for dot in main_list:
            plt.scatter(x=dot[0], y=dot[1], c=[give_colour(dot)], s=msize, edgecolors=None, linewidth=0, zorder=12,
                        marker=marker_dict[dot[main_idx[marker_col]]], alpha=alpha)
        legend_list = [Line2D([0], [0], marker=marker_dict[mark], color='black' if not marker_col == colour_col else colour_dict[mark], linestyle='None') for mark in marker_dict]
        source_legend = plt.legend(legend_list, list(marker_dict.keys()), markerscale=1.49,
                                   scatterpoints=1, fontsize=10, title=marker_col,
                                   bbox_to_anchor=None if len(legend_list) < 7 else (1.04, 1))
        source_legend.get_title().set_fontsize(11)
        plt.gca().add_artist(source_legend)
    else:
        plt.scatter(x=[x[0] for x in main_list], y=[x[1] for x in main_list], c=[give_colour(x) for x in main_list],
                    s=msize, edgecolors=None, linewidth=0, zorder=12, alpha=alpha)

    if colour_col and is_string_dtype(plot_df[colour_col]) and colour_col != marker_col:
        legend_list = [mpatches.Patch([0], [0], color=colour_dict[col], linestyle='None') for col in colour_dict]
        colour_legend = plt.legend(legend_list, list(colour_dict.keys()), markerscale=1,
                                   scatterpoints=1, fontsize=10, title=colour_col,
                                   bbox_to_anchor=None if len(legend_list) < 7 else (1.04, 1))
        colour_legend.get_title().set_fontsize(11)
        plt.gca().add_artist(colour_legend)

    ax.set_xlabel(score_cols[0], fontsize=16)
    ax.set_ylabel(score_cols[1], fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    if diagonal:
        ax.axline((0, 0), slope=1, linestyle='dotted', color='grey')
    if colour_col and not is_string_dtype(plot_df[colour_col]):
        cax = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.8)
        cax.set_label(colour_col, size=14)
    if line_plot:
        plt.plot([x[0] for x in line_plot], [y[1] for y in line_plot], linestyle='-', color='grey', zorder=12)
    if title:
        plt.title(title, fontsize=22, fontweight='bold', y=1.01)
    f.savefig((output_path + score_cols[0]+'Vs'+score_cols[1]+'_'+str(colour_col)+'_marker'+str(marker_col)+'.pdf').replace(' ', ''), bbox_inches='tight')
    plt.close()


def venn_from_list(plot_list, label_list, plot_path, blob_colours=ColoursAndShapes.tol_highcontrast, title='',
                   scaled=True, linestyle='', number_size=11):
    """Based on a list with the size of the sets and the respective labels, plot a non-scaled / scaled Venn diagram
     for up to three sets. If sets are given, the intersection will be done automatically.
     Choose non-scaled if the difference is too high.
     two sets: [a-b, b-a, a∩b]
     three sets: [a-b-c, b-a-c, a∩b-c, c-a-b, a∩c-b, b∩c-a, a∩b∩c]"""
    if sum([type(x) == set for x in plot_list]) == len(plot_list):
        if len(plot_list) == 2:
            a, b = plot_list
            plot_list = [len(a - b), len(b - a), len(a & b)]
        if len(plot_list) == 3:
            a, b, c = plot_list
            plot_list = [len(a - b - c), len(b - a - c), len(a & b - c), len(c-a-b), len(a & c - b), len(b & c - a), len(a & b & c)]
    f, ax = plt.subplots(figsize=(5, 5))
    if scaled and len(plot_list) == 3:
        v = matplotlib_venn.venn2(subsets=plot_list, set_labels=label_list, ax=ax, set_colors=blob_colours,
                                  normalize_to=0.8)
    elif not scaled and len(plot_list) == 3:
        v = matplotlib_venn.venn2_unweighted(subsets=plot_list, set_labels=label_list, ax=ax,
                                             set_colors=blob_colours, normalize_to=0.8)
    elif scaled and len(plot_list) == 7:
        v = matplotlib_venn.venn3(subsets=plot_list, set_labels=label_list, ax=ax, set_colors=blob_colours,
                                  normalize_to=0.8)
    elif not scaled and len(plot_list) == 7:
        v = matplotlib_venn.venn3_unweighted(subsets=plot_list, set_labels=label_list, ax=ax, set_colors=blob_colours,
                                  normalize_to=0.8)
    if len(plot_list) == 3:
        matplotlib_venn.venn2_circles(subsets=plot_list, linestyle=linestyle, linewidth=1, color="grey", normalize_to=0.8)
    elif len(plot_list) == 7:
        matplotlib_venn.venn3_circles(subsets=plot_list, linestyle=linestyle, linewidth=1, color="grey", normalize_to=0.8)
    for text in v.set_labels:
        if text:
            text.set_fontsize(14)
    for text in v.subset_labels:
        if text:
            text.set_fontsize(number_size)
    plt.title(title, size=16)
    plt.savefig((plot_path + '_Venn.pdf').replace(' ', ''), bbox_inches='tight')


def upset_plotter(inter_sets, max_groups=None, sort_by='cardinality', y_label='Intersection', title_tag='', plot_path=''):
    """Based on a dictionary with sets as values creates the intersection and an upsetplot.
    max_groups: defines the maximum number of intersections plotted, sorted descending by size
    
    notes for that hideous API:
    intersection_plot_elements: height
    totals_plot_elements: ~size of the horizontal bars for total size
    element_size: ~overall size and margins
    """
    intersection = upsetplot.from_contents(inter_sets)
    fig = plt.figure(figsize=(10 + int(len(inter_sets) / 2), 7 + int(len(inter_sets) / 2)))
    max_string = max([len(k) for k in inter_sets])
    if max_groups:
        upset = upsetplot.UpSet(intersection, show_counts=False, intersection_plot_elements=0, min_degree=0,
                                with_lines=True, totals_plot_elements=max_groups - 2, element_size=60, sort_by=None,
                                sort_categories_by=None)
        cut_off = upset.intersections.sort_values(ascending=False).tolist()[max_groups]
        filtered_sets = set(upset.intersections[upset.intersections > cut_off].index.values)
        base_totals = upset.totals
        filtered_intersection = intersection[intersection.index.isin(filtered_sets)]
        filtered_upset = upsetplot.UpSet(filtered_intersection, show_counts=True, min_degree=0, with_lines=True,
                                         intersection_plot_elements=5+len(inter_sets), facecolor="#010d4a",
                                         show_percentages=False, totals_plot_elements=max(len(base_totals) - 5, 3),
                                         element_size=len(inter_sets)*1.2+max_string+30, sort_by=sort_by,
                                         sort_categories_by=None)
        filtered_upset.totals = base_totals
    else:
        filtered_upset = upsetplot.UpSet(intersection, show_counts=True, min_degree=0, with_lines=True,
                                         intersection_plot_elements=5+len(inter_sets), facecolor="#010d4a",
                                         show_percentages=True, totals_plot_elements=max(len(inter_sets) - 5, 3),
                                         element_size=len(inter_sets)*1.2+max_string+30, sort_by=sort_by,
                                         sort_categories_by=None)

    with plt.rc_context({'axes.titlesize': 20,
                         'axes.labelsize': 14,
                         'xtick.labelsize': 11,  # Ensures enough space between horizontal bars and UpSet.
                         'ytick.labelsize': 12,
                         'font.size': 12,  # For whatever reason font.size is only for the counts on top of the bars.
                         'legend.fontsize': 14}):
        filtered_upset.plot(fig)

    ax = plt.gca()
    ax.set_ylabel(y_label)
    ax.set_title("UpSet " + title_tag + '\n#' + str(len(intersection)), fontsize=16, fontweight='bold', y=1.05)
    plt.savefig((plot_path + '_UpSet.pdf').replace(' ', ''), bbox_inches='tight')

