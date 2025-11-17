import copy
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import numpy as np
import math
from matplotlib.patches import Patch
import matplotlib_venn
from matplotlib.colors import to_hex
from itertools import chain
import upsetplot
from collections import Counter
from pandas.api.types import is_string_dtype
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from adjustText import adjust_text
import seaborn as sns
import scipy.stats
import itertools
import ColoursAndShapes


def basic_bars(plot_df, x_col, y_col, x_order=None, hue_col=None, hue_order=None, title=None, output_path='', y_label='',
               x_size=8, y_size=6, rotation=None, palette=None, legend=True, font_s=14, legend_out=False, xlim=None, ylim=None,
               hlines=[], vlines=[], flip_y=False, flip_pad=-30, colour='#2d63ad', edgecolour='black', formats=['pdf']):
    """
    Plots a basic barplot, allows to select hue levels.

    Args:
        y_col: Can be list, then the df will be transformed long format and var_name set to hue_col.
        y_label: to have an appropriate y-axis label.
        flip_y: Very specific case, where the y-ticklabels are flipped into the plot.
        flip_pad: If flip_y, the padding by which the y-ticklabels are moved.
        colour: Colour of the bars, only used if no hue_cols is given.
    """
    plot_df = copy.deepcopy(plot_df)
    if x_col not in plot_df.columns:  # Assumes the x_col is the index if the column doesn't exist.
        plot_df[x_col] = plot_df.index
    if palette and 'glasbey' in palette:
        palette = ColoursAndShapes.glasbey_palettes[palette]
    if type(y_col) == list:  # Assume that the columns should be stacked and hued.
        plot_df = pd.melt(plot_df, id_vars=x_col, value_vars=y_col, value_name=y_label, var_name=hue_col,
                          ignore_index=False)
        y_col = y_label

    f, ax = plt.subplots(figsize=(x_size, y_size))
    ax.set_axisbelow(True)
    ax.grid(True, axis='y', color='#f2f2f2', linewidth=1, which='major')
    sns.barplot(data=plot_df, x=x_col, y=y_col, order=x_order, hue=hue_col, hue_order=hue_order, ax=ax, alpha=1,
                errorbar=None, edgecolor=edgecolour, linewidth=1, color=colour,
                palette='tab10' if hue_col and not palette else palette)
    if ylim:
        ax.set_ylim(ylim)
    if xlim:
        ax.set_xlim(xlim)
    ax.tick_params(axis='both', labelsize=font_s)
    ax.set_ylabel(y_col if not y_label else y_label, fontsize=font_s+2)
    ax.set_xlabel(x_col, fontsize=font_s+2)
    if rotation:
        ax.tick_params(axis='x', rotation=rotation)
    if hue_col:
        if legend_out:
            ax.legend(prop={'size': 14, 'weight': 'bold'}, loc='upper right',
                      bbox_to_anchor=(2 if type(legend_out) == bool else legend_out, 1))
        else:
            ax.legend(prop={'size': 14, 'weight': 'bold'})
    for pos in hlines:
        plt.axhline(pos, color="#a7a8a7", linestyle="--")
    for pos in vlines:
        plt.axvline(pos, color="#a7a8a7", linestyle="--")
    if not legend and ax.get_legend():
        ax.get_legend().remove()
    if flip_y:
        ax.tick_params(axis='y', direction='in', pad=flip_pad, length=0)
        for label in ax.get_yticklabels():
            label.set_horizontalalignment('left')
        ax.set_axisbelow(False)
    plt.title(title, fontsize=font_s+4, fontweight='bold')

    if type(formats) != list:
        formats = [formats]
    for form in formats:
        f.savefig((output_path + '_'.join([str(x_col), str(y_col)]) + '_Bars.'+form).replace(' ', ''),
                  bbox_inches='tight', format=form)
    plt.close()


def jitter_bars(plot_df, x_col, y_col, hue_col=None, hue_order=None, title=None, output_path='', x_size=8, y_size=6,
               rotation=None, palette=None, font_s=14, ylim=None, numerate=True, numerate_break=True, x_order=None,
                jitter_colour='black', jitter_hue_col=None, legend_out=None, formats=['pdf']):
    """
    Plot the average across the x-column groups as bar plot and overlay the individual entries as jitter.
    """
    if palette and 'glasbey' in palette:
        palette = ColoursAndShapes.glasbey_palettes[palette]
    # TODO handle hue col for the means
    mean_vals = {x: np.mean(plot_df[plot_df[x_col] == x][y_col].values) for x in set(plot_df[x_col].values)}
    plot_df['mean x'] = [mean_vals[x] for x in plot_df[x_col].values]
    f, ax = plt.subplots(figsize=(x_size, y_size))
    ax.set_axisbelow(True)
    ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='major')
    bars = sns.barplot(data=plot_df, x=x_col, y='mean x', hue=hue_col, hue_order=hue_order, ax=ax, alpha=1, edgecolor='k',
                linewidth=1, color='#2d63ad', palette='tab10' if hue_col and not palette else palette, order=x_order)
    sns.stripplot(data=plot_df, x=x_col, y=y_col, jitter=True, ax=ax, hue=jitter_hue_col, order=x_order,
                  color=jitter_colour, dodge=True)
    ax.tick_params(axis='both', labelsize=font_s+4)
    ax.set_ylabel(y_col, fontsize=font_s+8)
    ax.set_xlabel(x_col, fontsize=font_s+8)
    if ylim:
        ax.set_ylim(ylim)
    if numerate:
        if not x_col:
            ax.set_xticklabels(['(#' + str((~plot_df[y_col].isna()).sum()) + ')' for x in ax.get_xmajorticklabels()])
        else:
            count_df = plot_df[[x_col, y_col]][~plot_df[y_col].isna()]
            x_counts = Counter(count_df[x_col].values)
            ax.set_xticklabels([x._text+'\n'*numerate_break+'(#'+str(x_counts[x._text])+')' for x in ax.get_xmajorticklabels()])
    if hue_col:
        plt.setp(bars.get_legend().get_texts(), fontsize=font_s)
        plt.setp(bars.get_legend().get_title(), fontsize=font_s+2)
    if rotation:
        plt.xticks(rotation=rotation, ha='center')
    if legend_out:
        ax.legend(prop={'size': font_s, 'weight': 'bold'}, loc='upper right',
                  bbox_to_anchor=(2 if type(legend_out) == bool else legend_out, 1))
    else:
        ax.legend(prop={'size': font_s, 'weight': 'bold'})
    plt.title(title, fontsize=22, fontweight='bold', y=1.01)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        f.savefig((output_path + str(x_col) + '_' + str(y_col) + '_' + str(hue_col) + '_JitterBars.'+form).replace(' ', ''),
                  bbox_inches='tight', format=form)
    plt.close()


def stacked_bars(plot_df, x_col, y_cols, y_label='', title=None, output_path='', x_size=8, y_size=6,
                 rotation=None, palette=None, legend=True, fraction=False, numerate=False, sort_stacks=True,
                 legend_out=False, width=0.8, vertical=False, hatches=None, font_s=14, formats=['pdf']):
    """
    Plots a stacked barplot, with a stack for each y_col.

    Args:
        x_col: Column name to use for splitting on the x-axis, if not an existent column will take the index.
        y_cols: List of the columns to use for the stacks.
        fraction: If True take all values as fraction of the row sum.
        numerate: Whether to add the total number per x-group to the x-label.
        sort_stacks: Whether to sort the stack groups by size.
        legend_out: False to let seaborn place the legend, otherwise a float by how much it should be shifted in x-direction.
        vertical: Whether to swap the layout of the plot from horizontal to vertical.
        hatches: If given assumes the colour list is meant for the x-axis.
    """
    plot_df = copy.deepcopy(plot_df)
    if vertical:  # To start the index at the top.
        plot_df = plot_df.iloc[::-1]
    count_dict = plot_df[y_cols].sum(axis=1).to_dict()
    if hatches and type(hatches) != list:
        hatches = ColoursAndShapes.hatches[:len(y_cols)]
        print(hatches)
    if fraction:
        y_label = 'Fraction of ' + y_label
        plot_df = plot_df.div(plot_df[y_cols].sum(axis=1).values, axis='rows')
    if sort_stacks:
        sorted_cols = list(plot_df.sum().sort_values(ascending=False).index)
        y_cols = sorted_cols
        plot_df.columns = pd.CategoricalIndex(plot_df.columns.values,
                                              ordered=True,
                                              categories=sorted_cols)
        plot_df = plot_df.sort_index(axis=1)
    if x_col not in plot_df.columns:  # Assumes the x_col is the index if the column doesn't exist.
        plot_df[x_col] = list(plot_df.index)
    if palette and 'glasbey' in palette:
        palette = ColoursAndShapes.glasbey_palettes[palette]
    if type(palette) == str:
        cmap = cm.get_cmap(palette, len(y_cols))
        palette = [cmap(i) for i in range(len(y_cols))]

    f, ax = plt.subplots(figsize=(x_size, y_size))
    ax.set_axisbelow(True)
    ax.grid(True, axis='y', color='#f2f2f2', linewidth=1, which='major')
    plot_df.plot(x=x_col, y=y_cols, kind='barh' if vertical else 'bar', stacked=True, ax=ax, alpha=1,
                 linewidth=2 if hatches else 0, color=palette, width=width, edgecolor='black' if hatches else None)
    if hatches:
        stacks = [thing for thing in ax.containers if isinstance(thing, mpl.container.BarContainer)]
        for s, stack in enumerate(stacks):  # For some reason the iteration is the stack layers.
            for p, patch in enumerate(stack):
                patch.set_hatch(hatches[s])
                patch.set_facecolor(palette[p])
    if legend_out:
        ax.legend(prop={'size': font_s, 'weight': 'bold'}, loc='upper right',
                  bbox_to_anchor=(2 if type(legend_out) == bool else legend_out, 1))
    else:
        ax.legend(prop={'size': font_s, 'weight': 'bold'})

    if hatches:  # TODO might not work when the y_cols are sorted in the function
        legend_elements = [Patch(facecolor='w', edgecolor='black', linewidth=2, label=c_label, hatch=hatches[p])
                           for p, c_label in enumerate(y_cols)]
        if legend_out:
            ax.legend(handles=legend_elements, prop={'size': 30}, loc='upper right',
                      bbox_to_anchor=(2 if type(legend_out) == bool else legend_out, 1))
        else:
            ax.legend(handles=legend_elements, prop={'size': 30})

    if not legend and ax.get_legend():
        ax.get_legend().remove()
    ax.tick_params(axis='both', labelsize=font_s)
    ax.set_ylabel(y_label if not vertical else x_col, fontsize=font_s+2)
    ax.set_xlabel(x_col if not vertical else y_label, fontsize=font_s+2)
    if numerate:
        if not vertical:
            ax.set_xticklabels([x.get_text()+'\n(#' + str(int(count_dict[x.get_text()])) + ')' for x in ax.get_xmajorticklabels()],
                               ha='right' if rotation != 90 and rotation != 0 else 'center')
        else:
            ax.set_yticklabels([y.get_text()+'\n(#' + str(int(count_dict[y.get_text()])) + ')' for y in ax.get_ymajorticklabels()])
    if rotation is not None:
        ax.tick_params(axis='x', rotation=rotation)
    plt.title(title, fontsize=font_s+4, fontweight='bold')
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        f.savefig((output_path + '_'.join([str(x_col), str(y_label)]) + '_'+"Frac"*fraction+'StackedBars.'+form).replace(' ', ''),
                  bbox_inches='tight', format=form)
    plt.close()


def basic_hist(plot_df, x_col, hue_col=None, hue_order=None, bin_num=None, title=None, output_path='', stat='count',
               cumulative=False, palette='tab10', binrange=None, xsize=12, ysize=8, colour='#2d63ad', font_s=14,
               ylabel=None, element='step', alpha=0.3, kde=False, legend_out=False, legend_title=True, fill=True,
               edgecolour=None, multiple='layer', shrink=1, hlines=[], vlines=[], discrete=False, grid=True, 
               linewidth=None, log_scale=False, formats=['pdf']):
    """
    Plots a basic layered histogram which allows for hue, whose order can be defined as well.
    If x_col is not a column in the df, it will be assumed that hue_col names all the columns which are supposed to be
    plotted.

    Args:
        stat: count: show the number of observations in each bin. frequency: show the number of observations divided
            by the bin width. probability or proportion: normalize such that bar heights sum to 1. percent: normalize
            such that bar heights sum to 100. density: normalize such that the total area of the histogram equals 1.
        element: {“bars”, “step”, “poly”}.
        multiple: {“layer”, “dodge”, “stack”, “fill”}
        cumulative: If a cumulative distribution should be plotted instead of a histogram.
        kde: Whether to plot the kernel density estimate.
        discrete: If True, each data point gets their own bar with binwidth=1 and bin_num is ignored.
        legend_out: False to let seaborn place the legend, otherwise a float by how much it should be shifted in x-direction.
        hlines: Plot horizontal dashed grey lines at all positions listed in hlines.
        vlines: Same as hlines but vertical.
        shrink: Float to shrink the size of the bar-width to.
    """

    if x_col not in plot_df.columns:  # Reformat so that seaborn can interpret x_col as hue.
        plot_df = pd.DataFrame(list(chain(*[[[c, x] for x in plot_df[c].to_numpy()] for c in hue_order])),
                               columns=[hue_col, x_col])
    if palette and 'glasbey' in palette:
        palette = ColoursAndShapes.glasbey_palettes[palette][: 1 if not hue_col else len(set(plot_df[hue_col]))]
    f, ax = plt.subplots(figsize=(xsize, ysize))
    ax.set_axisbelow(True)
    if grid:
        ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='major')
    hist = sns.histplot(data=plot_df, x=x_col, hue=hue_col, ax=ax, alpha=alpha, bins=bin_num if bin_num else 'auto',
                        linewidth=linewidth if linewidth is not None else (1 if fill else 3), color=colour, palette=palette if hue_col else None,
                        hue_order=hue_order, multiple=multiple, fill=fill, element=element, stat=stat,
                        discrete=discrete, log_scale=[False, log_scale],
                        common_norm=False, cumulative=cumulative, binrange=binrange, kde=kde, shrink=shrink)
    for patch in ax.patches:
        if not edgecolour:
            clr = patch.get_facecolor()
            patch.set_edgecolor(clr)
        else:
            patch.set_edgecolor(edgecolour)
    ax.tick_params(axis='both', labelsize=font_s)
    ax.set_ylabel("log "*log_scale + stat, fontsize=font_s+2)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=font_s+2)
    ax.set_xlabel(str(x_col), fontsize=font_s+2)
    if hue_col:
        plt.setp(hist.get_legend().get_texts(), fontsize=font_s)
        plt.setp(hist.get_legend().get_title(), fontsize=font_s+2)
        if not legend_title:
            hist.get_legend().set_title('')
        if legend_out:
            sns.move_legend(hist, prop={'size': font_s, 'weight': 'normal'}, loc='upper right',
                            bbox_to_anchor=(2 if type(legend_out) == bool else legend_out, 1))
        else:
            sns.move_legend(hist, prop={'size': font_s, 'weight': 'normal'}, loc='best')
    for pos in hlines:
        plt.axhline(pos, color="#a7a8a7", linestyle="--")
    for pos in vlines:
        plt.axvline(pos, color="#a7a8a7", linestyle="--")
    plt.title(title, fontsize=font_s+4, fontweight='bold')
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        f.savefig((output_path + log_scale*'log' + str(x_col) + '_' + str(hue_col) + '_Hist.'+form).replace(' ', ''),
                  bbox_inches='tight', format=form)
    plt.close()


def basic_2Dhist(plot_df, columns, hue_col=None, hue_order=None, bin_num=200, title=None, output_path='',
                 xsize=12, ysize=8, palette='tab10', cbar=False, cmap='mako', hlines=[], vlines=[], binrange=None,
                 diagonal=False, grid=True, font_s=14, formats=['pdf']):
    """
    Plots a basic 2D histogram as heatmap which allows for hue, whose order can be defined as well. Useful when
    a scatterplot would be just a big blob of dots.

    Args:
        columns: List with 2 entries representing the columns from plot_df for the x- and y-axis.
        bin_num: Number of bins to use, the higher the finer the resolution. It behaves it a odd though.
        cbar: Whether to plot the colorbar that shows the number of counts.
        hlines: Plot horizontal dashed grey lines at all positions listed in hlines.
        vlines: Same as hlines but vertical.
        binrange: Boundaries for the histogram.
        diagonal: Whether to add a line on the diagonal.
    """
    if len(columns) != 2:
        print("ERROR: the 2Dhist can only work with 2 columns")
        return
    if palette and 'glasbey' in palette:
        palette = ColoursAndShapes.glasbey_palettes[palette][:len(set(plot_df[hue_col]))]
    f, ax = plt.subplots(figsize=(xsize, ysize))
    ax.set_axisbelow(True)
    if grid:
        ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='major')
    hist = sns.histplot(data=plot_df, x=columns[0], y=columns[1], hue=hue_col, ax=ax, bins=bin_num, binrange=binrange,
                            palette=palette if hue_col else None, hue_order=hue_order, cbar=cbar,
                            cbar_kws=dict(shrink=.75), rasterized=True, cmap=cmap)
    for patch in ax.patches:
        clr = patch.get_facecolor()
        patch.set_edgecolor(clr)
    ax.tick_params(axis='both', labelsize=font_s)
    ax.set_ylabel(columns[1], fontsize=font_s+2)
    ax.set_xlabel(columns[0], fontsize=font_s+2)
    if hue_col:
        plt.setp(hist.get_legend().get_texts(), fontsize=font_s)
        plt.setp(hist.get_legend().get_title(), fontsize=font_s+2)
    for pos in hlines:
        plt.axhline(pos, color="#a7a8a7", linestyle="--")
    for pos in vlines:
        plt.axvline(pos, color="#a7a8a7", linestyle="--")
    if diagonal:
        ax.axline((0, 0), slope=-1 if diagonal == -1 else 1, linestyle='dotted', color='grey')
    plt.title(title, fontsize=18, fontweight='bold')
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        f.savefig((output_path + columns[0] + "_" + columns[1] + '_' + str(hue_col) + '_2DHist.'+form).replace(' ', ''),
                  bbox_inches='tight', format=form)
    plt.close()


def basic_lineplot(plot_df, x_col, y_col, hue_col=None, hue_order=None, title=None, output_path='', markers=False,
                   x_size=8, y_size=6, palette=None, font_s=14, legend_out=False, ylim=None, xlim=None, dashes=True,
                   linewidth=1, hlines=[], vlines=[], grid=True, formats=['pdf']):
    """
    Plots a basic lineplot, allows to select hue levels.

    Args:
        markers: Whether to add markers add each point, can be a list of markers.
        dashes: Whether to distinguish groups by linestyle, can be a list or dictionary of dashes, as defined by pyplot
            (https://matplotlib.org/stable/gallery/lines_bars_and_markers/line_demo_dash_control.html).
    """
    if x_col not in plot_df.columns:  # Assumes the x_col is the index if the column doesn't exist.
        plot_df[x_col] = plot_df.index
    if palette and 'glasbey' in palette:
        palette = ColoursAndShapes.glasbey_palettes[palette][:len(set(plot_df[hue_col]))]

    f, ax = plt.subplots(figsize=(x_size, y_size))
    ax.set_axisbelow(True)
    if grid:
        ax.grid(True, axis='y', color='#f2f2f2', linewidth=1, which='major')
    sns.lineplot(data=plot_df, x=x_col, y=y_col, hue=hue_col, hue_order=hue_order, ax=ax, alpha=1, markers=markers,
                 linewidth=linewidth, palette='tab10' if hue_col and not palette else palette, dashes=dashes, style=hue_col)
    ax.tick_params(axis='both', labelsize=font_s)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    ax.set_xlabel(x_col, fontsize=font_s+2)
    ax.set_ylabel(y_col, fontsize=font_s+2)
    if hue_col:
        if legend_out:
            ax.legend(prop={'size': 14, 'weight': 'bold'}, loc='upper right',
                      bbox_to_anchor=(2 if type(legend_out) == bool else legend_out, 1))
        else:
            ax.legend(prop={'size': 14, 'weight': 'bold'})
    for pos in hlines:
        plt.axhline(pos, color="#a7a8a7", linestyle="--")
    for pos in vlines:
        plt.axvline(pos, color="#a7a8a7", linestyle="--")
    plt.title(title, fontsize=font_s+4, fontweight='bold')
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        f.savefig((output_path + '_'.join([str(x_col), str(y_col)]) + '_Lineplot.'+form).replace(' ', ''),
                  bbox_inches='tight', format=form)
    plt.close()


def basic_violin(plot_df, y_col, x_col, x_order=None, hue_col=None, hue_order=None, title=None, output_path='',
                 numerate=False, ylim=None, palette=None, xsize=12, ysize=8, boxplot=False, boxplot_meanonly=False,
                 rotation=None, numerate_break=True, jitter=False, colour='#2d63ad', font_s=14, saturation=0.75,
                 jitter_colour='black', jitter_size=5, vertical_grid=False, legend_title=True, legend=True, grid=True,
                 formats=['pdf']):
    """
    Plots a basic violin plot which allows for hue, whose order can be defined as well. Optionally plot boxplot, or
    add jitter points for the individual data points.
    Use y_col=None and x_col=None for seaborn to interpret the columns as separate plots on the x-asis.

    Args:
        numerate: Whether to add the total number per x-group to the x-label.
        numerate_break: Whether to add a line break before writing the size of the x-group.
        jitter: Whether to add jitter points for each data point on top.
        jitter_colour: Colour for the jitter points.
        jitter_size: Dot size for the jitter points.
        saturation: Controls the saturation of the colours.
        boxplot: Whether to plot boxplots instead of violinplots.
        boxplot_meanonly: Remove all lines from the boxplot and show just the mean as horizontal line.
    """
    if palette and 'glasbey' in palette:
        palette = ColoursAndShapes.glasbey_palettes[palette]
    f, ax = plt.subplots(figsize=(xsize, ysize))
    ax.set_axisbelow(True)
    if grid:
        ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='major')
    if not boxplot:
        vio = sns.violinplot(data=plot_df, y=y_col, x=x_col, order=x_order, hue=hue_col, ax=ax,
                             color=colour if not palette else None, saturation=saturation,
                             palette='tab10' if hue_col and not palette else palette, hue_order=hue_order)
    else:
        if boxplot_meanonly:
            vio = sns.boxplot(data=plot_df, y=y_col, x=x_col, order=x_order, hue=hue_col, ax=ax, hue_order=hue_order,
                              showfliers=False, showbox=False, showcaps=False, showmeans=True, meanline=True,
                              meanprops={'color': 'k', 'ls': '-', 'lw': 2}, medianprops={'visible': False},
                              whiskerprops={'visible': False}, zorder=1, saturation=saturation,
                              palette='tab10' if hue_col and not jitter_colour else jitter_colour)
        else:
            vio = sns.boxplot(data=plot_df, y=y_col, x=x_col, order=x_order, hue=hue_col, ax=ax, saturation=saturation,
                              color=colour if not palette else None, showfliers=False if jitter else True,
                              palette='tab10' if hue_col and not palette else palette, hue_order=hue_order)
    if jitter:
        sns.stripplot(data=plot_df, x=x_col, y=y_col, jitter=True, ax=ax, hue=hue_col, hue_order=hue_order, zorder=10,
                      order=x_order, palette=jitter_colour, dodge=True, legend=False, edgecolor='black', linewidth=1,
                      size=jitter_size)
    ax.tick_params(axis='both', labelsize=font_s+4)
    ax.set_ylabel(y_col, fontsize=font_s+8)
    ax.set_xlabel(x_col, fontsize=font_s+8)
    if ylim:
        ax.set_ylim(ylim)
    if numerate:
        if not x_col:
            ax.set_xticklabels(['(#' + str((~plot_df[y_col].isna()).sum()) + ')' for x in ax.get_xmajorticklabels()])
        else:
            count_df = plot_df[[x_col, y_col]][~plot_df[y_col].isna()]
            x_counts = Counter(count_df[x_col].values)
            ax.set_xticklabels([x._text+'\n'*numerate_break+'(#'+str(x_counts[x._text])+')' for x in ax.get_xmajorticklabels()])
    if hue_col:
        plt.setp(vio.get_legend().get_texts(), fontsize=font_s)
        plt.setp(vio.get_legend().get_title(), fontsize=font_s+2)
        if not legend_title:
            vio.get_legend().set_title('')
        sns.move_legend(vio, prop={'size': 14, 'weight': 'normal'}, loc='best')
    if rotation:
        plt.xticks(rotation=rotation, ha='center')
    if vertical_grid:  # Fun part is minor ticks are always x5.
        for x in range(len(set(plot_df[x_col]))):
            plt.axvline(x+0.5, color='#f2f2f2', linewidth=1, zorder=0)
    if not legend:
        ax.get_legend().remove()
    plt.title(title, fontsize=22, fontweight='bold', y=1.02)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        f.savefig((output_path + str(x_col) + '_' + str(y_col) + '_' + str(hue_col) + '_Violin.'+form).replace(' ', ''),
                  bbox_inches='tight', format=form)
    plt.close()


def basic_pie(plot_df, title='', palette=None, numerate=True, legend_perc=True, output_path='', legend_title='',
              formats=['pdf']):
    """
    A pie chart.

    Args:
        plot_df: Either a DataFrame with each category as index or a dictionary with {category: count}.
        numerate: Whether to add the summed count across categories to the title.
        legend_perc: Whether the legend should write the percentages or absolute numbers.
    """
    if type(palette) == list:
        clrs = palette
    elif palette is not None and 'glasbey' not in palette:
        clrs = sns.color_palette(palette, n_colors=len(plot_df.index))
    else:
        clrs = ColoursAndShapes.glasbey_palettes[palette]
    if type(plot_df) != pd.core.frame.DataFrame:
        plot_df = pd.DataFrame.from_dict(plot_df, orient='index')
    count_list = plot_df[plot_df.columns[0]].values
    f, ax = plt.subplots()
    wedges, texts, autotexts = ax.pie(count_list, autopct="", shadow=False, startangle=0, colors=clrs,
                                      textprops={'fontsize': 11})
    if legend_perc:
        legend_labels = [str(round(count_list[i] / plot_df.sum().values[0] * 100, 2)) + '% ' + i_label for i, i_label in
                         enumerate(plot_df.index)]
    else:
        legend_labels = plot_df.index
    ax.legend(wedges, legend_labels,
              title=legend_title,
              loc="center left",
              bbox_to_anchor=(1, 0, 0.5, 1))
    ax.set_title(title + ('\n#' + str(plot_df.sum().values[0])) * numerate, fontsize=14, y=1)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        f.savefig(output_path + str(plot_df.columns[0]).replace(' ', '_') + "_PieChart."+form, bbox_inches='tight',
                  format=form)
    plt.close()


def multi_mod_plot(plot_df, score_cols, colour_col=None, marker_col=None, output_path='', diagonal=False, title=None,
                   colour_order=None, marker_order=None, line_plot=False, alpha=0.7, xsize=8, ysize=6, palette=None,
                   xlim=None, ylim=None, msize=30, vlines=[], hlines=[], add_spear=False, na_colour='black', grid=True,
                   label_dots=None, font_s=14, adjust_labels=True, formats=['pdf']):
    """
    Scatterplot that compares two scores. For each entry in plot_df plot one dot with [x,y] based on score_col and
    allows to colour all dots based on colour_col, and if marker_col is selected assigns each class a different
    marker.

    Args:
        score_cols: List of two column names to use for the x- and y-axis.
        line_plot: 2D list of dots which will be connected to a lineplot.
        hlines: Plot horizontal dashed grey lines at all positions listed in hlines.
        vlines: Same as hlines but vertical.
        diagonal: Whether to add a line on the diagonal.
        add_spear: CARE not properly tested. Whether to add the spearman correlation coefficient between the score_cols to the title.
        na_colour: How to colour dots where the colour_col is NA.
        label_dots: A pair of columns [do_label, label_col] with boolean do_label telling which entries should get a
            text label within the plot, and label_col giving the string of the label.
        adjust_labels: Whether to use the adjustText package to try to avoid overlap of text labels.
    """
    main_list = plot_df[[x for x in score_cols+[colour_col, marker_col] if x is not None]].values.tolist()
    main_idx = {x: i for i, x in enumerate([y for y in score_cols+[colour_col, marker_col] if y is not None])}
    if colour_col == marker_col:
        if colour_order:
            marker_order = colour_order
        elif marker_order:
            colour_order = marker_order

    if colour_col:
        if is_string_dtype(plot_df[colour_col]) or plot_df[colour_col].dtype == bool:
            if not palette:
                if len(set(plot_df[colour_col])) > len(ColoursAndShapes.tol_vibrant):
                    palette = ColoursAndShapes.glasbey_palettes['glasbey']
                else:
                    palette = ColoursAndShapes.tol_vibrant
            else:
                if 'glasbey' in palette:
                    palette = ColoursAndShapes.glasbey_palettes[palette]
            if colour_order:
                colour_dict = {c: palette[i] for i, c in
                               enumerate(colour_order)}
            else:
                colour_dict = {c: palette[i] for i, c in
                               enumerate(list(set(plot_df[colour_col].values)))}
        else:
            if plot_df[colour_col].min() < 0 and plot_df[colour_col].max() > 0:
                cmap = cm.get_cmap('coolwarm' if not palette else palette)
                bound = max([abs(plot_df[colour_col].min()), abs(plot_df[colour_col].max())])
                norm = plt.Normalize(-bound, bound)
            else:
                cmap = cm.get_cmap('viridis' if not palette else palette)
                norm = plt.Normalize(plot_df[colour_col].min(), plot_df[colour_col].max())

    def give_colour(entry):
        """Takes the entry that is supposed to be plotted, and checks whether we have categorical colours, or numerical
        ones."""
        if colour_col:
            if is_string_dtype(plot_df[colour_col]) or plot_df[colour_col].dtype == bool:
                return colour_dict[entry[main_idx[colour_col]]]
            else:
                if pd.isna(entry[2]):
                    return na_colour
                return to_hex(cmap(norm(entry[2])))
        else:
            return '#032e99'

    f, ax = plt.subplots(figsize=(xsize, ysize))
    ax.set_axisbelow(True)
    if grid:
        ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='major')
    if marker_col:
        if len(set(plot_df[marker_col].values)) == 2:
            marker_selection = ['o', 'd']
        else:
            marker_selection = ColoursAndShapes.marker_shapes * max(1, math.ceil(len(set(plot_df[marker_col])) / len(ColoursAndShapes.marker_shapes)))
        if marker_order:
            marker_dict = {c: marker_selection[i] for i, c in enumerate(marker_order)}
        else:
            marker_dict = {c: marker_selection[i] for i, c in enumerate(list(set(plot_df[marker_col].values)))}

        # Iterate over the marker groups to not have to plot every dot individually.
        for marker_group in marker_dict:
            marker_entries = [x for x in main_list if x[main_idx[marker_col]] == marker_group]
            plt.scatter(x=[x[0] for x in marker_entries], y=[x[1] for x in marker_entries], c=[give_colour(x) for x in marker_entries], s=msize, edgecolors=None, linewidth=0, zorder=12,
                        marker=marker_dict[marker_group], alpha=alpha)
        legend_list = [Line2D([0], [0], marker=marker_dict[mark], color='black' if not marker_col == colour_col else colour_dict[mark], linestyle='None') for mark in marker_dict]
        source_legend = plt.legend(legend_list, list(marker_dict.keys()), markerscale=1.49,
                                   scatterpoints=1, fontsize=font_s-4, title=marker_col, loc='best',
                                   bbox_to_anchor=None if len(legend_list) < 7 and not colour_col else (1.03, 1))
        source_legend.get_title().set_fontsize(11)
        plt.gca().add_artist(source_legend)
    else:
        plt.scatter(x=[x[0] for x in main_list], y=[x[1] for x in main_list], c=[give_colour(x) for x in main_list],
                    s=msize, edgecolors=None, linewidth=0, zorder=12, alpha=alpha)

    if colour_col and (is_string_dtype(plot_df[colour_col]) or plot_df[colour_col].dtype == bool) and colour_col != marker_col:
        legend_list = [mpatches.Patch([0], [0], color=colour_dict[col], linestyle='None') for col in colour_dict]
        y_offset = 0 if not marker_col else 0.08 * len(marker_dict)
        colour_legend = plt.legend(legend_list, list(colour_dict.keys()), markerscale=1,
                                   scatterpoints=1, fontsize=font_s-4, title=colour_col, loc='best',
                                   bbox_to_anchor=None if len(legend_list) < 7 and not marker_col else (1.03, 1 - y_offset))
        colour_legend.get_title().set_fontsize(font_s-3)
        plt.gca().add_artist(colour_legend)

    if label_dots:
        texts = []
        to_label_df = plot_df[plot_df[label_dots[0]]]
        for i, entry in to_label_df.iterrows():
            texts.append(plt.text(x=entry[score_cols[0]], y=entry[score_cols[1]], s=entry[label_dots[1]],
                                  fontsize=font_s-2, color=give_colour([None, None, entry[colour_col]]), zorder=42))
        if adjust_labels:
            adjust_text(texts, arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    ax.set_xlabel(score_cols[0], fontsize=font_s)
    ax.set_ylabel(score_cols[1], fontsize=font_s)
    ax.tick_params(axis='both', labelsize=font_s-2)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    if diagonal:
        ax.axline((0, 0), slope=1, linestyle='dotted', color='grey')
    for pos in hlines:
        plt.axhline(pos, color="#a7a8a7", linestyle="--")
    for pos in vlines:
        plt.axvline(pos, color="#a7a8a7", linestyle="--")
    if colour_col and not is_string_dtype(plot_df[colour_col]) and plot_df[colour_col].dtype != bool:
        cax = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.8 if not marker_col else 0.6)
        cax.set_label(colour_col, size=font_s)
    if line_plot is not False and line_plot is not None:
        plt.plot([x[0] for x in line_plot], [y[1] for y in line_plot], linestyle='-', color='grey', zorder=12)
    if title:
        plt.suptitle(title, fontsize=font_s+6, fontweight='bold', y=1.01)
    if add_spear:
        spear_r, pval = scipy.stats.spearmanr(a=[x[0] for x in main_list], b=[x[1] for x in main_list])
        plt.title('spear_r='+str(round(spear_r, 3)), fontsize=font_s+4)

    if type(formats) != list:
        formats = [formats]
    for form in formats:
        f.savefig((output_path + str(score_cols[0])+'Vs'+str(score_cols[1])+'_'+str(colour_col)+'_marker'+str(marker_col)+'.'+form).replace(' ', ''),
                  bbox_inches='tight', format=form)
    plt.close()

    if add_spear:
        return spear_r


def basic_venn(input_sets, plot_path, blob_colours=ColoursAndShapes.tol_highcontrast, title='',
               scaled=True, linestyle='', number_size=11, xsize=5, ysize=5, normalize_to=0.5, formats=['pdf']):
    """
    Based on a dictionary of {key: set} with either two or three entries do the intersection of sets and create a
    Venn diagram from it.

    Args:
        input_sets: {key: set} for all bubbles that should be plotted.
        scaled: If the bubble sizes should be scaled to the set sizes. Choose False if the size difference is too high.
        linestyle: Linestyle for the rim of the bubbles.
        normalize_to: Factor to scale the bubble size. If bubbles get cut, try decreasing it.
     """
    if 'glasbey' in blob_colours:
        blob_colours = ColoursAndShapes.glasbey_palettes[blob_colours][:len(input_sets)]
    if len(input_sets) > 3:
        print("ERROR, only making Venns for three or two sets")
        return
    key_order = list(input_sets.keys())
    if len(key_order) == 2:
        a, b = [input_sets[x] for x in key_order]
        plot_list = [len(a - b), len(b - a), len(a & b)]
    elif len(key_order) == 3:
        a, b, c = [input_sets[x] for x in key_order]
        plot_list = [len(a - b - c), len(b - a - c), len(a & b - c), len(c-a-b), len(a & c - b), len(b & c - a), len(a & b & c)]
    f, ax = plt.subplots(figsize=(xsize, ysize))
    if scaled and len(plot_list) == 3:
        v = matplotlib_venn.venn2(subsets=plot_list, set_labels=key_order, ax=ax, set_colors=blob_colours,
                                  normalize_to=normalize_to)
    elif not scaled and len(plot_list) == 3:
        v = matplotlib_venn.venn2_unweighted(subsets=plot_list, set_labels=key_order, ax=ax,
                                             set_colors=blob_colours, normalize_to=normalize_to)
    elif scaled and len(plot_list) == 7:
        v = matplotlib_venn.venn3(subsets=plot_list, set_labels=key_order, ax=ax, set_colors=blob_colours,
                                  normalize_to=normalize_to)
    elif not scaled and len(plot_list) == 7:
        v = matplotlib_venn.venn3_unweighted(subsets=plot_list, set_labels=key_order, ax=ax, set_colors=blob_colours,
                                  normalize_to=normalize_to)
    if len(plot_list) == 3:
        matplotlib_venn.venn2_circles(subsets=plot_list, linestyle=linestyle, linewidth=1, color="grey", normalize_to=normalize_to)
    elif len(plot_list) == 7:
        matplotlib_venn.venn3_circles(subsets=plot_list, linestyle=linestyle, linewidth=1, color="grey", normalize_to=normalize_to)
    for text in v.set_labels:
        if text:
            text.set_fontsize(14)
    for text in v.subset_labels:
        if text:
            text.set_fontsize(number_size)
    plt.title(title, size=16, y=1.15)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        plt.savefig((plot_path + title + '_Venn.'+form).replace(' ', ''), bbox_inches='tight', format=form)
    plt.close('All')


def upset_plotter(inter_sets, y_label='Intersection', max_groups=None, min_degree=0, show_percent=False,
                  sort_categories_by='cardinality', sort_by='cardinality', intersection_plot_elements=None,
                  element_size=None, title_tag='', plot_path='', font_s=14, formats=['pdf']):
    """
    Based on a dictionary with sets as values creates the intersection and an upsetplot. Uses the upsetplot package
    (https://upsetplot.readthedocs.io/en/stable/#) based on doi.org/10.1109/TVCG.2014.2346248. The size of the resulting
    plot isn't that easy to control with the flags, the package isn't very user-friendly in that regard.

    Args:
        inter_sets: Dictionary of {key: set}, each key will be plotted as one row.
        y_label: Label for the y-axis.
        max_groups: Defines the maximum number of intersections plotted, sorted descending by size.
        min_degree: Minimum overlap between groups to show.
        show_percent: Write the percentage of overlap on top of the bars.
        sort_categories_by: How to sort the rows, 'cardinality', 'degree' or 'input'.
        sort_by: How to order the columns, meaning the groups of intersections, 'cardinality' or 'degree'.
        intersection_plot_elements: Roughly the height of the plot.
        element_size: Roughly the overall size and margins.
    """
    # totals_plot_elements: ~size of the horizontal bars for total size
    if sort_categories_by == 'input':  # Reverse the order, to have it from top to bottom.
        inter_sets = {k: inter_sets[k] for k in list(inter_sets.keys())[::-1]}
        sort_categories_by = None
    intersection = upsetplot.from_contents(inter_sets)
    fig = plt.figure(figsize=(10 + int(len(inter_sets) / 2), 7 + int(len(inter_sets) / 2)))
    max_string = max([len(k) for k in inter_sets])
    if max_groups:
        upset = upsetplot.UpSet(intersection, show_counts=False, intersection_plot_elements=0, min_degree=min_degree,
                                with_lines=True, totals_plot_elements=max_groups - 2, element_size=60, sort_by=None,
                                sort_categories_by=sort_categories_by)
        cut_off = upset.intersections.sort_values(ascending=False).tolist()[max_groups]
        filtered_sets = set(upset.intersections[upset.intersections > cut_off].index.values)
        base_totals = upset.totals
        filtered_intersection = intersection[intersection.index.isin(filtered_sets)]
        filtered_upset = upsetplot.UpSet(filtered_intersection, show_counts=True, min_degree=min_degree, with_lines=True,
                                         intersection_plot_elements=5+len(inter_sets) if not intersection_plot_elements else intersection_plot_elements,
                                         facecolor="#010d4a",
                                         show_percentages=show_percent, totals_plot_elements=max(len(base_totals) - 5, 3),
                                         element_size=len(inter_sets)*1.2+max_string+30 if not element_size else element_size, sort_by=sort_by,
                                         sort_categories_by=sort_categories_by)
        filtered_upset.totals = base_totals
    else:
        filtered_upset = upsetplot.UpSet(intersection, show_counts=True, min_degree=min_degree, with_lines=True,
                                         intersection_plot_elements=5+len(inter_sets) if not intersection_plot_elements else intersection_plot_elements,
                                         facecolor="#010d4a",
                                         show_percentages=show_percent, totals_plot_elements=max(len(inter_sets) - 5, 3),
                                         element_size=len(inter_sets)*1.2+max_string+30 if not element_size else element_size, sort_by=sort_by,
                                         sort_categories_by=sort_categories_by)

    with plt.rc_context({'axes.titlesize': font_s+6,
                         'axes.labelsize': 14,
                         'xtick.labelsize': font_s-3,  # Ensures enough space between horizontal bars and UpSet.
                         'ytick.labelsize': font_s-2,
                         'font.size': font_s-3,  # For whatever reason font.size is only for the counts on top of the bars.
                         'legend.fontsize': font_s}):
        filtered_upset.plot(fig)

    ax = plt.gca()
    ax.set_ylabel(y_label)
    ax.set_title("UpSet " + title_tag + '\n#' + str(len(intersection)), fontsize=16, fontweight='bold', y=1.05)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        plt.savefig((plot_path + '_UpSet.'+form).replace(' ', ''), bbox_inches='tight', format=form)
    plt.close()


def overlap_heatmap(inter_sets, title="", plot_path='', xsize=12, ysize=8, annot=True, mode='Jaccard', annot_type='Jaccard',
                    font_s=14, matrix_only=False, cmap='mako', formats=['pdf']):
    """
    Create heatmap of pairwise overlap of a collection of sets. Either make a symmetric heatmap of the Jaccard index
    or plot the asymmetric fraction of overlap. Combinations are formed by the keys of the dictionary.

    Args:
        mode: 'Jaccard' to show the Jaccard index of the intersection, 'Fraction' to show the percentage.
        annot_type: 'Jaccard' to get the Jaccard index written into the cells, 'Fraction' for the fraction, and 'Absolute' to get the absolute number of shared items.
        matrix_only: Whether only the matrix should be given without creating a plot, will return two dfs, one with
            the shared metric and the other the cell labels. If False, only producing the plot without returning matrices.
    """
    if mode == 'Jaccard':
        cbar_label = "Jaccard index"
    elif mode == 'Fraction':
        cbar_label = "Fraction shared [ (row & column) / row ]"
    else:
        print("ERROR: Invalid mode", mode)
        return
    set_tags = list(inter_sets.keys())
    combos = list(itertools.combinations(set_tags, 2))
    shared_mat = np.zeros([len(inter_sets), len(inter_sets)])
    annot_mat = np.zeros([len(inter_sets), len(inter_sets)])

    for comb in combos:
        shared = len(inter_sets[comb[0]] & inter_sets[comb[1]])
        if shared == 0:
            continue
        if mode == 'Jaccard':
            shared_mat[set_tags.index(comb[0])][set_tags.index(comb[1])] = shared / len(inter_sets[comb[0]].union(inter_sets[comb[1]]))
            shared_mat[set_tags.index(comb[1])][set_tags.index(comb[0])] = shared / len(inter_sets[comb[0]].union(inter_sets[comb[1]]))
        elif mode == "Fraction":
            shared_mat[set_tags.index(comb[0])][set_tags.index(comb[1])] = shared / len(inter_sets[comb[0]])
            shared_mat[set_tags.index(comb[1])][set_tags.index(comb[0])] = shared / len(inter_sets[comb[1]])
        if annot_type == "Jaccard":
            annot_mat[set_tags.index(comb[0])][set_tags.index(comb[1])] = shared / len(inter_sets[comb[0]].union(inter_sets[comb[1]]))
            annot_mat[set_tags.index(comb[1])][set_tags.index(comb[0])] = shared / len(inter_sets[comb[0]].union(inter_sets[comb[1]]))
        else:
            annot_mat[set_tags.index(comb[0])][set_tags.index(comb[1])] = len(inter_sets[comb[0]] & inter_sets[comb[1]])
            annot_mat[set_tags.index(comb[1])][set_tags.index(comb[0])] = len(inter_sets[comb[0]] & inter_sets[comb[1]])
    if annot_type == 'Jaccard' or annot_type == 'Fraction':
        annot_mat = annot_mat.round(3)

    # Fill the diagonal.
    np.fill_diagonal(shared_mat, 1)
    if annot_type == 'Absolute':
        np.fill_diagonal(annot_mat, [len(inter_sets[set_tags[n]]) for n in range(len(set_tags))])
        annot_mat = annot_mat.astype(int)
    else:
        np.fill_diagonal(annot_mat, 1)

    if matrix_only:
        shared_df = pd.DataFrame(shared_mat, columns=set_tags, index=set_tags)
        annot_df = pd.DataFrame(annot_mat, columns=set_tags, index=set_tags)
        return shared_df, annot_df

    f, ax = plt.subplots(figsize=(xsize, ysize))
    shared_heat = sns.heatmap(shared_mat, cmap=cmap, ax=ax, square=True, vmin=0, vmax=1, rasterized=True,
                              yticklabels=set_tags, xticklabels=set_tags, fmt='',
                              annot=None if not annot else annot_mat,
                              annot_kws={'size': font_s-2}, cbar_kws={'label': cbar_label, 'pad': 0.01})
    shared_heat.set_xticklabels(shared_heat.get_xmajorticklabels(), fontsize=font_s-2, rotation=90, fontweight='bold')
    shared_heat.set_yticklabels(shared_heat.get_ymajorticklabels(), fontsize=font_s-2, rotation=0, fontweight='bold')
    shared_heat_cbar = shared_heat.collections[0].colorbar
    shared_heat_cbar.ax.tick_params(labelsize=font_s)
    shared_heat_cbar.ax.yaxis.label.set_fontsize(font_s)
    plt.suptitle(title, y=0.93, size=font_s+4, fontweight='bold')
    plt.subplots_adjust(wspace=0.02)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        plt.savefig((plot_path + "_SharedHeatmap."+form).replace(" ", "_"), bbox_inches='tight', format=form)
    plt.close()


def volcano_plot(plot_df, x_col, y_col, mark_groups=None, mark_indexcol=None, mark_colours=None, mark_labels=None,
                 output_path="", title="", xsize=12, ysize=8, dot_size=2, top_labels=10, label_col="Gene name",
                 label_s=12, vlines=[], hlines=[], base_colour='black', log_y=True, formats=['pdf']):
    """
    Volcanoplot of x_col versus -log10(y_col). Usual application results from a differential expression analysis.

    Args:
        plot_df: pandas DataFrame from which all entries will be plotted.
        x_col: Column that will be used for the x-axis.
        y_col: Column used for the y-axis that will be -log10-transformed.
        mark_groups: Optional dictionary with {k: list/set} of identifiers that will be coloured differently.
        mark_indexcol: Optional column holding the index where to look for the identifiers from mark_groups.
        mark_colours: Optional dict with the colours assigned to the mark_groups, has to have the same keys as
            mark_groups. An existing colour palette will be used if not given.
        mark_labels: Optional legend labels for the mark_groups, leave empty to skip the legend.
        output_path: Path to store the pdf.
        title: Optional title to add.
        dot_size: Size of the dots of the scatterplots.
        top_labels: How many dots will be labelled, sorted by abs(x_col), if mark_groups is given will take
            top_labels from each of them, if not once from all.
        label_col: In which columns the text for the labels are found.
        label_s: Size of the label text.
        vlines: Positions of vertical lines to add.
        hlines: Positions of horizontal lines to add.
        base_colour: Colour for the dots that are not additionally marked.
    """
    if mark_groups and not mark_colours:
        mark_colours = {k: ColoursAndShapes.tol_vibrant[i] for i, k in enumerate(mark_groups.keys())}
    if mark_groups and not mark_labels:
        mark_labels = {k: None for k in mark_groups.keys()}
    if label_col and label_col not in plot_df.columns:  # Assume it's the index then.
        plot_df[label_col] = plot_df.index

    f, ax = plt.subplots(figsize=(xsize, ysize))
    ax.set_axisbelow(True)
    ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='major')
    plt.scatter(x=plot_df[x_col], y=plot_df[y_col].apply(lambda x: -np.log10(x)) if log_y else plot_df[y_col],
                s=dot_size, color=base_colour)
    if mark_groups:
        for m in mark_groups:
            plt.scatter(x=plot_df[plot_df[mark_indexcol].isin(mark_groups[m])][x_col],
                        y=plot_df[plot_df[mark_indexcol].isin(mark_groups[m])][y_col].apply(lambda x: -np.log10(x)) if log_y else
                        plot_df[plot_df[mark_indexcol].isin(mark_groups[m])][y_col],
                        s=dot_size, color=mark_colours[m], label=mark_labels[m])

    if not mark_colours:
        label_entries = plot_df.loc[plot_df.sort_values(by=x_col, key=abs, ascending=False).index[:top_labels]]
    else:
        label_entries = pd.DataFrame(columns=plot_df.columns)
        for m in mark_groups:
            top_m = plot_df.loc[plot_df[plot_df[mark_indexcol].isin(mark_groups[m])].sort_values(by=x_col, key=abs, ascending=False).index[:top_labels]]
            label_entries = label_entries.append(top_m)

    texts = []
    for i, entry in label_entries.iterrows():
        texts.append(plt.text(x=entry[x_col], y=-np.log10(entry[y_col]) if log_y else entry[y_col], s=entry[label_col],
                              fontsize=label_s))
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    plt.xlabel(x_col, fontsize=16)
    plt.ylabel("-log10("+y_col+")" if log_y else y_col, fontsize=16)
    for pos in hlines:
        plt.axhline(pos, color="#a7a8a7", linestyle="--")
    for pos in vlines:
        plt.axvline(pos, color="#a7a8a7", linestyle="--")
    ax.tick_params(axis='both', labelsize=14)
    if mark_groups:
        if np.all([x is not None for x in mark_labels.values()]):
            plt.legend(fontsize=14)
    plt.title(title, fontsize=18, fontweight='bold', y=1.05)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        f.savefig((output_path + "Volcano_" + str(x_col) + '_minuslog10' + str(y_col) + '.'+form).replace(' ', ''),
                  bbox_inches='tight', format=form)


def fisher_test_table(fisher_table, fisher_rows, fisher_cols, title='', output_path='', xsize=6, ysize=4,
                      formats=['pdf']):
    """
    Args:
        fisher_table: Nested list, usually [[Target GroupA, NonTarget GroupA], [Target GroupB, NonTarget GroupB]].
        fisher_rows: Names for the rows of the table, for the example [Group A, GroupB]
        fisher_cols: Names for the columns of the table, for the example [Target, NonTarget]
        title: Title for the table, followed by a newline and the p-value and oddsratio.
        output_path: Output path, the fisher_rows and cols will be added.
    """

    fish_stat, pval = scipy.stats.fisher_exact(fisher_table)
    if fish_stat == 0:
        print("WARNING: Oddsratio is 0")
    else:
        fish_stat = math.log2(fish_stat)

    f, ax = plt.subplots(figsize=(xsize, ysize))
    ax.axis('off')
    ks_table = plt.table(cellText=np.asarray(fisher_table), rowLabels=fisher_rows, colLabels=fisher_cols, loc='center',
                         cellLoc='center', rowLoc='center', zorder=12, bbox=[0, 0, 1, 1])
    ks_table.auto_set_font_size(False)
    ks_table.set_fontsize(16)
    plt.title(title + '\nFisher p-value: ' + str(round(pval, 5)) + '\nlog2(oddsratio): ' + str(round(fish_stat, 5)),
              fontsize=16)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        f.savefig(output_path+('FisherTable_'+'_'.join(fisher_rows+fisher_cols)).replace(' ', '')+'.'+form,
                  bbox_inches='tight', format=form)
    plt.close()
    return fish_stat, pval


def ks_test(curr_df, grouping_col, value_col):
    """
    Runs a Kolmogorov smirnov test for each non-redundant pairwise comparison. Tailored to the cumulative distribution
    function.
    """
    samples = sorted(curr_df[grouping_col].drop_duplicates().tolist())
    order_dict = {x: i for i, x in enumerate(samples)}
    comparisons = list(itertools.combinations(samples, 2))
    p_list = []
    for comp in comparisons:
        sample1_list = curr_df[curr_df[grouping_col] == comp[0]][value_col].values.tolist()
        sample2_list = curr_df[curr_df[grouping_col] == comp[1]][value_col].values.tolist()
        p_list.append([comp, scipy.stats.ks_2samp(sample1_list, sample2_list)[1]])
    p_nested = [[0 for x in range(y + 1)] for y in reversed(list(range(len(samples) - 1)))]
    p_background = [['white' for x in range(y + 1)] for y in reversed(list(range(len(samples) - 1)))]
    print(p_list)
    for i, p in enumerate(p_list):
        p_nested[order_dict[p[0][0]]][abs(order_dict[p[0][1]] - (len(samples) - 1))] = round(p[1], 5)
        if p[1] <= 0.05:
            p_background[order_dict[p[0][0]]][abs(order_dict[p[0][1]] - (len(samples) - 1))] = '#ed5739'
    for i, p in enumerate(p_nested):
        p_nested[i] += ['' for x in range((len(samples) - 1) - len(p_nested[i]))]
        p_background[i] += ['white' for x in range((len(samples) - 1) - len(p_background[i]))]

    return p_nested, p_background, samples


def cumulative_plot(plot_df, x_col, hue_col, hue_order=None, output_path='', numerate=False, title=None, xlimit=None,
                    vertical_line=False, add_all=False, table_width=0.3, table_x_pos=1.2, palette=None, font_s=16,
                    grid=True, formats=['pdf']):
    """
    Plot the cumulative distribution of all sets of grouping names in the plotting df.
    Adds a table with a Kolmogorov-Smirnov test for each non-redunant pairwise comparison next to the plot. Cells
    with a p-value ≤ 0.05 will be coloured in red.

    Args:
        plot_df: Pandas DataFrame holding the data.
        x_col: Column name in plot_df which should be plotted on the x-axis.
        hue_col: Column name in the DataFrame used for different curves.
        hue_order: If the groups in hue_col should have a specific order.
        numerate: If True, show the number of elements in each hue_col group in parentheses.
        vertical_line: To plot a vertical line at the given position, e.g. 0.
        add_all: If True, add the whole stack of values in x_col again as separate distribution labelled 'All'.
        table_width: Width of the KS table.
        table_x_pos: X-position of the KS table, to avoid it overlapping the plot.
        """
    p_df = plot_df.copy()  # Otherwise, updates for the enumeration would be referenced back to the function caller.
    if add_all:
        all_copy = p_df.copy()
        all_copy[hue_col] = 'All'
        p_df = pd.concat([all_copy, p_df], ignore_index=True)

    sns.set_style('ticks')
    sns.set_context('talk', font_scale=1.15)
    if hue_order:
        gene_set_cols = ['All']*add_all + hue_order
    else:
        gene_set_cols = p_df[hue_col].drop_duplicates().to_list()
    gene_set_cols_num = []
    different_groups = len(gene_set_cols)
    if different_groups == 3:
        linestyles = ['solid', 'dashed', 'dashdot', 'solid']
    else:
        linestyles = ['solid', 'dashed', 'solid', 'dashdot']

    if numerate:
        hue_order = []
        for i, x in enumerate(gene_set_cols):
            gene_set_cols_num.append(str(gene_set_cols[i]) + ' (n=' + str(sum(p_df[hue_col] == x)) + ')')
            hue_order.append(str(gene_set_cols[i]) + ' (n=' + str(sum(p_df[hue_col] == x)) + ')')
        for i, x in enumerate(gene_set_cols):
            p_df.loc[p_df[hue_col] == x, hue_col] = gene_set_cols_num[i]
    if not palette:
        if different_groups == 2:
            palette = ColoursAndShapes.two_contrasts[0]
        elif different_groups <= len(ColoursAndShapes.tol_vibrant):
            palette = ColoursAndShapes.tol_vibrant[:different_groups]
        else:
            palette = [to_hex(c) for c in ColoursAndShapes.glasbey_palettes['glasbey']][:different_groups]
    else:
        if 'glasbey' in palette:
            palette = [to_hex(c) for c in ColoursAndShapes.glasbey_palettes[palette]][:different_groups]
    cumu = sns.displot(p_df, x=x_col, hue=hue_col, kind='ecdf', hue_order=hue_order,
                       palette=palette, height=8, aspect=1.3, zorder=12)
    cumu.fig.subplots_adjust(top=0.93)
    cumu.fig.suptitle('Cumulative distribution of ' + x_col if not title else title, size=font_s + 6, x=0.28, y=1.01)
    cumu.ax.set_xlabel(x_col, size=font_s+4)
    cumu.ax.set_ylabel('Proportion', size=font_s+4)
    if xlimit:
        cumu.set(xlim=(xlimit[0], xlimit[1]))
    colours_order = []
    for a in range(len(cumu.ax.lines)):
        colours_order.append(cumu.ax.get_lines()[a].get_color().lower())
        cumu.ax.lines[a].set_linestyle(linestyles[a % len(linestyles)])
        cumu.ax.lines[a].set_linewidth('4')
    colour_dict = {x: i for i, x in enumerate(colours_order)}
    print(colour_dict)
    print(different_groups)
    for a in range(len(cumu.ax.lines)):
        if to_hex(cumu.legend.get_lines()[a].get_color()) in colour_dict:  # Legend has it even if there was no data.
            cumu.legend.get_lines()[a].set_linestyle(linestyles[colour_dict[to_hex(cumu.legend.get_lines()[a].get_color())] % len(linestyles)])
            cumu.legend.get_lines()[a].set_linewidth('3')
    for spine in cumu.ax.spines.values():  # Increase spline width.
        spine.set_linewidth(3)
    p_nested, p_background, samples = ks_test(p_df, hue_col, x_col)
    samples = [(str(x).split('(n=')[0][:int(len(str(x).split('(#')[0]) / 2)] + '\n' + str(x).split('(n=')[0][int(len(str(x).split('(#')[0]) / 2):]) if len(str(x).split('(n=')[0]) > 8 else x.split('(n=')[0]
               for x in samples]
    if len(samples) > 1:
        ks_table = plt.table(cellText=np.asarray(p_nested).T, cellColours=np.asarray(p_background).T,
                             rowLabels=list(reversed(samples))[:-1], colLabels=samples[:-1], loc='bottom right',
                             cellLoc='center', rowLoc='center', zorder=12, bbox=[table_x_pos, -0.11, table_width, 0.35])
        ks_table.auto_set_font_size(False)
        ks_table.set_fontsize(font_s-6)

    if grid:
        cumu.ax.grid(True, axis='both', color='#c7c7c7', linewidth=1, which='major')
    cumu.ax.set_facecolor('white')
    if vertical_line is not None:
        cumu.ax.axvline(x=vertical_line, color='#6e6e6e', zorder=1)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        cumu.savefig((output_path + x_col.replace(' ', '') + '_' + hue_col + '_CumulativeDistribution.' + form).replace(' ', '_'), bbox_inches='tight',
                     format=form)
    plt.close()




