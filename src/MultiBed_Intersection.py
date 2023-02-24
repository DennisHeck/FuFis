from pybedtools import BedTool
import pandas as pd
import itertools
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import upsetplot


def inter_heatmap(bed_dict, metric, plot_path, annot_nums=True, inter_thresh=1e-9, x_size=16, y_size=10, annot_s=15):
    """
    Based on a dictionary {tag: bed_path/BedTools} creates an asymmetric heatmap of the intersection fo each feature,
    and adds an additional heatcolumn with the absolute number of regions in the file. annot_nums decides if the
    numbers of intersection should also be written in the cells.
    """
    bed_formatted = {t: BedTool(bed) if type(bed) is str else bed for t, bed in bed_dict.items()}
    bed_counts = pd.DataFrame([len(bed) for t, bed in bed_formatted.items()], index=bed_formatted.keys())
    bed_order = list(bed_formatted.keys())

    # Construct a matrix with the similarities of bed-regions as percentage.
    bed_shared_frac = np.ones([len(bed_formatted), len(bed_formatted)])
    bed_shared_abs = np.zeros([len(bed_formatted), len(bed_formatted)], dtype=object)  # Otherwise, seaborn still writes floats.
    for b, bed in enumerate(bed_order):
        bed_shared_abs[b][b] = len(bed_formatted[bed])

    bed_pairs = list(itertools.combinations(bed_formatted.keys(), 2))
    for p, pair in enumerate(bed_pairs):
        shared0 = len(bed_formatted[pair[0]].intersect(bed_formatted[pair[1]], u=True, f=inter_thresh))
        shared1 = len(bed_formatted[pair[1]].intersect(bed_formatted[pair[0]], u=True, f=inter_thresh))
        bed_shared_frac[bed_order.index(pair[0])][bed_order.index(pair[1])] = shared0 / len(bed_formatted[pair[0]])
        bed_shared_frac[bed_order.index(pair[1])][bed_order.index(pair[0])] = shared1 / len(bed_formatted[pair[1]])
        bed_shared_abs[bed_order.index(pair[0])][bed_order.index(pair[1])] = str(shared0)
        bed_shared_abs[bed_order.index(pair[1])][bed_order.index(pair[0])] = str(shared1)

    if annot_nums:
        annot_str = bed_shared_abs
    else:
        annot_str = np.full([len(bed_formatted), len(bed_formatted)], "", dtype=str)

    f, axes = plt.subplots(nrows=1, ncols=3, figsize=(x_size, y_size), gridspec_kw={'width_ratios': [0.2, 0.4, 8]})
    sns.heatmap(bed_counts, cmap="Blues", ax=axes[1], xticklabels=False, rasterized=True, yticklabels=False, cbar=False,
                fmt='', annot=bed_counts, annot_kws={'size': 12})
    cobar = f.colorbar(axes[1].get_children()[0], cax=axes[0], orientation="vertical", label='total number of '+metric)
    cobar.ax.yaxis.label.set_fontsize(14)
    cobar.ax.tick_params(labelsize=14)
    cobar.ax.yaxis.set_label_position('left')
    cobar.ax.yaxis.set_ticks_position('left')
    shared_heat = sns.heatmap(bed_shared_frac, cmap="Blues", ax=axes[2], square=True, vmin=0, vmax=1, rasterized=True,
                              yticklabels=bed_order, xticklabels=['shared w/ '+c for c in bed_order],
                              fmt='', annot=annot_str, annot_kws={'size': annot_s},
                              cbar_kws={'label': "Fraction "+metric.lower()+" shared", 'pad': 0.01})
    shared_heat.set_xticklabels(shared_heat.get_xmajorticklabels(), fontsize=14, rotation=90)
    shared_heat.set_yticklabels(shared_heat.get_ymajorticklabels(), fontsize=14, rotation=0)
    shared_heat_cbar = shared_heat.collections[0].colorbar
    shared_heat_cbar.ax.tick_params(labelsize=14)
    shared_heat_cbar.ax.yaxis.label.set_fontsize(14)
    plt.suptitle("Shared "+metric.lower(), y=0.93, size=20, fontweight='bold')
    plt.subplots_adjust(wspace=0.02)
    plt.savefig(plot_path + "_MultiIntersectHeat.pdf", bbox_inches='tight')
    plt.close()


def upset_to_reference(bed_files, ref_tag, y_label='Intersecting regions', title_tag='', plot_path=''):
    """
    Creates an upsetplot based on the bed_file defined by the ref_tag. That means that the bars will all show the
    number of reference regions, to avoid ambiguous many-to-one intersections. That also means that the horizontal bars
    that would normally show the total size of all other comparing bed-files will not be shown.
    """
    comparison_formatted = {t: BedTool(bed) if type(bed) is str else bed for t, bed in bed_files.items()}
    reference = comparison_formatted[ref_tag]
    # For the comparisons store the regions of the reference to do the upset intersection on.
    comparison_hits = {c: set(['\t'.join(x.fields[:3]) for x in reference.intersect(comparison_formatted[c], u=True)])
                       for c in
                       reversed(list(comparison_formatted.keys()))}  # Reverse to have insertion order in upset.

    intersection = upsetplot.from_contents(comparison_hits)
    fig = plt.figure(figsize=(10 + int(len(bed_files) / 2), 7 + int(len(bed_files) / 2)))
    upset = upsetplot.UpSet(intersection, show_counts=True, sort_categories_by=None, sort_by='cardinality',
                            element_size=None, intersection_plot_elements=int(len(bed_files) * 1.8),
                            totals_plot_elements=1,
                            facecolor="#010d4a", show_percentages=True)
    for c in range(len(comparison_hits)):  # The crude way of removing the horizontal bars.
        upset.totals[c] = 0
    with plt.rc_context({'axes.titlesize': 1,
                         'axes.labelsize': 14,
                         'xtick.labelsize': 18,  # Ensures enough space between horizontal bars and UpSet.
                         'ytick.labelsize': 12,
                         'font.size': 12}):  # For whatever reason font.size is only for the counts on top of the bars.
        upset.plot(fig)
    ax = plt.gca()
    ax.set_ylabel(y_label)
    ax.set_title("UpSet " + title_tag + '\n#' + str(len(intersection)), fontsize=16, fontweight='bold', y=1.15)
    plt.savefig(plot_path + '_UpSet_' + ref_tag + '.pdf', bbox_inches='tight')


