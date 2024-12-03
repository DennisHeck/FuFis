import pybedtools
import pandas as pd
import itertools
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import upsetplot
import pickle
import os


def inter_heatmap(bed_dict, region_label, plot_path, annot_nums=True, inter_thresh=1e-9, x_size=16, y_size=10,
                  annot_s=15, col_beds=None, row_beds=None, cmap='plasma', num_cmap='Blues', wspace=0.1, hspace=0.1,
                  pickle_path='', tmp_dir='', width_ratios=[0.01, 0.05, 0.96], height_ratios=[0.05, 0.97]):
    """
    Based on a dictionary {tag: bed_path/BedTools} creates an asymmetric heatmap of the intersection fo each feature,
    and adds a heatcolumn with the absolute number of regions in the file. Also accepts different
    bed files/BedTools for rows and columns.
    @param bed_dict: {tag: bed_path/BedTools}
    @param region_label: label for the regions, e.g. 'peaks' or 'promoter'
    @param annot_nums: if the numbers of intersection should also be written in the cells
    @param col_beds: a list of tags from bed_dict that will be shown in the columns. If not given will do all vs all.
    @param row_beds: same as col_beds but for the rows. Either both or none have to be given.
    @param cmap: colourmap for the intersection
    @param num_cmap: colourmap for the total number of regions.
    @param pickle_path: Path where a pickle object will be stored. If it already exists will be loaded instead to
    skip the intersection steps.If not given will store it in plot_path.
    @param width_ratios: control the width ratios of the sub-parts of the plot, meaning the colourbar for the total
    number of regions, the heatmap of the total number of regions and the heatmap of intersections.
    @param height_ratios: control the height ratios of the plot rows, the first one with the total number of regions
    of the column beds and the second with the rest.
    """
    if tmp_dir:
        pybedtools.helpers.set_tempdir(tmp_dir)

    if (not col_beds and row_beds) or (col_beds and not row_beds):
        print("ERROR: For different axes both col_beds and row_beds are required")
        return

    if not pickle_path:
        pickle_path = plot_path + "Intersection.pkl"

    if not os.path.isfile(pickle_path):
        print("No pickle file found, doing the intersection.")
        intersection_storage = {}
        bed_formatted = {t: pybedtools.BedTool(bed) if type(bed) is str else bed for t, bed in bed_dict.items()}

        if not col_beds and not row_beds:
            col_beds = list(bed_formatted.keys())
            row_beds = list(bed_formatted.keys())
        # Storing the pandas Df with pickle could lead to an error due to a version conflict.
        intersection_storage['bed_rowcounts'] = [len(bed_formatted[bed]) for bed in row_beds]
        intersection_storage['bed_colcounts'] = [len(bed_formatted[bed]) for bed in col_beds]

        # Construct a matrix with the similarities of bed-regions as percentage.
        bed_shared_frac = np.ones([len(row_beds), len(col_beds)])
        bed_shared_abs = np.zeros([len(row_beds), len(col_beds)], dtype=object)  # Otherwise, seaborn still writes floats.

        bed_pairs = list(itertools.product(row_beds, col_beds))
        for p, pair in enumerate(bed_pairs):
            if pair[0] == pair[1]:  # In case we have some appearing in the column and rows.
                shared = len(bed_formatted[pair[0]])
            else:  # Prevent issues with different number of columns.
                print(pair)
                bed1 = pybedtools.BedTool('\n'.join(['\t'.join(x.fields[:3]).replace('chr', '') for x in bed_formatted[pair[0]]]), from_string=True)
                bed2 = pybedtools.BedTool('\n'.join(['\t'.join(x.fields[:3]).replace('chr', '') for x in bed_formatted[pair[1]]]), from_string=True)
                shared = len(bed1.intersect(bed2, u=True, f=inter_thresh))
                if tmp_dir:
                    pybedtools.helpers.cleanup(verbose=False, remove_all=False)
            bed_shared_frac[row_beds.index(pair[0])][col_beds.index(pair[1])] = 0 if not len(bed_formatted[pair[0]]) else shared / len(bed_formatted[pair[0]])
            bed_shared_abs[row_beds.index(pair[0])][col_beds.index(pair[1])] = str(shared)
        intersection_storage['bed_shared_frac'] = bed_shared_frac
        intersection_storage['bed_shared_abs'] = bed_shared_abs
        intersection_storage['row_beds'] = row_beds  # To guarantee the same order from the object.
        intersection_storage['col_beds'] = col_beds
        pickle.dump(intersection_storage, open(pickle_path, 'wb'))
    else:
        print("Loading existing pickle file")
        intersection_storage = pickle.load(open(pickle_path, 'rb'))

    row_beds = intersection_storage['row_beds']
    col_beds = intersection_storage['col_beds']
    intersection_storage['bed_rowcounts'] = pd.DataFrame(intersection_storage['bed_rowcounts'], index=row_beds)
    intersection_storage['bed_colcounts'] = pd.DataFrame(intersection_storage['bed_colcounts'], index=col_beds).T

    if annot_nums:
        annot_str = intersection_storage['bed_shared_abs']
    else:
        annot_str = np.full([len(col_beds), len(row_beds)], "", dtype=str)

    f, axes = plt.subplots(nrows=2, ncols=3, figsize=(x_size, y_size), gridspec_kw={'width_ratios': width_ratios,
                                                                                    'height_ratios': height_ratios})
    # outer = gridsp ec.GridSpec(nrows=2, ncols=2, figure=f, width_ratios=[0.06, 0.96], wspace=0)
    # gs1 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec = outer[1, 0])
    # gs2 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec = outer[1, 1], hspace = .05)
    axes[0][0].remove()
    axes[0][1].remove()
    # Get the boundaries jointly for the row_ and col_counts.
    count_min = min([intersection_storage['bed_colcounts'].min().min(), intersection_storage['bed_rowcounts'].min().min()])
    count_max = max([intersection_storage['bed_colcounts'].max().max(), intersection_storage['bed_rowcounts'].max().max()])
    rowcount_heat = sns.heatmap(intersection_storage['bed_rowcounts'], cmap=num_cmap, ax=axes[1][1], xticklabels=False, rasterized=True,
                                yticklabels=False, cbar=False, fmt='', annot=intersection_storage['bed_rowcounts'], annot_kws={'size': 12},
                                square=False, vmin=count_min, vmax=count_max)
    cobar = f.colorbar(axes[1][1].get_children()[0], cax=axes[1][0], orientation="vertical",
                       label='total number of ' + region_label)
    cobar.ax.yaxis.label.set_fontsize(14)
    cobar.ax.tick_params(labelsize=14)
    cobar.ax.yaxis.set_label_position('left')
    cobar.ax.yaxis.set_ticks_position('left')
    cobar.ax.ticklabel_format(style='plain')
    if row_beds != col_beds:
        colcount_heat = sns.heatmap(intersection_storage['bed_colcounts'], cmap=num_cmap, ax=axes[0][2], xticklabels=False, rasterized=True,
                                    yticklabels=False, cbar=True, fmt='', annot=intersection_storage['bed_colcounts'], annot_kws={'size': 12},
                                    square=False, cbar_kws={'label': "", 'pad': 0.01, 'shrink': 0})
        colcount_heat.collections[0].colorbar.ax.remove()  # Only needed it to align to the intersection heatmap.
    else:
        axes[0][2].remove()
    shared_heat = sns.heatmap(intersection_storage['bed_shared_frac'], cmap=cmap, ax=axes[1][2], square=False, vmin=0, vmax=1, rasterized=True,
                              yticklabels=row_beds, xticklabels=col_beds,  # ['shared w/ ' + c for c in col_beds],
                              fmt='', annot=annot_str, annot_kws={'size': annot_s},
                              cbar_kws={'label': "Fraction " + region_label.lower() + " shared [(row & column) / row]",
                                        'pad': 0.01})
    shared_heat.set_xticklabels(shared_heat.get_xmajorticklabels(), fontsize=14, rotation=90)
    shared_heat.set_yticklabels(shared_heat.get_ymajorticklabels(), fontsize=14, rotation=0)
    shared_heat_cbar = shared_heat.collections[0].colorbar
    shared_heat_cbar.ax.tick_params(labelsize=14)
    shared_heat_cbar.ax.yaxis.label.set_fontsize(14)
    axes[1][2].xaxis.tick_top()
    plt.suptitle("Shared " + region_label.lower(), y=0.93, size=20, fontweight='bold')
    plt.subplots_adjust(wspace=wspace, hspace=hspace)
    plt.savefig(plot_path + "_MultiIntersectHeat.pdf", bbox_inches='tight')
    plt.close()


def upset_to_reference(bed_files, ref_tag, y_label='Intersecting regions', title_tag='', plot_path=''):
    """
    Creates an upsetplot based on the bed_file defined by the ref_tag. That means that the bars will all show the
    number of reference regions, to avoid ambiguous many-to-one intersections. That also means that the horizontal bars
    that would normally show the total size of all other comparing bed-files will not be shown.
    """
    comparison_formatted = {t: pybedtools.BedTool(bed) if type(bed) is str else bed for t, bed in bed_files.items()}
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


