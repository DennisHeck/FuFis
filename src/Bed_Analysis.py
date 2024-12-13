import pybedtools
import itertools
import numpy as np
import upsetplot
import pickle
import os
from pybedtools import BedTool
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import gzip
import src.GTF_Processing as GTF_Processing
import src.BasicPlotter as BasicPlotter

"""Collection of functions related to analysis of bed-files, their location or intersection with other files."""


def gene_location_bpwise(bed_dict, gtf_file, plot_path, tss_type='5', external_bed={}, palette='tab20',
                         formats=['pdf']):
    """
    Based on a gtf file builds bed-objects for Promoter (±200bp), Exons, UTRs and Introns, and then counts how many
    of the bp in the bed file(s) are located within those annotations and which are intergenic. All gene features are
    exclusive, overlaps are removed. Introns are gene bodies subtracted by all other features. The bed_dict can also
    be a list of bed files, to omit recreating the gtf-annotations each time. Creates a pie chart with
    the percentages in the given path. Also returns a dictionary with the bp-location of each bed-region and the
    total overlaps.

    Args:
        bed_dict: A dictionary with {title tag: bed-file or BedTools object}
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        tss_type: What to consider as promoter of genes, either '5' to use only the most 5' TSS, or 'all' to consider
            ±200bp around all annotated TSS of a gene.
        external_bed: An additional dictionary of bed-file or BedTools object which will be added as category for the
            intersection. This will be considered as highest priority, meaning the regions in there are removed from the
            gene-related features, and a bp overlapping external_bed will not be counted anywhere else. Multiple external
            bed-regions shouldn't overlap, that causes undefined outcomes.

    Returns:
        tuple:
            - **regions_locs**: For each entry in the bed_dict a dictionary with the bp-wise locations for each individual region.
            - **total_locs**: For each entry in the bed_dict the overall overlap of base pairs with each genomic feature.
    """
    if gtf_file.endswith('.gz'):
        bed_annotations = {'Promoter': GTF_Processing.gene_window_bed(gtf_file, 200, tss_type=tss_type).sort().merge(),
                           'Exons': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in
                                                       gzip.open(gtf_file, 'rt').readlines() if
                                                       not x.startswith('#') and x.split('\t')[2] == 'exon']),
                                            from_string=True).sort().merge(),
                           'UTR': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in
                                                     gzip.open(gtf_file, 'rt').readlines() if
                                                     not x.startswith('#') and x.split('\t')[2] == 'UTR']),
                                          from_string=True).sort().merge()}
    else:
        bed_annotations = {'Promoter': GTF_Processing.gene_window_bed(gtf_file, 200, tss_type=tss_type).sort().merge(),
                           'Exons': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in
                                                       open(gtf_file).readlines() if
                                                       not x.startswith('#') and x.split('\t')[2] == 'exon']),
                                            from_string=True).sort().merge(),
                           'UTR': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in
                                                     open(gtf_file).readlines() if
                                                     not x.startswith('#') and x.split('\t')[2] == 'UTR']),
                                          from_string=True).sort().merge()}
    # Remove the Promoter regions and UTRs from Exons and Promoter from UTRs, want to have exclusive annotations.
    bed_annotations['Exons'] = bed_annotations['Exons'].subtract(bed_annotations['Promoter']).subtract(
        bed_annotations['UTR'])
    bed_annotations['UTR'] = bed_annotations['UTR'].subtract(bed_annotations['Promoter'])

    introns_bed = BedTool(gtf_file).sort().merge()
    for remove in ['Promoter', 'Exons', 'UTR']:
        introns_bed = introns_bed.subtract(bed_annotations[remove])
    bed_annotations['Introns'] = introns_bed

    # Add potential external bed files and subtract it from all other regions.
    if external_bed:
        for external in external_bed:
            bed_annotations[external] = BedTool(external_bed[external]).sort().merge()
            for anno in bed_annotations:
                if anno not in external_bed:
                    bed_annotations[anno] = bed_annotations[anno].subtract(bed_annotations[external])

    regions_locs = {}
    total_locs = {}
    for tag, this_bed in bed_dict.items():
        if len(this_bed) == 0:
            print("Empty bed", tag)
            continue
        # Merge the bedfile to have unique bp, but keep track of what got merged and assemble it back later.
        bed = BedTool(this_bed).sort().merge(c=[1, 2, 3], o=['collapse'] * 3)
        org_beds = {'\t'.join(x.fields[:3]): [] for x in bed}
        for site in bed:
            org_beds['\t'.join(site.fields[:3])] += [
                '\t'.join([site.fields[3].split(',')[i], site[4].split(',')[i], site[5].split(',')[i]]) for i in
                range(site.fields[4].count(',') + 1)]
        bed_inter = {'\t'.join(x.fields[:3]): {a: 0 for a in bed_annotations.keys()} for x in bed}
        # Very broad marks rarely intersect only one annotation, so we count each bp.
        overall_inter = {a: 0 for a in list(bed_annotations.keys()) + ['Intergenic']}

        for annot, annot_bed in bed_annotations.items():
            intersection = bed.intersect(annot_bed, wo=True)
            for inter in intersection:
                bed_inter['\t'.join(inter.fields[:3])][annot] += int(inter.fields[-1])
                overall_inter[annot] += int(inter.fields[-1])

        # Intergenic is every base that was not assigned to any of the other annotations yet.
        for entry in bed_inter:
            missing_bp = abs(int(entry.split('\t')[2]) - int(entry.split('\t')[1])) - sum(bed_inter[entry].values())
            if missing_bp < 0:
                print("WARNING: There's a negative number of missing bp!")
            overall_inter['Intergenic'] += missing_bp
            bed_inter[entry]['Intergenic'] = missing_bp

        locations = {k: 0 for k in overall_inter.keys()}
        for annot in set(overall_inter.keys()):
            max_inter = len([x for x in bed_inter if max(bed_inter[x], key=bed_inter[x].get) == annot])
            print(annot, round(overall_inter[annot] / sum(overall_inter.values()) * 100, 2), 'max', max_inter,
                  round(max_inter / len(bed) * 100, 2))
            locations[annot] = round(overall_inter[annot] / sum(overall_inter.values()) * 100, 2)

        loc_df = pd.DataFrame([[annot, locations[annot]] for annot in overall_inter.keys()],
                              columns=['Location', 'Overlap']).set_index('Location')
        BasicPlotter.basic_pie(plot_df=loc_df, title=tag + '\n#' + str(len(BedTool(this_bed))), palette=palette,
                               numerate=False, legend_perc=True, formats=formats, legend_title='',
                               output_path=plot_path + (tag + "_GeneFeatureLocation_bpwise").replace(" ", '_'))

        # Now map the potentially merged regions back to all its original regions.
        org_inters = {}
        for inter, locs in bed_inter.items():
            for sub_i in org_beds[inter]:
                org_inters[sub_i] = locs
        regions_locs[tag] = org_inters
        total_locs[tag] = locations

    return regions_locs, total_locs


def inter_heatmap(bed_dict, region_label, plot_path, annot_nums=True, x_size=16, y_size=10, annot_s=15, col_beds=None,
                  row_beds=None, inter_thresh=1e-9, cmap='plasma', num_cmap='Blues', pickle_path='', tmp_dir='',
                  wspace=0.1, hspace=0.1, width_ratios=[0.01, 0.05, 0.96], height_ratios=[0.05, 0.97], font_s=16,
                  formats=['pdf']):
    """
    Based on a dictionary {tag: bed_path/BedTools} creates an asymmetric heatmap of the intersection fo each feature,
    and adds a heatcolumn with the absolute number of regions in the file. Also accepts different
    bed files/BedTools for rows and columns. The asymmetry of the heatmap overcomes the ambiguity of overlap when the
    regions are not complete subsets of one another. For example, one region from file A can contain multiple regions
    from file B, so one value for describing their overlap doesn't do it justice.

    Args:
        bed_dict: {tag: bed_path/BedTools}
        region_label: label for the regions, e.g. 'peaks' or 'promoter'.
        annot_nums: If the numbers of intersection should also be written in the cells, else hide.
        col_beds: A list of tags from bed_dict that will be shown in the columns. If not given will do all vs all.
        row_beds: Same as col_beds but for the rows. Either both or none have to be given.
        inter_thresh: Fraction of the region in the row that has to be overlapped by a region in the column to be counted
            as overlap.
        cmap: Colourmap for the intersection heatmap.
        num_cmap: Colourmap for the total number of regions.
        pickle_path: Path where a pickle object will be stored. If it already exists will be loaded instead to
            skip the intersection steps. If not given will store it in plot_path.
        tmp_dir: Optional path for temporary files created by pybedtools.
        wspace: Vertical whitespace between plot elements.
        hspace: Horizontal whitespace between plot elements.
        width_ratios: Control the width ratios of the sub-parts of the plot, meaning the colourbar for the total
            number of regions, the heatmap of the total number of regions and the heatmap of intersections.
        height_ratios: control the height ratios of the plot rows, the first one with the total number of regions
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
                                yticklabels=False, cbar=False, fmt='', annot=intersection_storage['bed_rowcounts'], annot_kws={'size': font_s-2},
                                square=False, vmin=count_min, vmax=count_max)
    cobar = f.colorbar(axes[1][1].get_children()[0], cax=axes[1][0], orientation="vertical",
                       label='Total number of ' + region_label)
    cobar.ax.yaxis.label.set_fontsize(font_s)
    cobar.ax.tick_params(labelsize=font_s)
    cobar.ax.yaxis.set_label_position('left')
    cobar.ax.yaxis.set_ticks_position('left')
    cobar.ax.ticklabel_format(style='plain')
    if row_beds != col_beds:
        colcount_heat = sns.heatmap(intersection_storage['bed_colcounts'], cmap=num_cmap, ax=axes[0][2], xticklabels=False, rasterized=True,
                                    yticklabels=False, cbar=True, fmt='', annot=intersection_storage['bed_colcounts'], annot_kws={'size': font_s-2},
                                    square=False, cbar_kws={'label': "", 'pad': 0.01, 'shrink': 0})
        colcount_heat.collections[0].colorbar.ax.remove()  # Only needed it to align to the intersection heatmap.
    else:
        axes[0][2].remove()
    shared_heat = sns.heatmap(intersection_storage['bed_shared_frac'], cmap=cmap, ax=axes[1][2], square=False, vmin=0, vmax=1, rasterized=True,
                              yticklabels=row_beds, xticklabels=col_beds,  # ['shared w/ ' + c for c in col_beds],
                              fmt='', annot=annot_str, annot_kws={'size': annot_s},
                              cbar_kws={'label': "Fraction " + region_label.lower() + " shared [(row & column) / row]",
                                        'pad': 0.01})
    shared_heat.set_xticklabels(shared_heat.get_xmajorticklabels(), fontsize=font_s, rotation=90)
    shared_heat.set_yticklabels(shared_heat.get_ymajorticklabels(), fontsize=font_s, rotation=0)
    shared_heat_cbar = shared_heat.collections[0].colorbar
    shared_heat_cbar.ax.tick_params(labelsize=font_s)
    shared_heat_cbar.ax.yaxis.label.set_fontsize(font_s)
    axes[1][2].xaxis.tick_top()
    plt.suptitle("Shared " + region_label.lower(), y=0.93, size=font_s+2, fontweight='bold')
    plt.subplots_adjust(wspace=wspace, hspace=hspace)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        plt.savefig(plot_path + "_MultiIntersectHeat."+form, bbox_inches='tight', format=form)
    plt.close('All')


def upset_to_reference(bed_files, ref_tag, y_label='Intersecting regions', title_tag='', plot_path='', font_s=16,
                       formats=['pdf']):
    """
    Creates an upsetplot based on the bed_file defined by the ref_tag. That means that the bars will all show the
    number of reference regions, to avoid ambiguous many-to-one intersections. That also means that the horizontal bars
    that would normally show the total size of all other comparing bed-files will not be shown (the ugly leftovers
    one the lower left can be cropped).

    Args:
        bed_dict: {tag: bed_path/BedTools}
        ref_tag: Key in bed_dict to use as the reference.
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
                         'axes.labelsize': font_s-2,
                         'xtick.labelsize': font_s+2,
                         'ytick.labelsize': font_s-4,
                         'font.size': font_s-4}):  # For whatever reason font.size is only for the counts on top of the bars.
        upset.plot(fig)
    ax = plt.gca()
    ax.set_ylabel(y_label)
    ax.set_title("UpSet " + title_tag + '\n#' + str(len(intersection)), fontsize=font_s, fontweight='bold', y=1.15)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        plt.savefig(plot_path + ('_UpSet_' + ref_tag + '.'+form).replace(' ', '_'), bbox_inches='tight', format=form)
    plt.close('All')

