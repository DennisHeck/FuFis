import pybedtools
import itertools
import numpy as np
import upsetplot
import pickle
import os
import scipy.stats
from pybedtools import BedTool
from matplotlib import pyplot as plt
from multiprocessing import Pool
import seaborn as sns
import pandas as pd
import gzip
import src.GTF_Processing as GTF_Processing
import src.BasicPlotter as BasicPlotter
import src.Various as Various

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
        print(tag)
        if len(BedTool(this_bed)) == 0:
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


def intersection_heatmap_helper(args):
    """
    Helper function for intersection_heatmap that needs to be on the outermost level to allow parallelization.
    """
    pair, bed_formatted, inter_thresh, fisher, background_formatted, pseudocount = args
    if pair[0] == pair[1]:  # In case we have some appearing in the column and rows.
        shared = len(bed_formatted[pair[0]])
        return pair, shared, None, None
    else:
        print(pair)
        bed1 = bed_formatted[pair[0]]
        bed2 = bed_formatted[pair[1]]
        shared = len(bed1.intersect(bed2, u=True, f=inter_thresh))

    if fisher:
        bed1_background = background_formatted[pair[0]]
        background_shared = len(bed1_background.intersect(bed2, u=True, f=inter_thresh))
        fisher_table = [[shared + pseudocount, len(bed1) - shared + pseudocount],
                        [background_shared + pseudocount, len(bed1_background) - background_shared + pseudocount]]
        fish_stat, pval = scipy.stats.fisher_exact(fisher_table, alternative='two-sided')
    else:
        fish_stat = pval = None
    return pair, shared, fish_stat, pval


def intersection_heatmap(bed_dict, region_label, plot_path, fisher=False, fisher_background=None, pseudocount=1, annot_nums=True,
                         x_size=16, y_size=10, annot_s=15, col_beds=None, row_beds=None, inter_thresh=1e-9, vmin=None, vmax=None,
                         cmap='plasma', num_cmap='Blues', pickle_path='', tmp_dir='', wspace=0.1, hspace=0.1, n_cores=1,
                         width_ratios=[0.01, 0.05, 0.96], height_ratios=[0.05, 0.97], font_s=16, formats=['pdf']):
    """
    Based on a dictionary {tag: bed_path/BedTools} creates an asymmetric heatmap of the intersection fo each feature,
    and adds a heatcolumn with the absolute number of regions in the file. Also accepts different
    bed files/BedTools for rows and columns. The asymmetry of the heatmap overcomes the ambiguity of overlap when the
    regions are not complete subsets of one another. For example, one region from file A can contain multiple regions
    from file B, so one value for describing their overlap doesn't do it justice. The flag fisher allows to additionally run a
    Fisher's exact test for the enrichment of the regions versus the columns, but requires a background for each row to
    compare to. The cells are then coloured by the log2(oddsratio), if the test is significant (p-value ≤ 0.05).

    Args:
        bed_dict: {tag: bed_path/BedTools}
        region_label: label for the regions, e.g. 'peaks' or 'promoter'.
        fisher: If True, swap the mode of the function to calculate a two-sided Fisher's exact test instead of showing the
            fraction of overlap. Requires to set a background for each
            of the bed-dict entries used as row:
            [[row-regions with column overlap, row-regions without column overlap],
            [background-regions with column overlap, background-regions without column overlap]]
        fisher_background: A dictionary of {row_key: BedTool/bed-file} with the background on which the Fisher's exact
            test will be calculated for. Any regions in the background overlapping the original row regions is removed
            before running the test to ensure exclusivity of the groups.
        pseudocount: Pseudocount to be added to each cell in the contingency table for the Fisher's exact test.
        annot_nums: If the numbers of intersection should also be written in the cells, else hide.
        col_beds: A list of tags from bed_dict that will be shown in the columns. If not given will do all vs all.
        row_beds: Same as col_beds but for the rows. Either both or none have to be given.
        inter_thresh: Fraction of the region in the row that has to be overlapped by a region in the column to be counted
            as overlap.
        cmap: Colourmap for the intersection heatmap.
        num_cmap: Colourmap for the total number of regions.
        pickle_path: Path where a pickle object will be stored. If it already exists will be loaded instead to
            skip the intersection steps. If not given will store it in plot_path.
        tmp_dir: Optional path for temporary files created by pybedtools. CARE: If your outer function already uses
            the same tmp_dir for pybedtools, the objects in there will be deleted.
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
    if fisher and not fisher_background:
        print("ERROR: Can't run Fisher's exact test without a defined background.")
        return

    if not pickle_path:
        pickle_path = plot_path + "Intersection" + '_Fisher' * fisher + ".pkl"

    if not os.path.isfile(pickle_path):
        print("No pickle file found, doing the intersection.")
        intersection_storage = {}
        bed_formatted = {t: BedTool('\n'.join(['\t'.join(x.fields[:3]).replace('chr', '') for x in pybedtools.BedTool(val)]), from_string=True)
                         for t, val in bed_dict.items()}

        if not col_beds and not row_beds:
            col_beds = list(bed_formatted.keys())
            row_beds = list(bed_formatted.keys())

        if fisher:
            if not np.all([r in fisher_background for r in row_beds]):
                print("ERROR: Not all rows have a background in fisher_background")
                return
            # Prepare the background if needed, and remove any regions that overlap the respective row regions.
            background_formatted = {t: BedTool('\n'.join(['\t'.join(x.fields[:3]).replace('chr', '') for x in pybedtools.BedTool(val)]), from_string=True).intersect(bed_formatted[t], v=True)
                                    for t, val in fisher_background.items()}
        else:
            background_formatted = None

        # Storing the pandas Df with pickle could lead to an error due to a version conflict.
        intersection_storage['bed_rowcounts'] = [len(bed_formatted[bed]) for bed in row_beds]
        intersection_storage['bed_colcounts'] = [len(bed_formatted[bed]) for bed in col_beds]

        # Construct a matrix with the similarities of bed-regions as percentage.
        bed_shared_frac = np.ones([len(row_beds), len(col_beds)])
        bed_shared_abs = np.zeros([len(row_beds), len(col_beds)], dtype=object)  # Otherwise, seaborn still writes floats.
        bed_fisher_logodds = np.zeros([len(row_beds), len(col_beds)], dtype=float)

        bed_pairs = list(itertools.product(row_beds, col_beds))

        process_pool = Pool(processes=n_cores)
        pool_inter = process_pool.map(intersection_heatmap_helper, [[p, bed_formatted, inter_thresh, fisher, background_formatted, pseudocount] for p in bed_pairs])
        process_pool.close()

        for inter in pool_inter:
            pair, shared, fish_stat, pval = inter
            bed_shared_frac[row_beds.index(pair[0])][col_beds.index(pair[1])] = 0 if not len(bed_formatted[pair[0]]) else shared / len(bed_formatted[pair[0]])
            bed_shared_abs[row_beds.index(pair[0])][col_beds.index(pair[1])] = str(shared)
            if fish_stat:
                bed_fisher_logodds[row_beds.index(pair[0])][col_beds.index(pair[1])] = np.nan if (
                        pval > 0.05 or fish_stat == 0) else np.log2(fish_stat)
                print(pair, fish_stat, pval)
        intersection_storage['bed_shared_frac'] = bed_shared_frac
        intersection_storage['bed_shared_abs'] = bed_shared_abs
        intersection_storage['bed_fisher_logodds'] = bed_fisher_logodds
        intersection_storage['row_beds'] = row_beds  # To guarantee the same order from the object.
        intersection_storage['col_beds'] = col_beds
        pickle.dump(intersection_storage, open(pickle_path, 'wb'))
        if tmp_dir:
            pybedtools.helpers.cleanup(verbose=False, remove_all=False)
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

    for fisher_iteration in {False, fisher}:  # Only plot the fisher version when indicated.
        f, axes = plt.subplots(nrows=2, ncols=3, figsize=(x_size, y_size), gridspec_kw={'width_ratios': width_ratios,
                                                                                        'height_ratios': height_ratios})
        axes[0][0].remove()
        axes[0][1].remove()
        # Get the boundaries jointly for the row_ and col_counts.
        count_min = min([intersection_storage['bed_colcounts'].min().min(), intersection_storage['bed_rowcounts'].min().min()])
        count_max = max([intersection_storage['bed_colcounts'].max().max(), intersection_storage['bed_rowcounts'].max().max()])
        rowcount_heat = sns.heatmap(intersection_storage['bed_rowcounts'], cmap=num_cmap, ax=axes[1][1], xticklabels=False, rasterized=False,
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
            colcount_heat = sns.heatmap(intersection_storage['bed_colcounts'], cmap=num_cmap, ax=axes[0][2], xticklabels=False, rasterized=False,
                                        yticklabels=False, cbar=True, fmt='', annot=intersection_storage['bed_colcounts'], annot_kws={'size': font_s-2},
                                        square=False, cbar_kws={'label': "", 'pad': 0.01, 'shrink': 0})
            colcount_heat.collections[0].colorbar.ax.remove()  # Only needed it to align to the intersection heatmap.
        else:
            axes[0][2].remove()
        if fisher_iteration:
            shared_heat = sns.heatmap(intersection_storage['bed_fisher_logodds'], cmap='bwr', ax=axes[1][2], square=False,
                                      rasterized=True, center=0, yticklabels=row_beds, xticklabels=col_beds, fmt='',
                                      # annot=annot_str, annot_kws={'size': annot_s},
                                      cbar_kws={'label': "log$_{2}$(oddsratio) of row to column overlap", 'pad': 0.01})
        else:
            shared_heat = sns.heatmap(intersection_storage['bed_shared_frac'], cmap=cmap, ax=axes[1][2], square=False, vmin=vmin, vmax=vmax, rasterized=True,
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
        if not fisher_iteration:
            plt.suptitle("Shared " + region_label.lower(), y=0.93, size=font_s+2, fontweight='bold')
        else:
            plt.suptitle("Oddsratio "+ region_label.lower()+" compared to background", y=0.93, size=font_s + 2, fontweight='bold')
        plt.subplots_adjust(wspace=wspace, hspace=hspace)
        if type(formats) != list:
            formats = [formats]
        for form in formats:
            plt.savefig((plot_path + "_MultiIntersectHeat"+"_Fisher"*fisher_iteration+"."+form).replace(" ", '_'), bbox_inches='tight', format=form)
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


def peaks_peaks_overlap(peak_file, other_peak_file):
    """
    Intersects two bed-files or BedTools object to return a dictionary with {chr\tstart\tend: {other peaks}}.

    Args:
        peak_file: Path to a bed-file or BedTools object of the regions that will be the keys in the mapping dictionary.
        other_peak_file: Path or BedTools object to be mapped to peak_file and listed as values in the dictionary.

    Returns:
        dict:
            - **peak_dict**: A dictionary with {chr\tstart\tend: {other peaks}}. Values are empty if there was no intersection.
    """
    if type(peak_file) == str:
        peak_dict = {'\t'.join(x.strip().split('\t')[:3]): set() for x in open(peak_file).readlines()}
    else:
        peak_dict = {'\t'.join(x.fields[:3]): set() for x in peak_file}
    peaks_inter = BedTool(peak_file).intersect(BedTool(other_peak_file), wo=True)
    other_start = len(BedTool(peak_file)[0].fields)
    for inter in peaks_inter:
        peak_dict['\t'.join(inter.fields[:3])].add('\t'.join(inter.fields[other_start:other_start+3]))
    return peak_dict


def peaks_promoter_overlap(peak_file, gtf_file, tss_type='all', gene_set=(), extend=200):
    """
    Based on a bed-file path or BedTools object returns a dictionary with {chr\tstart\tend: {genes whose promoter overlap}}
    and one with {gene: {peaks}}.

    Args:
        peak_file: Path to a bed-file or BedTools object of the regions that will be used for the intersection.
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        tss_type: "5" to do the overlap only for the 5' TSS or "all" to do it for all unique TSS of all transcripts in the gtf-file.
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
        extend: Number of base pairs to extend the TSS in each direction. 200 means a window of size 401.

    Returns:
        tuple:
            - **peak_dict**: Dictionary mapping peaks to gene promoters {chr\tstart\tend: {genes whose promoter overlap}}.
            - **gene_dict**: Dictionary doing the mapping the other way around with {gene: {peaks}}.
    """
    promoter_bed = GTF_Processing.gene_window_bed(gtf_file=gtf_file, extend=extend, tss_type=tss_type, merge=True,
                                                  gene_set=gene_set)
    gene_dict = {x.fields[3]: set() for x in promoter_bed}
    peak_file = BedTool(peak_file)
    peak_dict = {'\t'.join(x.fields[:3]): set() for x in peak_file}
    promoter_inter = peak_file.intersect(promoter_bed, wo=True)
    for inter in promoter_inter:
        peak_dict['\t'.join(inter.fields[:3])].add(inter.fields[-4])
        gene_dict[inter.fields[-4]].add('\t'.join(inter.fields[:3]))
    return peak_dict, gene_dict


def peaks_fetch_col(base_regions, pattern, same_peaks=False, fetch_col='log2FC'):
    """
    Take a bed-file or BedTools object with regions as peaks and intersect it with other peak files
    defined with a filesystem pattern to get their fetch_col values in the base_regions. Useful for example, if one has
    multiple differential peak calls and want to map them to a base set of regions. If multiple regions overlap a base
    region, their average signal is taken. Entries of base peaks without overlap will be filled with NaN.

    Args:
        base_regions: BedTool's object or path to a bed file with the regions on which the intersection is centred on.
        pattern: File path pattern e.g., DiffPeaks/DiffBind_*.bed, or the path to just one individual file.
            All files matching that pattern will be used for the intersection. The string at the asterisk will be
            used as identifier. E.g., DiffPeaks/DiffBind_Macrophages.bed will be identified as Macrophages.
            The files need to have a header starting with # with fetch_col as column name. If it's not a path pattern
            but an individual file, the suffix after the last '/' is used as identifier.
        same_peaks: Boolean if the base regions are the same regions as the ones found in the pattern files.
        fetch_col: The column name in the files found by pattern to identify from which column the value should be taken.

    Returns:
        tuple:
            - **fill_dict**: Nested dictionary of the base_regions as keys and the average of the fetch_col in each of the matched files.
            - **matched_files**: List of the identifiers derived from the files matching the pattern.
    """
    fill_bed = BedTool(base_regions)
    fill_dict = {'\t'.join(x.fields[:3]): {} for x in fill_bed}
    if '*' in pattern:
        comparison_files = Various.fn_patternmatch(pattern)
    else:  # Assume we have just one file.
        comparison_files = {pattern: pattern.split('/')[-1]}
    print(comparison_files)
    # Store whether an enhancer was differential in any condition.
    for comp_file, comp in comparison_files.items():
        print(comp)
        comp_head = {x: i for i, x in enumerate(open(comp_file).readline().strip().split('\t'))}
        # Get a bed-object of the differential peaks to intersect with the enhancers, then get the average in case of
        # multiple overlaps. Works for both versions from DiffBind, with and w/o recentering on the summits.
        if not same_peaks:
            comp_peaks = []
            for entry in open(comp_file).readlines()[1:]:
                entry = entry.strip().split('\t')
                comp_peaks.append('\t'.join(entry))
            comp_peaks_inter = fill_bed.intersect(BedTool('\n'.join(comp_peaks), from_string=True), wo=True)
            # Collect the fetch_col of all peaks that intersect the enhancers.
            comp_hits = {x: [] for x in fill_dict}
            for inter in comp_peaks_inter:
                comp_hits['\t'.join(inter.fields[:3])] += [
                    float(inter.fields[len(fill_bed[0].fields) + comp_head[fetch_col]])]
            # And now take the average.
            for hit in comp_hits:
                fill_dict[hit][comp] = np.mean(comp_hits[hit]) if comp_hits[hit] else np.nan
        else:
            for entry in open(comp_file).readlines()[1:]:
                entry = entry.strip().split('\t')
                fill_dict['\t'.join(entry[:3])][comp] = float(entry[comp_head[fetch_col]])
    return fill_dict, list(comparison_files.values())


def promoter_fetch_col(pattern, gtf_file, tss_type='all', extend=200, gene_set=(), fetch_col='log2FC'):
    """
    Gets the promoter regions of genes and intersects them with peak files defined with a filesystem pattern to get
    their fetch_col values in the promoters. Useful for example, if one has
    multiple differential peak calls and wants to know the log2FC in promoters. If multiple regions overlap a promoter,
    their average signal is taken. Entries of promoters without overlap will be filled with NaN.

    Args:
        pattern: File path pattern e.g., DiffPeaks/DiffBind_*.bed, or the path to just one individual file.
            All files matching that pattern will be used for the intersection. The string at the asterisk will be
            used as identifier. E.g., DiffPeaks/DiffBind_Macrophages.bed will be identified as Macrophages.
            The files need to have a header starting with # with fetch_col as column name. If it's not a path pattern
            but an individual file, the suffix after the last '/' is used as identifier.
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        tss_type: "5" to get only the 5' TSS or "all" to get all unique TSS of all transcripts in the gtf-file. When 'all',
            values are averaged across multiple promoters.
        extend: Number of base pairs to extend the TSS in each direction. 200 means a window of size 401.
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
        fetch_col: The column name in the files found by pattern to identify from which column the value should be taken.

    Returns:
        tuple:
            - **gene_values**: Nested dictionary of the Ensembl IDs as keys and the average of the fetch_col in each of the matched files.
            - **matched_files**: List of the identifiers derived from the files matching the pattern.
        """
    promoter_bed = GTF_Processing.gene_window_bed(gtf_file=gtf_file, extend=extend, tss_type=tss_type, merge=True,
                                                  gene_set=gene_set)
    prom_gene_map = {'\t'.join(x.fields[:3]): x.fields[3] for x in promoter_bed}
    fill_dict, fill_cols = peaks_fetch_col(promoter_bed, pattern, same_peaks=False, fetch_col=fetch_col)
    gene_values = {g: {c: [] for c in fill_cols} for g in set([x.fields[3] for x in promoter_bed])}
    for prom, val in fill_dict.items():
        for f_col in fill_cols:
            hit_val = float(val[f_col])
            if not np.isnan(hit_val):  # If not all promoter had a value we can still form the mean after excluding NaNs.
                gene_values[prom_gene_map[prom]][f_col].append(hit_val)
    gene_values = {g: {c: np.mean(val[c]) if val[c] else np.nan for c in fill_cols} for g, val in gene_values.items()}
    return gene_values, fill_cols


def peaks_genebody_overlap(peak_file, gtf_file, gene_set=()):
    """
    Based on a bed-file path or BedTools object returns a dictionary with {gene: fraction of gene body overlapping with peak_file}.

    Args:
        peak_file: BedTool's object or path to a bed file with the regions on which the intersection is centred on.
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
   """
    genebody_bed = GTF_Processing.gene_body_bed(gtf_file=gtf_file, gene_set=gene_set)
    genebody_overlap = {x.fields[3]: 0 for x in genebody_bed}
    genebody_lengths = {x.fields[3]: x.length for x in genebody_bed}
    genebody_inter = genebody_bed.intersect(BedTool(peak_file).sort().merge(), wo=True)
    for inter in genebody_inter:
        genebody_overlap[inter.fields[3]] += int(inter.fields[-1])
    gene_dict = {g: genebody_overlap[g] / genebody_lengths[g] for g in genebody_overlap}
    return gene_dict



