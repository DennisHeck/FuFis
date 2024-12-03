from pybedtools import BedTool
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import gzip
import TSS_Fetcher
import BasicPlotter


def gene_location_bpwise(gtf_file, bed_dict, plot_path, tss_type='5', external_bed={}, palette='tab20',
                         formats=['pdf']):
    """
    Based on a gtf file builds bed-objects for Promoter (±200bp), Exons, UTRs and Introns, and then counts how many
    of the bp in the bed file(s) are located within those annotations and which are intergenic. All gene features are
    exclusive, overlaps are removed. Introns are gene bodies subtracted by all other features. The bed_dict can also
    be a list of bed files, to omit recreating the gtf-annotations each time. Creates a pie chart with
    the percentages in the given path.
    Also returns a dictionary with the bp-location of each bed-region.
    @param bed_dict: A dictionary with {title tag: bed-file or BedTools object}
    @param external_bed: An additional dictionary of bed-file or BedTools object which will be added as category for the
    intersection. This will be considered as highest priority, meaning the regions in there are removed from the
    gene-related features, and a bp overlapping external_bed will not be counted anywhere else. Multiple external
    bed-regions shouldn't overlap, that causes undefined outcomes.
    @param: palette: Note with ax.pie this behaves weird and the list of colours has to be predefined. Gives unexpected
    palettes.
    """
    if gtf_file.endswith('.gz'):
        bed_annotations = {'Promoter': TSS_Fetcher.gene_window_bed(gtf_file, 200, tss_type=tss_type).sort().merge(),
                       'Exons': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in gzip.open(gtf_file, 'rt').readlines() if not x.startswith('#') and x.split('\t')[2] == 'exon']), from_string=True).sort().merge(),
                       'UTR': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in gzip.open(gtf_file, 'rt').readlines() if not x.startswith('#') and x.split('\t')[2] == 'UTR']), from_string=True).sort().merge()}
    else:
        bed_annotations = {'Promoter': TSS_Fetcher.gene_window_bed(gtf_file, 200, tss_type=tss_type).sort().merge(),
                       'Exons': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in open(gtf_file).readlines() if not x.startswith('#') and x.split('\t')[2] == 'exon']), from_string=True).sort().merge(),
                       'UTR': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in open(gtf_file).readlines() if not x.startswith('#') and x.split('\t')[2] == 'UTR']), from_string=True).sort().merge()}
    # Remove the Promoter regions and UTRs from Exons and Promoter from UTRs, want to have exclusive annotations.
    bed_annotations['Exons'] = bed_annotations['Exons'].subtract(bed_annotations['Promoter']).subtract(bed_annotations['UTR'])
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

    bed_inters = {}
    bed_locations = {}
    for tag, this_bed in bed_dict.items():
        if len(this_bed) == 0:
            print("Empty bed", tag)
            continue
        # Merge the bedfile to have unique bp, but keep track of what got merged and assemble it back later.
        bed = BedTool(this_bed).sort().merge(c=[1, 2, 3], o=['collapse']*3)
        org_beds = {'\t'.join(x.fields[:3]): [] for x in bed}
        for site in bed:
            org_beds['\t'.join(site.fields[:3])] += ['\t'.join([site.fields[3].split(',')[i], site[4].split(',')[i], site[5].split(',')[i]]) for i in range(site.fields[4].count(',')+1)]
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
            print(annot, round(overall_inter[annot] / sum(overall_inter.values()) * 100, 2), 'max', max_inter, round(max_inter / len(bed) * 100, 2))
            locations[annot] = round(overall_inter[annot] / sum(overall_inter.values()) * 100, 2)

        loc_df = pd.DataFrame([[annot, locations[annot]] for annot in overall_inter.keys()],
                              columns=['Location', 'Overlap']).set_index('Location')
        BasicPlotter.basic_pie(plot_df=loc_df, title=tag + '\n#' + str(len(BedTool(this_bed))), palette=palette,
                               numerate=False, legend_perc=True, formats=formats,
                               output_path=plot_path+tag+"_GeneFeatureLocation_bpwise", legend_title='')

        # Now map the potentially merged regions back to all its original regions.
        org_inters = {}
        for inter, locs in bed_inter.items():
            for sub_i in org_beds[inter]:
                org_inters[sub_i] = locs
        bed_inters[tag] = org_inters
        bed_locations[tag] = locations

    return bed_inters, bed_locations
