from pybedtools import BedTool
from matplotlib import pyplot as plt
import seaborn as sns
import gzip
import TSS_Fetcher


def gene_location_bpwise(gtf_file, bed_file, plot_path, tss_type='5'):
    """
    Based on a gtf file builds bed-objects for promoter (±200bp), exons, UTRs and introns, and then counts how many
    of the bp in the bed file(s) are located within those annotations and which are intergenic. All gene features are
    exclusive, overlaps are removed. Introns are gene bodies subtracted by all other features. The bed_file can also
    be a list of bed files, to omit recreating the gtf-annotations each time. Creates a pie chart with
    the percentages in the given path.
    """
    if gtf_file.endswith('.gz'):
        bed_annotations = {'promoter': TSS_Fetcher.gene_window_bed(gtf_file, 200, tss_type=tss_type).sort().merge(),
                       'exons': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in gzip.open(gtf_file, 'rt').readlines() if not x.startswith('#') and x.split('\t')[2] == 'exon']), from_string=True).sort().merge(),
                       'utr': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in gzip.open(gtf_file, 'rt').readlines() if not x.startswith('#') and x.split('\t')[2] == 'UTR']), from_string=True).sort().merge()}
    else:
        bed_annotations = {'promoter': TSS_Fetcher.gene_window_bed(gtf_file, 200, tss_type=tss_type).sort().merge(),
                       'exons': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in open(gtf_file).readlines() if not x.startswith('#') and x.split('\t')[2] == 'exon']), from_string=True).sort().merge(),
                       'utr': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in open(gtf_file).readlines() if not x.startswith('#') and x.split('\t')[2] == 'UTR']), from_string=True).sort().merge()}
    # Remove the promoter regions and UTRs from exons and promoter from UTRs, want to have exclusive annotations.
    bed_annotations['exons'] = bed_annotations['exons'].subtract(bed_annotations['promoter']).subtract(bed_annotations['utr'])
    bed_annotations['utr'] = bed_annotations['utr'].subtract(bed_annotations['promoter'])

    introns_bed = BedTool(gtf_file).sort().merge()
    for remove in ['promoter', 'exons', 'utr']:
        introns_bed = introns_bed.subtract(bed_annotations[remove])
    bed_annotations['introns'] = introns_bed

    if type(bed_file) is not list:
        bed_file = [bed_file]
    for this_bed in bed_file:
        bed = BedTool(this_bed).sort().merge()
        bed_inter = {'\t'.join(x.fields[:3]): {a: 0 for a in bed_annotations.keys()} for x in bed}
        # Very broad marks rarely intersect only one annotation, so we count each bp.
        overall_inter = {a: 0 for a in list(bed_annotations.keys()) + ['intergenic']}

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
            overall_inter['intergenic'] += missing_bp
            bed_inter[entry]['intergenic'] = missing_bp

        locations = {k: 0 for k in overall_inter.keys()}
        for annot in set(overall_inter.keys()):
            max_inter = len([x for x in bed_inter if max(bed_inter[x], key=bed_inter[x].get) == annot])
            print(annot, round(overall_inter[annot] / sum(overall_inter.values()) * 100, 2), 'max', max_inter, round(max_inter / len(bed) * 100, 2))
            locations[annot] = round(overall_inter[annot] / sum(overall_inter.values()) * 100, 2)

        clrs = sns.color_palette('mako', n_colors=len(overall_inter)+2)[2:]  # List of RGB tuples, offset to skip darkest.
        f, ax = plt.subplots()
        if max(locations.values()) > 80:
            connector = ' '
        else:
            connector = '\n'
        ax.pie([locations[annot] for annot in overall_inter.keys()], labels=[annot + connector + str(locations[annot])+'%' for annot in overall_inter.keys()],
               autopct='', pctdistance=0.8, shadow=False, startangle=0, colors=clrs,
               labeldistance=1.1, textprops={'fontsize': 12})
        ax.axis('equal')
        if type(this_bed) == str:
            ax.set_title(this_bed.split('/')[-1].split('.')[0] + '\n#' + str(len(BedTool(this_bed))), fontsize=14, y=1.14)
            f.savefig(plot_path + this_bed.split('/')[-1].split('.')[0] + "_GeneFeatureLocation_bpwise.pdf",
                      bbox_inches='tight')
        else:
            ax.set_title('#' + str(len(BedTool(this_bed))), fontsize=14, y=1.12)
            f.savefig(plot_path + "_GeneFeatureLocation_bpwise.pdf",
                      bbox_inches='tight')
        plt.close()
