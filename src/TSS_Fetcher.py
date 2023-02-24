from pybedtools import BedTool
import gzip
from itertools import chain


def gene_window_bed(gtf_file, extend=200, gene_set=set(), tss_type='5', dict_only=False, merge=False,
                    open_regions=False):
    """
    Based on a gtf file fetches all or the most 5' TSS for all genes, and returns a BedTool object with windows
    around the TSS, expanding by 'extend' in each direction, resulting in a total window size of 2*'extend'+1.
    Alternatively gives a dictionary with the TSS.
    The BedTools intervals will be 0-based, the TSS in the dictionary still 1-based like in the gtf-file.
    Care: removes the .-suffixes from all gene IDs.
    @param gtf_file: gtf-file in GENCODE's format, either .gz or .gtf
    @param extend: number of base pairs to extend the TSS in each direction
    @param gene_set: Limits the output to the given gene set, leave empty to get all.
    @param tss_type: "5" to get only the 5' TSS or "all" to get all unique TSS of all transcripts in the gtf-file
    @param dict_only: Returns a dictionary instead of a BedTool's object.
    @param merge: If True, merges all intersecting promoter of the same gene into one row in the BedTool's object.
    @param open_regions: Optional bed file or BedTools' object, only overlapping parts of promoters will be kept for the
                         BedTool's object.
    """
    if tss_type == '5':
        identifier = 'gene'
    elif tss_type == 'all':
        identifier = 'transcript'
    if gtf_file.endswith('.gz'):
        file_lines = gzip.open(gtf_file, 'rt').readlines()
    else:
        file_lines = open(gtf_file).readlines()

    if gene_set:
        gene_set = set([g.split('.')[0] for g in gene_set])

    tss_locs = {x.split('\t')[8].split('gene_id "')[-1].split('"; ')[0].split('.')[0]: {'chr': None, 'tss': set()}
                   for x in file_lines if not x.startswith('#') and x.split('\t')[2] == identifier}
    for line in [x.strip('\n').split('\t') for x in file_lines if not x.startswith('#')]:
        if line[2] == identifier:
            this_gene = line[8].split('gene_id "')[-1].split('";')[0].split('.')[0]
            gene_name = line[8].split('gene_name "')[-1].split('";')[0].split('.')[0]
            if gene_set:  # Skip and remove if not present.
                if this_gene not in gene_set and gene_name not in gene_set:
                    # Newest annotation has a duplicate entry for whatever reason, we might have removed it already.
                    if this_gene in tss_locs:
                        del tss_locs[this_gene]
                    continue
            tss_locs[this_gene]['chr'] = line[0]
            if line[6] == '+':
                if identifier == 'gene' and (not tss_locs[this_gene]['tss'] or list(tss_locs[this_gene]['tss'])[0] > int(line[3])):
                    tss_locs[this_gene]['tss'] = {int(line[3])}
                elif identifier == 'transcript':
                    tss_locs[this_gene]['tss'].add(int(line[3]))
            if line[6] == '-':
                if identifier == 'gene' and (not tss_locs[this_gene]['tss'] or list(tss_locs[this_gene]['tss'])[0] < int(line[4])):
                    tss_locs[this_gene]['tss'] = {int(line[4])}
                elif identifier == 'transcript':
                    tss_locs[this_gene]['tss'].add(int(line[4]))

    if dict_only:
        return tss_locs

    promoter_bed = BedTool('\n'.join(chain(*[[vals['chr'] + '\t' + str(max([0, tss - int(extend) - 1])) + '\t' +
                                             str(tss + int(extend)) + '\t' + g for tss in vals['tss']]
                                             for g, vals in tss_locs.items()])), from_string=True)

    if open_regions and str(open_regions).lower() != "false":
        promoter_bed = promoter_bed.intersect(open_regions)

    if merge:  # Flip the chr and geneID column to merge promoter of the same gene, and afterwards flip again.
        promoter_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0]]) for x in promoter_bed]), from_string=True).sort().merge(c=4, o='distinct')
        promoter_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0]]) for x in promoter_bed]), from_string=True)

    return promoter_bed


def gene_body_bed(gtf_file, gene_set=set()):

    if gene_set:
        gene_set = set([g.split('.')[0] for g in gene_set])

    hits = []
    if gtf_file.endswith('.gz'):
        file_lines = gzip.open(gtf_file, 'rt').readlines()
    else:
        file_lines = open(gtf_file).readlines()
    for line in [x.strip('\n').split('\t') for x in file_lines if not x.startswith('#') if x.strip().split('\t')[2] == 'gene']:
        gene = line[8].split('gene_id "')[-1].split('"; ')[0].split('.')[0]
        if not gene_set or gene in gene_set:
            hits.append([line[0], line[3], line[4], gene, line[5], line[6]])

    return BedTool('\n'.join(['\t'.join(x) for x in hits]), from_string=True)

