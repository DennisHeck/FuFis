from pybedtools import BedTool
import gzip
from itertools import chain


def bed_to_length_dict(bedobj):
    """
    Takes a bedtools object with an identifier (Ensembl ID) in the fourth column, and sums the length of all
    regions of that gene.
    """
    length_dict = {}
    for entry in bedobj:
        if entry.fields[3] not in length_dict:
            length_dict[entry.fields[3]] = 0
        length_dict[entry.fields[3]] += entry.length
    return length_dict


def bed_to_feature_dict(bedobj):
    """
    Takes a bedtools object and returns a dictionary with the fourth column as key and the bedtools fields as value.
    Multiple occurrences of the identifier are appended.
    """
    feature_dict = {}
    for entry in bedobj:
        if entry.fields[3] not in feature_dict:
            feature_dict[entry.fields[3]] = []
        feature_dict[entry.fields[3]].append(entry.fields)
    return feature_dict


def gene_window_bed(gtf_file, extend=200, gene_set=set(), tss_type='5', dict_only=False, merge=False,
                    open_regions=False):
    """
    Based on a gtf file fetches all or the most 5' TSS for all genes, and returns a BedTool object with windows
    around the TSS, expanding by 'extend' in each direction, resulting in a total window size of 2*'extend'+1.
    Alternatively gives a dictionary with the TSS.
    The BedTools intervals will be 0-based, the TSS in the dictionary still 1-based like in the gtf-file.
    Care: removes the .-suffixes from all gene IDs.

    Args:
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        extend: number of base pairs to extend the TSS in each direction
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
        tss_type: "5" to get only the 5' TSS or "all" to get all unique TSS of all transcripts in the gtf-file
        dict_only: Returns a dictionary instead of a BedTool's object.
        merge: If True, merges all intersecting promoter of the same gene into one row in the BedTool's object.
        open_regions: Optional bed file or BedTools' object, only overlapping parts of promoters will be kept for the
            BedTool's object.
    """
    if tss_type == '5':
        identifier = 'gene'
    elif tss_type == 'all':
        identifier = 'transcript'
    if gtf_file.endswith('.gz'):
        file_opener = gzip.open(gtf_file, 'rt')
    else:
        file_opener = open(gtf_file)

    if gene_set:
        gene_set = set([g.split('.')[0] for g in gene_set])

    tss_locs = {}
    with file_opener as gtf_in:
        for entry in gtf_in:
            if not entry.startswith('#') and entry.split('\t')[2] == identifier:
                line = entry.strip().split('\t')
                # Some gene IDs are non-unique if they have a _PAR_Y version.
                if not line[8].split('gene_id "')[-1].split('";')[0].endswith("_PAR_Y"):
                    this_gene = line[8].split('gene_id "')[-1].split('";')[0].split('.')[0]
                    gene_name = line[8].split('gene_name "')[-1].split('";')[0]
                    if not gene_set or this_gene in gene_set or gene_name in gene_set:
                        if this_gene not in tss_locs:
                            tss_locs[this_gene] = {'chr': None, 'tss': set(), '#transcripts': 0}

                        tss_locs[this_gene]['chr'] = line[0]
                        tss_locs[this_gene]['name'] = gene_name
                        if line[6] == '+':
                            if identifier == 'gene' and (not tss_locs[this_gene]['tss'] or list(tss_locs[this_gene]['tss'])[0] > int(line[3])):
                                tss_locs[this_gene]['tss'] = {int(line[3])}
                            elif identifier == 'transcript':
                                tss_locs[this_gene]['tss'].add(int(line[3]))
                                tss_locs[this_gene]['#transcripts'] += 1
                            tss_locs[this_gene]['strand'] = '+'
                        if line[6] == '-':
                            if identifier == 'gene' and (not tss_locs[this_gene]['tss'] or list(tss_locs[this_gene]['tss'])[0] < int(line[4])):
                                tss_locs[this_gene]['tss'] = {int(line[4])}
                            elif identifier == 'transcript':
                                tss_locs[this_gene]['tss'].add(int(line[4]))
                                tss_locs[this_gene]['#transcripts'] += 1
                            tss_locs[this_gene]['strand'] = '-'

    if dict_only:
        return tss_locs

    promoter_bed = BedTool('\n'.join(chain(*[[vals['chr'] + '\t' + str(max([0, tss - int(extend) - 1])) + '\t' +
                                             str(tss + int(extend)) + '\t' + g + '\t.\t' + vals['strand'] for tss in vals['tss']]
                                             for g, vals in tss_locs.items()])), from_string=True)

    if open_regions and str(open_regions).lower() != "false":
        if type(open_regions) == str:
            open_regions = BedTool('\n'.join(['\t'.join(x.strip().split('\t')[:3]) for x
                                              in open(open_regions).readlines() if not x.startswith('#')]), from_string=True)
        promoter_bed = promoter_bed.intersect(open_regions)

    if merge:  # Flip the chr and geneID column to merge promoter of the same gene, and afterwards flip again.
        promoter_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0], x.fields[4], x.fields[5]]) for x in promoter_bed]), from_string=True).sort().merge(c=[4, 5, 6], o='distinct')
        promoter_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0], x.fields[4], x.fields[5]]) for x in promoter_bed]), from_string=True)

    return promoter_bed


def gene_body_bed(gtf_file, gene_set=set(), dict_only=False):
    """
    From a gtf-file fetches the gene bodies, meaning start-end for the entries labelled as 'gene'.

    Args:
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
        dict_only: Returns a dict {Ensembl ID: [chr, start, end, Ensembl ID, '.', strand]. Else, returns a bedtool-
            object with the coordinates per gene.
    """
    if gene_set:
        gene_set = set([g.split('.')[0] for g in gene_set])

    hits = []
    if gtf_file.endswith('.gz'):
        file_opener = gzip.open(gtf_file, 'rt')
    else:
        file_opener = open(gtf_file)

    with file_opener as gtf_in:
        for line in gtf_in:
            if not line.startswith('#') and line.split('\t')[2] == 'gene':
                # Some gene IDs are non-unique if they have a _PAR_Y version.
                if not line.strip().split('\t')[8].split('gene_id "')[-1].split('";')[0].endswith("_PAR_Y"):
                    line = line.strip().split('\t')
                    gene = line[8].split('gene_id "')[-1].split('"; ')[0].split('.')[0]
                    gene_name = line[8].split('gene_name "')[-1].split('"; ')[0].split('.')[0]
                    if not gene_set or gene in gene_set or gene_name in gene_set:
                        hits.append([line[0], line[3], line[4], gene, line[5], line[6]])

    if dict_only:
        return {x[3]: x for x in hits}
    else:
        return BedTool('\n'.join(['\t'.join(x) for x in hits]), from_string=True)


def gene_feature_bed(gtf_file, feature, gene_set=set(), dict_only=False, merge=True, length_only=False,
                     keep_strand=False):
    """
    On gene-level fetch a specific feature of a gene, e.g. all exons, matching the 3rd column in the gtf-file.

    Args:
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        feature: one of gtf's annotated regions: 'CDS', 'Selenocysteine', 'UTR', 'exon', 'gene', 'start_codon',
            'stop_codon', 'transcript'.
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
        merge: Merge the features per gene.
        dict_only: Return a dict instead.
        length_only: Return a dict with {gene: feature_len}, only makes sense with merge=True.
        keep_strand: Whether the strand should be kept as information. CARE: merging and keep_strand together is not
            properly tested.
    """
    if not merge and length_only:
        print("WARNING: Getting length without merging")

    if gene_set:
        gene_set = set([g.split('.')[0] for g in gene_set])

    hits = []
    if gtf_file.endswith('.gz'):
        file_opener = gzip.open(gtf_file, 'rt')
    else:
        file_opener = open(gtf_file)

    with file_opener as gtf_in:
        for line in gtf_in:
            if not line.startswith('#') and line.split('\t')[2] == feature:
                # Some gene IDs are non-unique if they have a _PAR_Y version.
                if not line.strip().split('\t')[8].split('gene_id "')[-1].split('";')[0].endswith("_PAR_Y"):
                    line = line.strip().split('\t')
                    gene = line[8].split('gene_id "')[-1].split('"; ')[0].split('.')[0]
                    gene_name = line[8].split('gene_name "')[-1].split('"; ')[0].split('.')[0]
                    if not gene_set or gene in gene_set or gene_name in gene_set:
                        if keep_strand:
                            hits.append([line[0], line[3], line[4], gene, line[5], line[6]])
                        else:
                            hits.append([line[0], line[3], line[4], gene])

    feature_bed = BedTool('\n'.join(['\t'.join(x) for x in hits]), from_string=True)

    if merge:  # Flip the chr and geneID column to merge features of the same gene, and afterwards flip again.
        if keep_strand:
            feature_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0], x.fields[4], x.fields[5]]) for x in feature_bed]), from_string=True).sort().merge(c=4, o='distinct')
            feature_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0], x.fields[4], x.fields[5]]) for x in feature_bed]), from_string=True)
        else:
            feature_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0]]) for x in feature_bed]), from_string=True).sort().merge(c=4, o='distinct')
            feature_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0]]) for x in feature_bed]), from_string=True)

    if length_only:
        return bed_to_length_dict(feature_bed)

    if dict_only:
        return bed_to_feature_dict(feature_bed)

    else:
        return feature_bed


def gene_introns(gtf_file, gene_set=set(), dict_only=False, length_only=False):
    """
    Gets the introns of annotated genes in a gtf-file, by taking the gene bodies and subtracting exons and UTRs.

    Args:
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
        dict_only: Return a dict instead.
        length_only: Return a dict with {gene: feature_len}, only makes sense with merge=True.
    """
    gene_bodies = gene_body_bed(gtf_file=gtf_file, gene_set=gene_set, dict_only=False)
    exons = gene_feature_bed(gtf_file=gtf_file, feature='exon', gene_set=gene_set, merge=True)
    utrs = gene_feature_bed(gtf_file=gtf_file, feature='UTR', gene_set=gene_set, merge=True)

    # Format the bed files to have the Ensembl ID in the first column and subtract exons and utrs.
    gene_bodies_pergene = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0]] + x.fields[4:]) for x in gene_bodies]), from_string=True)
    exons_pergene = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0]] + x.fields[4:]) for x in exons]), from_string=True)
    utrs_pergene = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0]] + x.fields[4:]) for x in utrs]), from_string=True)
    introns = gene_bodies_pergene.subtract(exons_pergene)
    introns = introns.subtract(utrs_pergene)
    introns = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0]] + x.fields[4:]) for x in introns]), from_string=True)

    if length_only:
        return bed_to_length_dict(introns)

    if dict_only:
        return bed_to_feature_dict(introns)

    return introns
