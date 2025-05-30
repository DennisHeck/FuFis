from liftover import get_lifter
import os
import gzip
import numpy as np
import subprocess
from pybedtools import BedTool
import Various
import GTF_Processing


def genome_lifter(regions, input_version, output_version, same_chr=True, size_change=2):
    """
    Takes a list of regions and returns a list of all successful lifted regions,
    based on the converter object given. Strand is not considered. The chr can be w/ or w/o the 'chr' prefix,
    the original chr will be retained.

    Args:
        regions: Either the path to a bed-file (without header or with the first line starting with #), a BedTool's
            objects or a list of regions with [chr, start, end, ...]. Any additional columns are kept.
        input_version: Genome version of the entries in region_list. E.g., "hg19".
        output_version: Genome version region_list should be lifted to. E.g., "hg38".
        same_chr: Boolean if regions should only be lifted if they remain on the same chromosome. Yes, lifting may
            change the chromosome.
        size_change: Multiplier by how much the original region size is allowed to change. E.g., 2 means that a 100 bp
            region is allowed to be 200 bp after lifting, the region will be skipped otherwise.

    Returns:
        tuple:
            - **lifted_regions**: List of regions lifted to output_version.
            - **unliftable**: List of regions that could not be lifted with their coordinates still in the input genome version.
    """

    if type(regions) != list:
        region_list = [x.fields for x in BedTool(regions)]
    else:
        region_list = regions

    converter = get_lifter(input_version, output_version)
    lifted_regions = []
    unliftable = []
    for region in region_list:
        chro, start, end = region[:3]
        chro = str(chro)
        if '_' in chro:  # For those fun scaffold like chrUn_KI270467v1.
            unliftable += region
            continue
        try:
            new_start = converter[chro][int(float(start))]
            new_end = converter[chro][int(float(end))]
        except KeyError:  # In cases where we got the weird scaffolds.
            unliftable += region
            continue
        # The conversion failed if one position is not liftable and when the strands differ. We add as condition
        # that the lifted region does not exceed a size of twice the original region, and that it's on the same chr.
        if not (len(new_start) > 0 and len(new_end) > 0 and new_start[0][2] == new_end[0][2]):
            unliftable += region
            continue
        if same_chr:
            if not (chro.replace('chr', '') == new_start[0][0].replace('chr', '')
                    and chro.replace('chr', '') == new_end[0][0].replace('chr', '')):
                unliftable += region
                continue
        if abs(new_end[0][1] - new_start[0][1]) <= (size_change * (int(float(end)) - int(float(start)))):
            if new_end[0][1] > new_start[0][1]:  # The conversion sometimes maps to the minus strand, mixing locs.
                lifted_regions.append([chro, str(new_start[0][1]), str(new_end[0][1])] + region[3:])
            else:
                # On the minus strand we have to +1 and reverse start and end to match the UCSC lift.
                lifted_regions.append([chro, str(new_end[0][1] + 1), str(new_start[0][1] + 1)] + region[3:])
        else:
            unliftable += region

    return lifted_regions, unliftable


def abc_lifter(abc_folder, lift_folder, input_version, output_version, output_gtf, tss_mode="5"):
    """
    Takes a folder with ABC-called interactions and lifts them to another genome version. Restricted to genes
    that are present in both annotations, and where the interaction distance doesn't change >10,000.
    Expects naming scheme and format of STARE's ABCpp. Note that any version suffixes .X will be removed.
    @param abc_folder: directory with the ABC-scored interactions
    @param lift_folder: folder to write the output to
    @param input_version: genome annotation version of the ABC interactions, e.g. hg19
    @param output_version: genome annotation the interactions should be lifted to
    @param output_gtf: gtf annotation file of the version to lift to
    @param tss_mode: use the 5' TSS for distance ('5') or average the distance to all annotated TSS ('all')
    """

    if not os.path.isdir(lift_folder):
        os.mkdir(lift_folder)
    already_lifted = set([x.split('.txt.gz')[0] for x in os.listdir(lift_folder) if x.endswith('_'+output_version+'.txt.gz')])
    to_lift = set([x.split('.txt.gz')[0] for x in os.listdir(abc_folder) if not x.startswith('.') and '_ABCpp_scoredInteractions' in x and x.endswith(".txt.gz")]) - already_lifted

    # Get the 5' TSS in the output annotation to check how much the distance of interactions changed.
    output_tss = GTF_Processing.gene_window_bed(output_gtf, extend=1, tss_type=tss_mode, dict_only=True)

    for tag in to_lift:
        abc_file = abc_folder + '/' + tag + '.txt.gz'
        out_abc_file = lift_folder + '/' + tag + '_'+output_version+'.txt'
        print(abc_file)

        # Now we lift the interactions back to hg38, while keeping track of the TSS-dist.
        inter_head = {x: i for i, x in enumerate(gzip.open(abc_file, 'rt').readline().strip().split('\t'))}
        input_inter = [x.strip().split('\t') for x in gzip.open(abc_file, 'rt').readlines()[1:]]
        output_tagged = genome_lifter(input_inter, input_version=input_version, output_version=output_version)

        misses = 0
        gene_misses = 0
        missing_genes = set()
        with open(out_abc_file, 'w') as output:
            output.write(gzip.open(abc_file, 'rt').readline())
            for entry in output_tagged:
                this_gene = entry[inter_head['Ensembl ID']].split('.')[0]
                try:
                    new_dist = np.mean([Various.get_distance_to_one(start=int(entry[1]), end=int(entry[2]), other=tss) for tss in output_tss[this_gene]['tss']])
                    if abs(new_dist - float(entry[inter_head['TSS-dist']])) <= 10000 and output_tss[this_gene]['chr'].replace('chr', '') == entry[0].replace('chr', ''):
                        entry[inter_head['Ensembl ID']] = this_gene
                        entry[inter_head['TSS-dist']] = str(new_dist)
                        output.write('\t'.join(entry) + '\n')
                    else:
                        misses += 1
                except KeyError:
                    gene_misses += 1
                    missing_genes.add(this_gene)
        print('Distance mismatch', misses / len(output_tagged) * 100, '%')
        print('Gene mismatch', gene_misses / len(output_tagged) * 100, '%')
        subprocess.call('gzip ' + out_abc_file, shell=True)


