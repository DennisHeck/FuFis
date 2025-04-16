from collections import Counter
import pandas as pd
import numpy as np
import random
from timeit import default_timer as clock
import pybedtools
import sys
from pybedtools import BedTool
import gzip
import copy
import scipy.stats
import statsmodels.stats.multitest
import os
import fnmatch
import re
from itertools import chain
import GTF_Processing


def df_column_binner(df, col, num_bins, string_precision=0, lower_bound=None, upper_bound=None, tag='', numerate=True):
    """Takes a pandas DataFrame and returns it with a column with a string indicating to which bin the
    chosen column belongs to, including the size of that bin. When specifying lower and upper_bound, everything
    beyond those bounds will be put into a general group of < / > bound, and only the range between bounds will be
    binned. The tag is added as prefix to new column name."""
    df = copy.deepcopy(df)
    to_bin = df[col].values.tolist()
    if not lower_bound:
        lower_bound = min(to_bin)
    if not upper_bound:
        upper_bound = max(to_bin)
    bin_size = (upper_bound - lower_bound) / num_bins

    bin_idx = []
    for entry in to_bin:
        if entry < lower_bound:
            bin_idx.append([entry, 'lower'])
        elif entry >= upper_bound:
            bin_idx.append([entry, 'upper'])
        else:
            bin_idx.append([entry, (entry-lower_bound)//bin_size])

    bin_occs = Counter([x[1] for x in bin_idx])
    bin_strings = {k: str(round(k*bin_size+lower_bound, string_precision)) + '≤ ' + col + ' <' +
                      str(round((k+1)*bin_size+lower_bound, string_precision)) + (' (#' + str(val) + ')') * numerate
                   for k, val in bin_occs.items() if k not in ['lower', 'upper']}
    bin_strings['lower'] = str(round(lower_bound, string_precision)) + '> ' + col + (' (#' + str(len([x for x in to_bin if x < lower_bound])) + ')')*numerate
    bin_strings['upper'] = str(round(upper_bound, string_precision)) + '≤ ' + col + (' (#' + str(len([x for x in to_bin if x >= upper_bound])) + ')')*numerate

    bin_out = [[x[0], bin_strings[x[1]]] for x in bin_idx]
    df[tag + 'binned ' + col] = [x[1] for x in bin_out]
    return df


def df_fixed_binner(df, col, bin_borders, tag='', numerate=True, lower_open=False, upper_open=True):
    """Takes a pandas DataFrame and returns it with a column with a string indicating to which bin the
    chosen column belongs to, including the size of that bin.
    @param bin_borders: As opposed to the function above, this one expects all bin borders as list. E.g [0, 2, 4] will
    give the bins 0 ≤ x < 2, 2 ≤ x < 4; x ≥ 4. Must cover all entries in the df column.
    @param tag: The tag is added as prefix to new column name.
    @param numerate: Add the number of entries in each bin in brackets."""
    df = copy.deepcopy(df)
    if len(bin_borders) == 1:
        print("ERROR: requires more than one bin border.")
        sys.exit()
    if min(df[col]) < bin_borders[0]:
        print("ERROR: lowest bin border is larger then the minimum of the df entries.")
        return
    to_bin = df[col].values.tolist()
    bin_idx = []
    for entry in to_bin:
        for n in range(len(bin_borders)):
            if n == len(bin_borders) - 1:
                bin_idx.append([entry, n])
            elif bin_borders[n] <= entry < bin_borders[n + 1]:
                bin_idx.append([entry, n])
                break

    bin_occs = Counter([x[1] for x in bin_idx])
    bin_strings = {k: str(bin_borders[k]) + " ≤ " + col + (' (#' + str(bin_occs[k]) + ')') * numerate if k == len(bin_borders) - 1 else
                   str(bin_borders[k]) + ' ≤ ' + col + " < " + str(bin_borders[k+1]) + (' (#' + str(bin_occs[k]) + ')') * numerate
                   for k in range(len(bin_borders))}

    bin_out = [[x[0], bin_strings[x[1]]] for x in bin_idx]
    df[tag + 'binned ' + col] = [x[1] for x in bin_out]
    return df


def get_distance_to_one(start, end, other):
    """ Assumes that start - end is one region, and gives you the distance to other. other is supposed to be
     1-based, usually a TSS from a TSS file. Returns 0 if other is within the region."""
    distance = min([abs(other - 1 - start), abs(other - end)])  # For start subtract 1 to match 0-based.
    if start < other < end or end < other < start:
        distance = 0
    return distance


def tpm_norm(count_df, gtf_file):
    """
    Converts a count df into its TPM-normalized version. Genes for which no matching name/id can be found are removed
    before normalization.
    @param count_df: Pandas DataFrame with row indices as gene names or IDs and columns as samples.
    @param gtf_file: gtf-annotation file, either zipped or not.
    @return: tpm_df: The df with normalized counts
    @return: missed_genes: Indices that could not be found in the gtf file.
    """
    # Fetch the gene IDs and lengths.
    if gtf_file.endswith('.gz'):
        gtf_rows = gzip.open(gtf_file, 'rt').readlines()
    else:
        gtf_rows = open(gtf_file).readlines()

    name_id_map = {x.strip().split('\t')[8].split('gene_name "')[-1].split('";')[0]: x.strip().split('\t')[8].split('gene_id "')[-1].split('";')[0].split('.')[0]
                   for x in gtf_rows if not x.startswith("#") and x.split('\t')[2] == 'gene'}
    id_name_map = {x.strip().split('\t')[8].split('gene_id "')[-1].split('";')[0].split('.')[0]: x.strip().split('\t')[8].split('gene_name "')[-1].split('";')[0]
                   for x in gtf_rows if not x.startswith("#") and x.split('\t')[2] == 'gene'}
    # For gene length get the non-overlapping exons for each gene.
    exon_fetcher = ""
    gene_lengths = {}
    for row in [x.strip().split('\t') for x in gtf_rows if not x.startswith('#') and x.strip().split('\t')[2] == 'exon']:
        gene_id = row[8].split('gene_id "')[-1].split('";')[0].split('.')[0]
        gene_lengths[gene_id] = 0
        exon_fetcher += gene_id + '\t' + str(int(row[3]) - 1) + '\t' + row[4] + '\n'  # Gtf is 1-based, but bed 0-based.

    exon_bed = BedTool(exon_fetcher, from_string=True).sort().merge()
    for gene in exon_bed:
        gene_lengths[gene.fields[0]] += gene.length

    # First reduce to the genes we have the matching name/ID for.
    missed_genes = set(count_df.index[(~count_df.index.isin(id_name_map.keys())) & ~(count_df.index.isin(name_id_map.keys()))])
    count_df = count_df[(count_df.index.isin(id_name_map.keys())) | (count_df.index.isin(name_id_map.keys()))]

    gene_length_vec = [gene_lengths[name_id_map[g]] if g not in gene_lengths else gene_lengths[g] for g in count_df.index]
    read_scalings = {c: 0 for c in count_df.columns}
    for c in count_df.columns:
        read_scalings[c] = 10**6 / sum(count_df[c].values / gene_length_vec)

    tpm_df = pd.DataFrame(index=count_df.index)
    for c in count_df.columns:
        len_counts = count_df[c].values / gene_length_vec
        norm_counts = len_counts * read_scalings[c]
        tpm_df[c] = norm_counts

    return tpm_df, missed_genes


def open_genes_multibed(bed_dict, annotation):
    """For each of the bed files in bed_dict {tag: bed file} returns the set of genes that intersect with at least
    one of their promoters with a peak."""
    promoter_bed = GTF_Processing.gene_window_bed(annotation, extend=200, tss_type='all')
    open_genes = {}
    for tag, bed in bed_dict.items():
        print(bed)
        if type(bed) == str:
            bed = BedTool(bed)
        open_genes[tag] = set([x.fields[3] for x in promoter_bed.intersect(bed, u=True)])

    return open_genes


def open_genes(activity_file, annotation, first_col, tmp_dir):
    """
    NOTE: The GTF_Processing.gene_window_bed fulfils this function when not working with an activity column but just a bed file.
    Takes a bed-styled file with the activities of regions across multiple samples/conditions in columns.
    Takes all annotated promoters ±200bp for genes and checks which genes have an intersecting peak that has an
    activity > 0. The activity file requires a header preceded by #.
    @param open_genes: A dictionary where the keys are the column headers from the activity file and the
    values the set of genes with a nonzero peak at one of their promoters in that column.
    """
    pybedtools.set_tempdir(tmp_dir)
    start_prom = clock()
    promoter_bed = GTF_Processing.gene_window_bed(annotation, extend=200, tss_type='all', open_regions=activity_file)
    print(clock() - start_prom, "Got promoter bed")
    open_genes = {c: set() for c in open(activity_file).readline().strip().split('\t')[first_col:]}
    # Get the index of the columns in the bed intersection. The promoter columns are the first 4 entries.
    col_idx = {i: c for i, c in enumerate(open(activity_file).readline().strip().split('\t')[first_col:])}

    helper = 0
    for chro in ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']:
        print(chro)
        peak_nonzeros = {}
        start_chro = clock()
        bed_chro = ''
        with open(activity_file, 'r') as bed_in:
            bed_in.readline()  # Skip header.
            for entry in bed_in:
                if entry[:entry.find('\t')] == chro:
                    entry = entry.split('\t')
                    peak_nonzeros['\t'.join(entry[1:3])] = [c for c, val in enumerate(entry[first_col:]) if
                                                            float(val) > 0]
                    bed_chro += '\t'.join(entry[:3]) + '\n'
        print(clock() - start_chro, "bed fetched")
        chr_promoter = BedTool(''.join([str(x) for x in promoter_bed if x.fields[0] == chro]), from_string=True)
        for inter in chr_promoter.intersect(BedTool(bed_chro, from_string=True), wo=True):
            helper += 1
            if helper % 10000 == 0:
                print(helper)
            gene = inter.fields[3]
            for c in peak_nonzeros['\t'.join(inter.fields[5:7])]:
                open_genes[col_idx[c]].add(gene)
        print(clock() - start_chro)


def concat_df_dict(df_dict, new_col):
    """
    Takes a dictionary of {key: pd.DataFrame()} and concatenates them into one joint DataFrame. The key will be
    set as entry in a newly added column new_col.
    """
    joint_df = pd.DataFrame()
    for df_key, df_vals in df_dict.items():
        if df_vals is not None:
            df_vals[new_col] = df_key
            joint_df = pd.concat([joint_df, df_vals], axis=0, ignore_index=True)
    joint_df = joint_df[[new_col] + list(joint_df.columns[:-1])]
    return joint_df


def split_row_till_col(row_string, till_col, sep='\t'):
    """
    CARE: Test if it's faster
    Takes a string and splits it according to sep until a defined number of cols is reached. Uses .find() which is
    faster than .split() for long lines.
    @param row_string: Input string, should be one line and contain sep.
    @param till_col: Until which column index the row_string should be split (including till_col), e.g. 2 will give the
     first three columns. Everything after that is discarded.
    @param sep: The separater to be split on.
    @return: A list of split strings of length till_col.
    """

    split_cols = []
    for c in range(till_col+1):
        sep_idx = row_string.find(sep)
        split_cols.append(row_string[:row_string.find(sep)])
        row_string = row_string[sep_idx + 1:]

    return split_cols


def fn_patternmatch(pattern):
    """
    Grabs all files in the file system that match the pattern and returns a dictionary with {file: wildcard}.
    Only works if the wildcard is in the file name and not in a directory name.
    Files that start with ._ (e.g. temporary or system files) will be skipped.
    E.g. for a directory BirdCollection/ that contains Bird_Kakapo.txt and Bird_Kea.txt:
    fn_patternmatch("BirdCollection/*.txt" = {"BirdCollection/Bird_Kakapo.txt": "Kakao", "BirdCollection/Bird_Kea.txt": "Kea"}
     """
    parent_folder = '/'.join(pattern.split('/')[:-1]) + '/'
    children_pattern = pattern.split('/')[-1]
    re_pattern = re.compile(children_pattern.replace('*', '(.*)'))
    matched_files = {parent_folder + x: re_pattern.search(x).group(1)
                     for x in os.listdir(parent_folder)
                     if fnmatch.fnmatch(x, children_pattern) and not x.startswith('._')}
    return matched_files


def random_sampler_percentiled(df, var_col, idx_col, hit_list, iterations=1, perc_size=20, seed_int=1234):
    """
    Sample random sets of entries from a pandas DataFrame with a similar distribution of a second variable as a given
    list of entries.
    @param var_col: Numeric column in the DataFrame whose distribution should be matched.
    @param idx_col: A column with an identifier to sample from. If it's the index set to 'index'.
    @param hit_list: The list of identifiers which should be matched in their distribution of var_col by random
    identifiers.
    @param iterations: Number of random identifier sets.
    @param perc_size: Percentile size into which the distribution of var_col will be binned. The number of the hit_list
    in each list will be matched by random identifiers of the same bin.
    @return: Dictionary of {random iteration: set(random identifiers)}.
    """
    random.seed(seed_int)
    df = copy.deepcopy(df)  # Otherwise, we mutate the original df.

    if 100 % perc_size != 0:
        print("ERROR choose a percentile size that allows even splitting.")
        return
    if idx_col not in df.columns:  # Assume it's the index and paste it to a column.
        df[idx_col] = df.index

    perc_borders = [np.percentile(df[var_col].values, q=i*perc_size) for i in range(int(100/perc_size))]
    perc_borders += [df[var_col].max() + 1]  # The binning function would assume an open upper bound otherwise.

    df = df_fixed_binner(df, col=var_col, bin_borders=perc_borders, tag='', numerate=False)

    hit_brackets = Counter(df.loc[list(hit_list)]['binned '+var_col])

    random_sets = {}
    for rand_i in range(iterations):
        random_idx = []
        for bracket in hit_brackets:
            random_idx += random.sample(list(df[df['binned '+var_col] == bracket][idx_col]), k=hit_brackets[bracket])
        random_sets[rand_i] = random_idx

    return random_sets


def gene_cre_overlap_fisher(interaction_df, overlap_col, peak_id_col='peak_id'):
    """CARE: the input df is quite specific.
    Takes a pandas df of interactions (gene in Ensembl ID and a peak identifier in peak_id_col) and tests for each gene
    in the Df if its CREs more often have a hit in overlap_col (e.g. ChIP-seq peak overlap) than compared to the whole
    set of CREs that form interactions."""
    inter_genes = set(interaction_df['Ensembl ID'])
    inter_cres = len(set(interaction_df[peak_id_col]))
    gene_num_cres = Counter(interaction_df['Ensembl ID'])
    hits_only = interaction_df[interaction_df[overlap_col] > 0]
    cres_w_hits = len(set(hits_only[peak_id_col]))
    gene_cre_hits = Counter(hits_only['Ensembl ID'])
    gene_pvals = []
    for gene in inter_genes:
        fisher_table = [[gene_cre_hits[gene], gene_num_cres[gene] - gene_cre_hits[gene]],
                        [cres_w_hits - gene_cre_hits[gene], inter_cres - cres_w_hits - (gene_num_cres[gene] - gene_cre_hits[gene])]]
        fish_stat, pval = scipy.stats.fisher_exact(fisher_table, alternative='greater')
        gene_pvals.append([gene, pval])
    fisher_df = pd.DataFrame(gene_pvals, columns=['Ensembl ID', overlap_col+' p-value']).set_index('Ensembl ID')
    fisher_df[overlap_col+' FDR'] = statsmodels.stats.multitest.fdrcorrection(fisher_df[overlap_col+' p-value'],
                                                                              alpha=0.05, method='indep', is_sorted=False)[1]
    return fisher_df


