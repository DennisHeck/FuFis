from pybedtools import BedTool
import numpy as np
import Various

"""Collections of functions for integrating different data on top of gABC interaction files."""


def peaks_fetch_col(base_regions, pattern, same_peaks=False, fetch_col='log2FC'):
    # CARE Currently only allows one fetch_col and not a list.
    """Take a bed-file or BedTools object with regions as peaks and intersect it with other peak files
    defined with pattern to get their fetch_col values in the dictionary regions."""
    fill_bed = BedTool(base_regions)
    fill_dict = {'\t'.join(x.fields[:3]): {} for x in fill_bed}
    if '*' in pattern:
        diff_files = Various.fn_patternmatch(pattern)
    else:  # Assume we have just one file.
        diff_files = {pattern: pattern.split('/')[-1]}
    print(diff_files)
    # Store whether an enhancer was differential in any condition.
    for diff_file, comp in diff_files.items():
        print(comp)
        diff_head = {x: i for i, x in enumerate(open(diff_file).readline().strip().split('\t'))}
        # Get a bed-object of the differential peaks to intersect with the enhancers, then get the average in case of
        # multiple overlaps. Works for both versions from DiffBind, with and w/o recentering on the summits.
        if not same_peaks:
            diff_peaks = []
            for entry in open(diff_file).readlines()[1:]:
                entry = entry.strip().split('\t')
                diff_peaks.append('\t'.join(entry))
            diff_peaks_inter = fill_bed.intersect(BedTool('\n'.join(diff_peaks), from_string=True), wo=True)
            # Collect the fetch_col of all peaks that intersect the enhancers.
            diff_hits = {x: [] for x in fill_dict}
            for inter in diff_peaks_inter:
                diff_hits['\t'.join(inter.fields[:3])] += [
                    float(inter.fields[len(fill_bed[0].fields) + diff_head[fetch_col]])]
            # And now take the average.
            for hit in diff_hits:
                fill_dict[hit][comp] = str(np.mean(diff_hits[hit]))
        else:
            for entry in open(diff_file).readlines()[1:]:
                entry = entry.strip().split('\t')
                fill_dict['\t'.join(entry[:3])][comp] = str(entry[diff_head[fetch_col]])
    return fill_dict, list(diff_files.values())


def peaks_peaks_overlap(peak_file, other_peak_file):
    """Based on two bed-file path or BedTools object returns a dictionary with
       {chr\tstart\tend: {other peaks}}."""
    if type(peak_file) == str:
        peak_dict = {'\t'.join(x.strip().split('\t')[:3]): set() for x in open(peak_file).readlines()}
    else:
        peak_dict = {'\t'.join(x.fields[:3]): set() for x in peak_file}
    peaks_inter = BedTool(peak_file).intersect(BedTool(other_peak_file), wo=True)
    other_start = len(BedTool(peak_file)[0].fields)
    for inter in peaks_inter:
        peak_dict['\t'.join(inter.fields[:3])].add('\t'.join(inter.fields[other_start:other_start+3]))
    return peak_dict


def peaks_promoter_overlap(peak_file, gtf_file, tss_type='all', gene_set=()):
    """Based on a bed-file path or BedTools object returns a dictionary with
       {chr\tstart\tend: {genes whose promoter overlap}} and one with {gene: {peaks}}."""
    promoter_bed = TSS_Fetcher.gene_window_bed(gtf_file=gtf_file, extend=200, tss_type=tss_type, merge=True,
                                               gene_set=gene_set)
    gene_dict = {x.fields[3]: set() for x in promoter_bed}
    if type(peak_file) == str:
        peak_dict = {'\t'.join(x.strip().split('\t')[:3]): set() for x in open(peak_file).readlines()}
    else:
        peak_dict = {'\t'.join(x.fields[:3]): set() for x in peak_file}
    promoter_inter = BedTool(peak_file).intersect(promoter_bed, wo=True)
    for inter in promoter_inter:
        peak_dict['\t'.join(inter.fields[:3])].add(inter.fields[-4])
        gene_dict[inter.fields[-4]].add('\t'.join(inter.fields[:3]))
    return peak_dict, gene_dict


def promoter_fetch_col(pattern, gtf_file, tss_type='all', gene_set=(), fetch_col='log2FC'):
    # CARE Currently only allows one fetch_col and not a list.
    """Based on a bed-file path or BedTools object returns a dictionary with
       {g: average value of fetch_col of intersecting regions.}.
       Uses peaks_fetch_col but maps the regions back to the genes and takes care of multiple promoters."""
    promoter_bed = TSS_Fetcher.gene_window_bed(gtf_file=gtf_file, extend=200, tss_type=tss_type, merge=True,
                                               gene_set=gene_set)
    prom_gene_map = {'\t'.join(x.fields[:3]): x.fields[3] for x in promoter_bed}
    fill_dict, fill_cols = peaks_fetch_col(promoter_bed, pattern, same_peaks=False, fetch_col=fetch_col)
    gene_fetch = {g: [] for g in set([x.fields[3] for x in promoter_bed])}
    for prom, val in fill_dict.items():
        hit_val = float(val[fill_cols[0]])
        if not np.isnan(hit_val):  # If not all promoter had a value we can still form the mean after excluding NaNs.
            gene_fetch[prom_gene_map[prom]].append(hit_val)
    gene_fetch = {g: np.mean(val) for g, val in gene_fetch.items()}
    return gene_fetch, fill_cols


def peaks_genebody_overlap(peak_file, gtf_file, gene_set=()):
    """Based on a bed-file path or BedTools object returns a dictionary with
       {gene: fraction of gene body overlapping with peak_file}."""
    genebody_bed = TSS_Fetcher.gene_body_bed(gtf_file=gtf_file, gene_set=gene_set)
    genebody_overlap = {x.fields[3]: 0 for x in genebody_bed}
    genebody_lengths = {x.fields[3]: x.length for x in genebody_bed}
    genebody_inter = genebody_bed.intersect(BedTool(peak_file).sort().merge(), wo=True)
    for inter in genebody_inter:
        genebody_overlap[inter.fields[3]] += int(inter.fields[-1])
    gene_dict = {g: genebody_overlap[g] / genebody_lengths[g] for g in genebody_overlap}
    return gene_dict

