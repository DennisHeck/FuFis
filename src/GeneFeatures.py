import pandas as pd
from pybedtools import BedTool
from collections import Counter
import TSS_Fetcher
import Various


def gtf_gene_features(gtf_file, gene_set=(), extend=500000):
    """
    Based on a gtf-file fetch several features on gene-level.
    @param gtf_file: gtf-styled annotation, gzipped or uncompressed.
    @param gene_set: A set of genes to limit the output too, otherwise will use all present in the gtf-file.
    @param extend: How much the gene windows should be extended to each side of the TSS to calculate gene density.
    @return: A pandas DataFrame with the features, with the indices being either the gene_set or the gtf-file genes.
    """
    print("Fetching transcripts, size and biotype")
    tss_annot = TSS_Fetcher.gene_window_bed(gtf_file=gtf_file, gene_set=gene_set, tss_type='all', dict_only=True)
    gene_bodies = TSS_Fetcher.gene_body_bed(gtf_file=gtf_file, gene_set=gene_set, dict_only=True)
    gene_biotypes = Various.gene_biotypes(gtf_file=gtf_file, gene_set=gene_set)
    gene_exons = TSS_Fetcher.gene_feature_bed(gtf_file, feature='exon', gene_set=set(), merge=True, length_only=True)

    # For the gene density we need the whole annotation.
    print("Getting gene density")
    gene_windows = TSS_Fetcher.gene_window_bed(gtf_file=gtf_file, tss_type='5', extend=extend, dict_only=False)
    tss_bed = TSS_Fetcher.gene_window_bed(gtf_file=gtf_file, tss_type='5', dict_only=False, extend=0)
    gene_inter_counter = Counter([x.fields[3] for x in gene_windows.intersect(tss_bed, wa=True)])
    gene_inter_counter = {g: val - 1 for g, val in gene_inter_counter.items()}  # To correct the hit with itself.

    if not gene_set:
        df_genes = list(tss_annot.keys())
    else:
        df_genes = list(gene_set)

    gene_df = pd.DataFrame()
    gene_df.index = df_genes
    gene_df.index.name = 'Ensembl ID'
    gene_df['Gene name'] = [None if g not in tss_annot else tss_annot[g]['name'] for g in df_genes]
    gene_df['#TSS'] = [None if g not in tss_annot else len(tss_annot[g]['tss']) for g in df_genes]
    gene_df['#Transcripts'] = [None if g not in tss_annot else tss_annot[g]['#transcripts'] for g in df_genes]
    gene_df['Gene length'] = [None if g not in gene_bodies else int(gene_bodies[g][2]) - int(gene_bodies[g][1]) for g in df_genes]
    gene_df['Exons length'] = [None if g not in gene_exons else gene_exons[g] for g in df_genes]
    gene_df['Biotype gtf'] = [None if g not in gene_biotypes else gene_biotypes[g]['gtf'] for g in df_genes]
    gene_df['Biotype general'] = [None if g not in gene_biotypes else gene_biotypes[g]['general'] for g in df_genes]
    gene_df['Gene density'] = [None if g not in gene_inter_counter else gene_inter_counter[g] for g in df_genes]

    return gene_df


