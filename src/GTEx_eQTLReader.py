import gzip
import pybedtools
from timeit import default_timer as clock
import TSS_Fetcher


"""Functions to read-in the GTEx files. Is fixed for the specific GTEx paths."""


def get_eqtls(gtex_folder, gtex_tissues, hg38_annotation, max_distance=None):
    """
    From the three fine-mapping methods provided on GTEx, go through the fixed file paths, and get the 
    eQTL-gene pairs matching the tissues queried. Optionally, set a maximum distance between eQTL and 
    the gene's TSS.

    Args:
        gtex_tissues: List of GTEx tissues (must match exactly), one entry can be a comma-separated list of GTEx tissues.

    Returns:
        tuple:
            - **eqtl_beds**: Dict of {tissues: {eqtl_type: bed object}} with the bed-object being all eQTL-gene pairs in format EnsemblID|POS-1|POS
    """
    start = clock()
    hg38_tss = TSS_Fetcher.gene_window_bed(hg38_annotation, tss_type='5', dict_only=True)
    eqtl_types = {'CAVIAR': gtex_folder+"/GTEx_v8_finemapping_CAVIAR/CAVIAR_Results_v8_GTEx_LD_HighConfidentVariants.gz",
                  'CaVEMaN': gtex_folder+"/GTEx_v8_finemapping_CaVEMaN/GTEx_v8_finemapping_CaVEMaN.txt.gz",
                  'DAP-G': gtex_folder+"/GTEx_v8_finemapping_DAPG/GTEx_v8_finemapping_DAPG.CS95.txt.gz"}

    eqtl_lists = {t: {e: [] for e in eqtl_types} for t in gtex_tissues}

    for eqtl in eqtl_types:
        print(eqtl)

        if eqtl == "CAVIAR":
            with gzip.open(eqtl_types[eqtl], 'rt') as caviar:
                caviar_head = {x: i for i, x in enumerate(caviar.readline().strip().split('\t'))}
                for entry in caviar:
                    entry = entry.strip().split('\t')
                    if entry[caviar_head['TISSUE']] in gtex_tissues:
                        this_eqtl = 'chr' + entry[caviar_head['eQTL']]
                        this_gene = entry[caviar_head['GENE']].split('.')[0]
                        if max_distance and not (this_gene in hg38_tss and \
                                    abs(int(this_eqtl.split('_')[1]) - list(hg38_tss[this_gene]['tss'])[0]) <= max_distance):
                            continue
                        eqtl_lists[entry[caviar_head['TISSUE']]][eqtl].append(this_gene + '\t' + this_eqtl)

        if eqtl == "CaVEMaN":
            with gzip.open(eqtl_types[eqtl], 'rt') as caveman:
                cave_head = {x: i for i, x in enumerate(caveman.readline().strip().split('\t'))}
                for entry in caveman:
                    entry = entry.strip().split('\t')
                    if entry[cave_head['TISSUE']] in gtex_tissues:
                        this_eqtl = '_'.join(entry[cave_head['eQTL']].split('_')[:2])
                        this_gene = entry[cave_head['GENE']].split('.')[0]
                        if max_distance and not (this_gene in hg38_tss and \
                                abs(int(this_eqtl.split('_')[1]) - list(hg38_tss[this_gene]['tss'])[0]) <= max_distance):
                            continue
                        eqtl_lists[entry[cave_head['TISSUE']]][eqtl].append(this_gene + '\t' + this_eqtl)

        if eqtl == "DAP-G":
            with gzip.open(eqtl_types[eqtl], 'rt') as dapg:
                for entry in dapg:
                    entry = entry.strip().split('\t')
                    this_eqtl = '_'.join(entry[2].split('_')[:2])
                    for sub_entry in entry[5].split('|'):
                        if sub_entry and sub_entry.split('@')[1].split('=')[0] in gtex_tissues:
                            this_gene = sub_entry.split('.')[0]
                            if max_distance and not (this_gene in hg38_tss and \
                                    abs(int(this_eqtl.split('_')[1]) - list(hg38_tss[this_gene]['tss'])[0]) <= max_distance):
                                continue
                            eqtl_lists[sub_entry.split('@')[1].split('=')[0]][eqtl].append(this_gene + '\t' + this_eqtl)


    eqtl_beds = {t: {} for t in gtex_tissues}
    tissue_genes = {t: {e: set() for e in eqtl_types} for t in gtex_tissues}  # To know the genes with eQTL in a tissue.
    for tissue in gtex_tissues:
        for eqtl in eqtl_types:
            eqtl_bed = pybedtools.BedTool('\n'.join(set([x.split('\t')[0] + '\t' + str(int(x.split('_')[1])-1) + '\t' + x.split('_')[1] for x in eqtl_lists[tissue][eqtl]])), from_string=True)
            eqtl_beds[tissue][eqtl] = eqtl_bed
            tissue_genes[tissue][eqtl] |= set([x.fields[0] for x in eqtl_bed])
    print('eQTL-beds collected', clock() - start)
    
    return eqtl_beds, tissue_genes

