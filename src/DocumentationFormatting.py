import pandas as pd
import Various

"""Snippets of code to format files used in the documentation as examples."""

# ---------------------------------------------------------------------------------------------------
# Bed_Analysis.intersection_heatmap
# ---------------------------------------------------------------------------------------------------
inter_heat_out = '/Users/dennis/GitHub/FuFis/ExampleData/MultiBedFisher/'
efo_map = {"EFO_0003777": 'heart disease',
           'EFO_0004326': 'heart rate'}
gwas_files = Various.fn_patternmatch("/Users/dennis/GitHub/FuFis/Unformatted_ExampleData/gwas-association-downloaded_2025-01-08-*-withChildTraits.tsv.gz")

for gwas in gwas_files:
    gwas_df = pd.read_table(gwas, sep='\t', header=0)
    gwas_df = gwas_df[~gwas_df['CHR_ID'].isna()]
    if gwas_df['CHR_ID'].dtype == float:
        gwas_df['CHR_ID'] = gwas_df['CHR_ID'].astype(int)
    gwas_df = gwas_df[gwas_df['CHR_ID'].astype(str) == '1']
    gwas_df['chr'] = ['chr' + str(int(val)) for val in gwas_df['CHR_ID'].values]
    gwas_df['start'] = gwas_df['CHR_POS'].astype(int)
    gwas_df['end'] = gwas_df['start'] + 1
    gwas_df[['chr', 'start', 'end']].to_csv(inter_heat_out + 'GWAS_'+gwas_files[gwas] + '_' + efo_map[gwas_files[gwas]].replace(' ', '-') + '.txt', header=False, sep='\t', index=False)

hocker_df = pd.read_table('/Users/dennis/Dev/STARE_GAZE/DataHocker/GSE165837_Union_peak_set_RPKM.tsv', sep='\t', header=0)
hocker_df = hocker_df[hocker_df['# chr'] == 'chr1']
hocker_df.columns = [c.replace('_RPKM', '') for c in hocker_df.columns]
hocker_df[['# chr', 'start', 'end']].to_csv(inter_heat_out + "Hocker2021_ATAC_Union_chr1.bed", sep='\t', header=False, index=False)
for cell in ['vCM', 'MAC', 'LC']:
    cell_df = hocker_df[hocker_df[cell] >= 2]
    cell_df[['# chr', 'start', 'end']].to_csv(inter_heat_out + "Hocker2021_ATAC_" + cell+ "_chr1.bed", sep='\t',
                                                header=False, index=False)

# ---------------------------------------------------------------------------------------------------
# GO_Enrichment - translate the Ensembl IDs
# ---------------------------------------------------------------------------------------------------
import GTF_Processing
annotation = '/Users/dennis/Desktop/DataEssentials/Annotations/gencode.v38.annotation.gtf.gz'
gene_sets = {"CAD-van der Harst": set(open('ExampleData/VanderHarst2018.txt').read().strip().split('\n')),
             'CAD-Schnitzler ': set(open('ExampleData/Schnitzler2024.txt').read().strip().split('\n'))}
for g_set in gene_sets:
    mapped_ids, missed_ids = GTF_Processing.match_gene_identifiers(gene_sets[g_set], gtf_file=annotation, species='human',
                                                               fields="symbol")
    out_name = 'ExampleData/VanderHarst2018_names.txt' if'Harst' in g_set else 'Exampledata/Schnitzler2024_names.txt'
    open(out_name, 'w').write('\n'.join([val['symbol'] for val in mapped_ids.values()]))


