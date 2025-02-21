import pandas as pd
pd.options.display.max_columns = None
pd.options.display.max_rows = None
import sys
sys.path.append("src/")

"""Gather the code that creates plots and output to be shown in the documentation. The code blocks also serve
as template for the code in the documentation itself. This is why the imports are above the clode blocks directly."""

# _________________________________________________________________________________________________________
# GTF_Processing
# _________________________________________________________________________________________________________

# ***GTF_Processing.base_code
import GTF_Processing
annotation = 'ExampleData/gencode.v38.annotation_Mini.gtf'
gene_list = ['ENSG00000160294', 'ENSG00000279493', 'ENSG00000279720']
# ---

# ***GTF_Processing.gene_window_bed
promoter_regions = GTF_Processing.gene_window_bed(gtf_file=annotation, extend=200, gene_set=gene_list, tss_type='5')
print(promoter_regions)
# ---
open("docs/gallery/src.GTF_Processing.gene_window_bed.txt", 'w').write(str(promoter_regions))

# ***GTF_Processing.gene_body_bed
gene_bodies = GTF_Processing.gene_body_bed(gtf_file=annotation, gene_set=gene_list, dict_only=False)
print(gene_bodies)
# ---
open("docs/gallery/src.GTF_Processing.gene_body_bed.txt", 'w').write(str(gene_bodies))

# ***GTF_Processing.gene_feature_bed
gene_exons = GTF_Processing.gene_feature_bed(gtf_file=annotation, feature='exon', gene_set=gene_list, dict_only=False,
                                             merge=False, keep_strand=False)
# This is list has 118 entries, so let's only print the first 5.
print(''.join([str(x) for x in gene_exons[:5]]))
# ---
open("docs/gallery/src.GTF_Processing.gene_feature_bed.txt", 'w').write(''.join([str(x) for x in gene_exons[:5]]))

# ***GTF_Processing.gene_introns_bed
# Getting the introns is a bit more tricky, as they are not explicitly annotated. We get them by subtracting all
# other annotations from the gene bodies.
gene_introns = GTF_Processing.gene_introns_bed(gtf_file=annotation, gene_set=gene_list)
print(''.join([str(x) for x in gene_introns[:5]]))
# ---
open("docs/gallery/src.GTF_Processing.gene_introns_bed.txt", 'w').write(''.join([str(x) for x in gene_introns[:5]]))

# ***GTF_Processing.match_gene_identifiers.ensembl
# The best task of all, map different identifiers. Let's start with Ensembl IDs to gene symbols, which can often
# be successfully done by a gtf-file alone.
gene_ids = ['ENSG00000160294', 'ENSG00000279493', 'ENSG00000279720']
mapped_ids, missed_ids = GTF_Processing.match_gene_identifiers(gene_ids, gtf_file=annotation, species='human',
                                                               fields="symbol")
print(mapped_ids)
# ---
open("docs/gallery/src.GTF_Processing.match_gene_identifiers.ensembl.txt", 'w').write(str(mapped_ids))

# ***GTF_Processing.match_gene_identifiers.symbols
# Next we start from symbols, which is more prone to failing. We can also query for the Entrez ID.
gene_symbols = ['AXL', 'MYO10', 'ATP5SL']
mapped_symbols, missed_symbols = GTF_Processing.match_gene_identifiers(gene_symbols, gtf_file=annotation, species='human',
                                                                       fields="ensembl,entrezgene")
print(mapped_symbols)
# ---
open("docs/gallery/src.GTF_Processing.match_gene_identifiers.symbols.txt", 'w').write(str(mapped_symbols))

# ***GTF_Processing.match_gene_identifiers.entrez
# If we want to lookup Entrez IDs, we have to add entrezgene to the scope in which mygene looks.
gene_entrez = [4651, 558]
mapped_entrez, missed_entrez = GTF_Processing.match_gene_identifiers(gene_entrez, gtf_file=annotation, species='human',
                                                                     scopes='symbol,entrezgene', fields="ensembl")
print(mapped_entrez)
# ---
open("docs/gallery/src.GTF_Processing.match_gene_identifiers.entrez.txt", 'w').write(str(mapped_entrez))


# ***GTF_Processing.bed_to_length_dict
# Convenient, for example, to get the total exon length.
gene_exons = GTF_Processing.gene_feature_bed(gtf_file=annotation, feature='exon', gene_set=gene_list, dict_only=False,
                                             merge=False, keep_strand=False)
exon_lengths = GTF_Processing.bed_to_length_dict(gene_exons)
print(exon_lengths)
# ---
open("docs/gallery/src.GTF_Processing.bed_to_length_dict.txt", 'w').write(str(exon_lengths))

# ***GTF_Processing.bed_to_feature_dict
# Admittedly, has quite specific use cases. It can be handy to collect all locations on gene level or on the level of
# any other identifier.
gene_exons = GTF_Processing.gene_feature_bed(gtf_file=annotation, feature='exon', gene_set=gene_list, dict_only=False,
                                             merge=False, keep_strand=False)
exon_dict = GTF_Processing.bed_to_feature_dict(gene_exons)
print(exon_dict['ENSG00000279493'])
# ---
open("docs/gallery/src.GTF_Processing.bed_to_feature_dict.txt", 'w').write(str(exon_dict['ENSG00000279493']))


# ***GTF_Processing.gene_biotypes
# Gets the biotypes as annotated in a gtf-file, and also provides a high-level annotation.
biotype_dict = GTF_Processing.gene_biotypes(gtf_file=annotation, gene_set=gene_list)
print(biotype_dict)
# ---
open("docs/gallery/src.GTF_Processing.gene_biotypes.txt", 'w').write(str(biotype_dict))


# ***GTF_Processing.gene_feature_table
# Gets the biotypes as annotated in a gtf-file, and also provides a high-level annotation.
feature_table = GTF_Processing.gene_feature_table(gtf_file=annotation, gene_set=gene_list)
print(feature_table)
# ---
open("docs/gallery/src.GTF_Processing.gene_feature_table.txt", 'w').write(str(feature_table))


# _________________________________________________________________________________________________________
# Bed_Analysis
# _________________________________________________________________________________________________________
# ***Bed_Analysis.base_code
import Bed_Analysis
from pybedtools import BedTool
out_dir = 'docs/gallery/'  # Replace with wherever you want to store it.
# ---

# ***Bed_Analysis.gene_location_bpwise
example_bed_file = "ExampleData/H3K27acPeaks_chr21.narrowPeak"
bed_dict = {'Example peaks': example_bed_file}  # Can have multiple entries, producing output for each.
annotation = 'ExampleData/gencode.v38.annotation_chr21Genes.gtf'
# Use a small set of peaks as example.
region_locs, total_locs = Bed_Analysis.gene_location_bpwise(bed_dict=bed_dict, gtf_file=annotation,
                                                            plot_path=out_dir, tss_type='5', palette='glasbey_cool', formats=['png'])
# We can add overlap with an additional bed-file as pie piece. In this example we use a few peaks from the original
# bed-file itself.
external_bed = 'ExampleData/H3K27acPeaks_chr21_subset.narrowPeak.txt'
region_locs, total_locs = Bed_Analysis.gene_location_bpwise(bed_dict=bed_dict, gtf_file=annotation, external_bed={"External bed": external_bed},
                                                            plot_path=out_dir+"InclExternal", tss_type='5', palette='glasbey_cool', formats=['png'])
# ---

# ***Bed_Analysis.intersection_heatmap
# For the example, let's create three BedTool object that match the image above. The first one with three large regions,
# the second repeating two of those regions, and the last having multiple small regions inside one of those.
large_regions = BedTool('\n'.join(['chr1\t1\t1000', 'chr1\t2000\t3000', 'chr1\t4000\t5000']), from_string=True)
subset_regions = BedTool('\n'.join(['chr1\t1\t1000', 'chr1\t2000\t3000']), from_string=True)
small_regions = BedTool('\n'.join(['chr1\t100\t300', 'chr1\t400\t600', 'chr1\t700\t900']), from_string=True)
multi_bed_dict = {'Large regions': large_regions,
                  'Subset regions': subset_regions,
                  'Small regions': small_regions}
Bed_Analysis.intersection_heatmap(multi_bed_dict, region_label='peaks', plot_path=out_dir, annot_nums=True,  x_size=10, y_size=7,
                           wspace=1.3, hspace=0.7, width_ratios=[0.05, 0.05, 0.96], height_ratios=[0.05, 0.97], formats=['png'])
# ---

# ***Bed_Analysis.intersection_heatmap_fisher
# Now we repeat the analysis above, but instead of just doing the pairwise-intersection, we make use of the fisher mode
# of the function. For that, we need data with more entries (limited to chr1 though). We use the scATAC-seq data from Hocker 2021 (10.1126/sciadv.abf1444)
# which was aggregated on cell type-level. They defined a union set of ATAC peaks and calculated the RPKM per cell type cluster.
# For four cell types (vCM: ventricular cardiomyocytes, MAC: macrophages, LC: lymphocytes) we subset this union to those peaks with RPKM ≥ 2 in
# the respective cell type. These peaks we then intersect with GWAS SNPs from the GWAS catalog for the traits heart rate
# (EFO_0004326) and immune system disease (EFO_0000540). Since we need a group for comparison for the Fisher's exact test,
# we take all chr1 ATAC peaks as background (the function will remove the foreground peaks from the background).
fisher_bed_dict = {"vCM": 'ExampleData/MultiBedFisher/Hocker2021_ATAC_vCM_chr1.bed',
                   "MAC": 'ExampleData/MultiBedFisher/Hocker2021_ATAC_MAC_chr1.bed',
                   "LC": 'ExampleData/MultiBedFisher/Hocker2021_ATAC_LC_chr1.bed',
                   'GWAS\nHeartRate': 'ExampleData/MultiBedFisher/GWAS_EFO_0004326_heart-rate.txt',
                   'GWAS\nHeartDisease': 'ExampleData/MultiBedFisher/GWAS_EFO_0003777_heart-disease.txt'}
fisher_background = {c: 'ExampleData/MultiBedFisher/Hocker2021_ATAC_Union_chr1.bed' for c in ['vCM', 'MAC', 'LC']}
Bed_Analysis.intersection_heatmap(fisher_bed_dict, region_label='peaks', fisher=True, fisher_background=fisher_background,
                                  row_beds=['vCM', 'MAC', 'LC'], col_beds=['GWAS\nHeartRate', 'GWAS\nHeartDisease'],
                                  plot_path=out_dir+"HockerATAC", annot_nums=True,  x_size=10, y_size=8, n_cores=2,
                                  wspace=0.4, hspace=0.9, width_ratios=[0.05, 0.2, 0.96], height_ratios=[0.08, 0.97], formats=['png'])
# ---

# ***Bed_Analysis.upset_to_reference
# Similar to the asymmetric above, this UpSet plot is meant for checking the overlap of bed-files with different sizes.
# However, this one shows the overlap only with respect to one selected bed-file. Let's use the same example bed objects
# as for the previous functions, and show the overlap only with respect to the second set of regions. One of the two
# peaks overlaps with both the large regions and the small regions, and one only with the large regions.
large_regions = BedTool('\n'.join(['chr1\t1\t1000', 'chr1\t2000\t3000', 'chr1\t4000\t5000']), from_string=True)
subset_regions = BedTool('\n'.join(['chr1\t1\t1000', 'chr1\t2000\t3000']), from_string=True)
small_regions = BedTool('\n'.join(['chr1\t100\t300', 'chr1\t400\t600', 'chr1\t700\t900']), from_string=True)
multi_bed_dict = {'Large regions': large_regions,
                  'Subset regions': subset_regions,
                  'Small regions': small_regions}
Bed_Analysis.upset_to_reference(bed_files=multi_bed_dict, ref_tag='Subset regions', y_label='Intersecting regions',
                                plot_path=out_dir, formats=['png'])
# ---

# ***Bed_Analysis.peaks_peaks_overlap
# Let's get a mapping of two small example sets of regions where all regions from the second set are located in
# one region of the first set.
large_regions = BedTool('\n'.join(['chr1\t1\t1000', 'chr1\t2000\t3000', 'chr1\t4000\t5000']), from_string=True)
small_regions = BedTool('\n'.join(['chr1\t100\t300', 'chr1\t400\t600', 'chr1\t700\t900']), from_string=True)
peak_peak_map = Bed_Analysis.peaks_peaks_overlap(peak_file=large_regions, other_peak_file=small_regions)
print(peak_peak_map)
# ---
open("docs/gallery/src.Bed_Analysis.peaks_peaks_overlap.txt", 'w').write(str(peak_peak_map))


# ***Bed_Analysis.peaks_promoter_overlap
# Now instead of intersecting multiple bed-files we get the intersection with promoter regions of genes (±200bp around the TSS).
example_bed_file = "ExampleData/H3K27acPeaks_chr21.narrowPeak"
annotation = 'ExampleData/gencode.v38.annotation_chr21Genes.gtf'
peak_promoter_map, promoter_peak_map = Bed_Analysis.peaks_promoter_overlap(peak_file=example_bed_file, gtf_file=annotation, tss_type='5')
print(list(peak_promoter_map.items())[:2])
print(list(promoter_peak_map.items())[:2])
# ---
open("docs/gallery/src.Bed_Analysis.peaks_promoter_overlap.txt", 'w').write(str(list(peak_promoter_map.items())[:2]) + '\n' + str(list(promoter_peak_map.items())[:2]))

# ***Bed_Analysis.peaks_fetch_col1
# Let's assume we have a set of ATAC peaks and found differential ATAC peaks for different conditions.
base_peaks_file = 'ExampleData/BaseATAC_peaks.txt'
# We have a directory with two bed files with a header from which we want to get the average log2FC value in all the
# base peaks they overlap.
diff_atac_pattern = 'ExampleData/DiffATAC/DiffATAC_*.txt'
# ---
# ***Bed_Analysis.peaks_fetch_col2
base_peaks_log2fc, matched_ids = Bed_Analysis.peaks_fetch_col(base_regions=base_peaks_file, pattern=diff_atac_pattern, same_peaks=False, fetch_col='log2FC')
print(base_peaks_log2fc)
# ---
open("docs/gallery/src.Bed_Analysis.peaks_fetch_col.txt", 'w').write(str(base_peaks_log2fc))

# ***Bed_Analysis.promoter_fetch_col
# This one is very similar to the previous, but starts from genes to get their promoter regions to then get the
# values from other bed files in those promoters.
diff_atac_pattern = 'ExampleData/DiffATAC/DiffATAC_*.txt'
annotation = 'ExampleData/gencode.v38.annotation_chr21Genes.gtf'
promoter_log2FC, matched_ids = Bed_Analysis.promoter_fetch_col(pattern=diff_atac_pattern, gtf_file=annotation, tss_type='5',
                                                               gene_set={'ENSG0000MOCK1', 'ENSG0000MOCK2'}, fetch_col='log2FC')
print(promoter_log2FC)
# ---
open("docs/gallery/src.Bed_Analysis.promoter_fetch_col.txt", 'w').write(str(promoter_log2FC))


# ***Bed_Analysis.peaks_genebody_overlap
# Based on an example peak file and gtf-file, we get the fraction of the gene body that is covered by the peaks.
peak_file = 'ExampleData/H3K27acPeaks_chr21.narrowPeak'
annotation = 'ExampleData/gencode.v38.annotation_chr21Genes.gtf'
gene_body_fractions = Bed_Analysis.peaks_genebody_overlap(peak_file=peak_file, gtf_file=annotation,
                                                          gene_set=['ENSG00000275895', 'ENSG00000273840'])
print(gene_body_fractions)
# ---
open("docs/gallery/src.Bed_Analysis.peaks_genebody_overlap.txt", 'w').write(str(gene_body_fractions))

# ***Bed_Analysis.possible_interactions
# For a small example, get all possible region-gene combinations within 500 base pairs.
peak_file = 'ExampleData/H3K27acPeaks_chr21.narrowPeak'
annotation = 'ExampleData/gencode.v38.annotation_chr21Genes.gtf'
all_interactions = Bed_Analysis.possible_interactions(peak_file=peak_file, gtf_file=annotation, extend=500, tss_type='5')
print(list(all_interactions)[:3])
# ---
open("docs/gallery/src.Bed_Analysis.possible_interactions.txt", 'w').write(str(list(all_interactions)[:3]))


# _________________________________________________________________________________________________________
# GOEnrichment
# _________________________________________________________________________________________________________
# ***GOEnrichment.go_enrichment1
import GOEnrichment
out_dir = 'docs/gallery/'  # Replace with wherever you want to store it.
# We're using a handful of genes from studies finding genes related to Coronary Artery Disease (CAD).
# 10.1161/CIRCRESAHA.117.312086 and 10.1038/s41586-024-07022-x
gene_sets = {"CAD-van der Harst": set(open('ExampleData/VanderHarst2018.txt').read().strip().split('\n')),
             'CAD-Schnitzler ': set(open('ExampleData/Schnitzler2024.txt').read().strip().split('\n'))}

# First let's run only one of them.
go_dict = GOEnrichment.go_enrichment({"CAD-van der Harst": gene_sets['CAD-van der Harst']}, title_tag='van der Harst 2018',
                                   out_tag=out_dir+'VanDerHarst', max_terms='all', organism='hsapiens',
                                   numerate=True, wanted_sources=['GO:BP'], rotation=45, font_s=16, formats='png')
print(go_dict['CAD-van der Harst'].head(n=3))  # The print is a bit ugly with so many columsn.
# ---
open("docs/gallery/src.GOEnrichment.go_enrichment1.txt", 'w').write(str(go_dict['CAD-van der Harst'].head(n=3)))

# ***GOEnrichment.go_enrichment2
# Next, we can see how both gene sets look like. The second gene set has a huge amount of terms enriched, so we limit
# it to the top 5 from each gene set.
go_dict = GOEnrichment.go_enrichment(gene_sets, title_tag='CAD gene sets',
                                   out_tag=out_dir+'BothSets', max_terms=5, organism='hsapiens',
                                   numerate=True, wanted_sources=['GO:BP'], rotation=45, font_s=16, formats='png')
# ---

# ***GOEnrichment.go_enrichment3
# Instead of limiting the terms to the most enriched ones, an alternative is to filter for specific keywords.
go_dict = GOEnrichment.go_enrichment(gene_sets, title_tag='CAD gene sets', keywords={"GO:BP": ['angio', 'vasc', 'circ', 'muscle']},
                                   out_tag=out_dir+'BothSetsFiltered', max_terms='all', organism='hsapiens',
                                   numerate=True, wanted_sources=['GO:BP'], rotation=45, font_s=16, formats='png')
# ---


# _________________________________________________________________________________________________________
# GenomeLifter
# _________________________________________________________________________________________________________
# ***GenomeLifter.genome_lifter
import GenomeLifter
from pybedtools import BedTool

# Lift an example list of hg38 regions to hg19.
hg38_bed_file = "ExampleData/H3K27acPeaks_chr21.narrowPeak"
hg38_regions = BedTool(hg38_bed_file)
print('hg38 coordinates:')
print(''.join([str(x) for x in hg38_regions[:3]]))
hg19_regions, unliftable = GenomeLifter.genome_lifter(hg38_regions, input_version='hg38', output_version='hg19')
print('hg19 coordinates:')  # Note the output is now a list.
print(hg19_regions[:3])
# ---
open("docs/gallery/src.GenomeLifter.genome_lifter_hg38.txt", 'w').write('hg38 coordinates:\n'+''.join([str(x) for x in hg38_regions[:3]]))
open("docs/gallery/src.GenomeLifter.genome_lifter_hg19.txt", 'w').write('hg19 coordinates:\n'+str(hg19_regions[:3]))


# _________________________________________________________________________________________________________
# FIMO TFBS
# _________________________________________________________________________________________________________
# ***FIMO_TFBS_inRegions.cmd
import subprocess
import pandas as pd

# We use subprocess here to run the bash command for easier documentation. It's also possible to run it directly
# in the terminal.
bed_file = 'ExampleData/chr21_MockPeaks.bed'
PWMs = 'ExampleData/Jaspar_Hocomoco_Kellis_human_meme.txt'
fasta = 'ExampleData/chr21_MiniMock.fa'  # For the example it's a random sequence from chr21.
fimo_src = "fimo"  # If it's on the PATH, otherwise full path to the executable.
out_dir = 'docs/gallery/FIMO_inRegions'

fimo_cmd = 'python3 src/FIMO_TFBS_inRegions.py --bed_file {} --PWMs {} --fasta {} --fimo_src {} --out_dir {} --write_sequence False'.format(bed_file, PWMs, fasta, fimo_src, out_dir)
print(fimo_cmd)
subprocess.call(fimo_cmd, shell=True)
# ---
open("docs/gallery/src.FIMO_TFBS_inRegions.cmd.txt", 'w').write(fimo_cmd)

# ***FIMO_TFBS_inRegions.matrix
# The output matrix will be stored in the following path.
fimo_matrix_file = 'docs/gallery/FIMO_inRegions/Out_Fimo_TFBSMatrix.txt.gz'
fimo_matrix = pd.read_table(fimo_matrix_file, sep='\t', header=0, index_col=0)
# With the tiny example most TF have no binding site, so let's pick a few that have some.
print(fimo_matrix[['FOXH1', 'STAT1', 'STAT3']])
# ---
open("docs/gallery/src.FIMO_TFBS_inRegions.matrix.txt", 'w').write(str(fimo_matrix[['FOXH1', 'STAT1', 'STAT3']]))

# ***FIMO_TFBS_inPromoter
import subprocess
import pandas as pd

# The run is similar to the previous, but we use a gtf-file instead of a bed-file.
gtf_file = 'ExampleData/gencode.v38.annotation_chr21Genes.gtf'  # With tow mock genes covered by the mock fasta.
PWMs = 'ExampleData/Jaspar_Hocomoco_Kellis_human_meme.txt'
fasta = 'ExampleData/chr21_MiniMock.fa'  # For the example it's a random sequence from chr21.
fimo_src = "fimo"  # If it's on the PATH, otherwise full path to the executable.
out_dir = 'docs/gallery/FIMO_inPromoter'

fimo_cmd = 'python3 src/FIMO_TFBS_inPromoter.py --gtf {} --PWMs {} --fasta {} --fimo_src {} --out_dir {} --write_sequence False'.format(gtf_file, PWMs, fasta, fimo_src, out_dir)
subprocess.call(fimo_cmd, shell=True)

fimo_matrix_file = 'docs/gallery/FIMO_inPromoter/Out_Fimo_TFBSMatrix.txt.gz'
fimo_matrix = pd.read_table(fimo_matrix_file, sep='\t', header=0, index_col=0)
# With the tiny example most TF have no binding site, so let's pick a few that have some.
print(fimo_matrix[['FOXH1', 'REST', 'PLAG1']])
# ---
open("docs/gallery/src.FIMO_TFBS_inPromoter.txt", 'w').write(str(fimo_matrix[['FOXH1', 'REST', 'PLAG1']]))


# _________________________________________________________________________________________________________
# Bigwig_Counter
# _________________________________________________________________________________________________________
# ***Bigwig_Counter
import BigWig_Counter
# Take a mini bed-file and get the signal from two chr21 bigwig files from IHEC (https://ihec-epigenomes.org/epiatlas/data/).
bed_file = "ExampleData/H3K27acPeaks_chr21.narrowPeak"
bigwigs = ['ExampleData/IHECRE00000013_chr21.bigwig', 'ExampleData/IHECRE00000017_chr21.bigwig']
bed_counts, errors = BigWig_Counter.bigwig_counts(bed_file, bigwigs, n_cores=1)
print(bed_counts.head())
# ---
open("docs/gallery/src.BigWig_Counter.txt", 'w').write(str(bed_counts.head()))


# _________________________________________________________________________________________________________
# BasicPlotter
# _________________________________________________________________________________________________________
# ***BasicPlotter.base_code
# Block that has to be executed for all.
import BasicPlotter
import pandas as pd
import seaborn as sns
out_dir = 'docs/gallery/'
penguin_df = sns.load_dataset('penguins')   # Example data from seaborn.
print(penguin_df.head())
# ---
open("docs/gallery/src.BasicPlotter_base_code.txt", 'w').write(str(penguin_df.head()))


# Barplot with one bar per group.
# ***BasicPlotter.basic_bars
# Make a basic barplot with the average flipper length per species.
avg_flipper_length = pd.DataFrame(penguin_df.groupby('species')['flipper_length_mm'].mean())
BasicPlotter.basic_bars(avg_flipper_length, x_col='species', y_col='flipper_length_mm', formats=['png'],
                        x_order=['Chinstrap', 'Adelie', 'Gentoo'], title='Example bar plot',
                        output_path=out_dir, y_label='Flipper length [mm]', rotation=None, palette='glasbey_cool')

# Make a second version where we additionally split by sex per species.
avg_flipper_length_sex = pd.DataFrame(penguin_df.groupby(['species', 'sex'])['flipper_length_mm'].mean()).reset_index()
BasicPlotter.basic_bars(avg_flipper_length_sex, x_col='species', y_col='flipper_length_mm', formats=['png'],
                        hue_col='sex', x_order=['Chinstrap', 'Adelie', 'Gentoo'], title='Example bar plot with hue',
                        output_path=out_dir + "SexHue", y_label='Flipper length [mm]', rotation=None,
                        palette='glasbey_cool')
# ---

# ***BasicPlotter.stacked_bars
# For multiple groups per bar, let's look at from which island the penguins came. Plot it once as absolute numbers
# and again as fraction.
species_island = penguin_df.groupby(['species', 'island']).size().reset_index().rename(columns={0: 'count'}).pivot(index='species', columns='island', values='count')
for do_fraction in [True, False]:
    BasicPlotter.stacked_bars(species_island, x_col='species', y_cols=species_island.columns, y_label='count', sort_stacks=False,
                              title='Stacked bars '+('fraction' if do_fraction else 'absolute'), output_path=out_dir, legend_out=1.33,
                              rotation=0, palette='glasbey_cool', fraction=do_fraction, formats=['png'])
# ---
# ***BasicPlotter.basic_pie
# An alternative to stacked barplots when focusing on one group are pie charts.
species_island_adelie = pd.DataFrame(penguin_df.groupby(['species', 'island']).size().reset_index().rename(columns={0: 'count'}).pivot(index='species', columns='island', values='count').loc['Adelie'])
BasicPlotter.basic_pie(species_island_adelie, title='Adelie islands', palette='glasbey_cool', numerate=True,
                       output_path=out_dir, legend_title='island', formats=['png'])
# ---

# ***BasicPlotter.basic_hist
# Look at the whole distribution of flipper length by using a histogram and split by species.
BasicPlotter.basic_hist(penguin_df, x_col='flipper_length_mm', hue_col='species', bin_num=20, title='Flipper length per species',
                        output_path=out_dir, stat='percent', palette='glasbey_cool', element='step', alpha=0.7, formats=['png'])
# ---

# ***BasicPlotter.basic_violin
# An alternative is to separate the species along the x-axis and do violin plots.
BasicPlotter.basic_violin(penguin_df, y_col='flipper_length_mm', x_col='species', title='Flipper length per species',
                          output_path=out_dir, numerate=True, palette='glasbey_cool', formats=['png'])
# Alternatively do boxplots instead and also add the individual data points as jitter.
BasicPlotter.basic_violin(penguin_df, y_col='flipper_length_mm', x_col='species', title='Flipper length per species',
                          output_path=out_dir+"BoxplotJitter", numerate=True, palette='glasbey_cool', formats=['png'],
                          boxplot=True, jitter=True)
# ---

# ***BasicPlotter.basic_2Dhist
# Plot the distribution of two features as 2D histogram. This example has very few points, so a scatter would work better.
BasicPlotter.basic_2Dhist(penguin_df, columns=['flipper_length_mm', 'body_mass_g'], bin_num=20, title='Flipper length vs body mass',
                          output_path=out_dir, cbar=True, formats=['png'])
# ---

# ***BasicPlotter.multi_mod_plot
# This one is a scatterplot with a lot of additional options.
# Start with a scatterplot where we colour the dots by the species, each dot being one penguin.
BasicPlotter.multi_mod_plot(penguin_df, score_cols=['flipper_length_mm', 'body_mass_g'], colour_col='species',
                            output_path=out_dir, title='#1: Flipper length vs body mass', alpha=1, palette='glasbey_cool',
                            msize=25, formats=['png'])
# Next, let's add markers to show the island where the penguin was measured.
BasicPlotter.multi_mod_plot(penguin_df, score_cols=['flipper_length_mm', 'body_mass_g'], colour_col='species',
                            marker_col='island', output_path=out_dir, title='#2 Flipper length vs body mass',
                            alpha=1, palette='glasbey_cool', msize=25, formats=['png'])
# Alternatively to having the colours categorical, we can also add a third continuous feature as colour.
BasicPlotter.multi_mod_plot(penguin_df, score_cols=['flipper_length_mm', 'body_mass_g'], colour_col='bill_length_mm',
                            marker_col='island', output_path=out_dir, title='#3 Flipper length vs body mass',
                            alpha=1, msize=25, formats=['png'])
# In addition, label the five heaviest for which we need an additional boolean column saying
# whether a dot should be labelled. The text for the label could be anything, here we write the body mass itself.
top5_index = penguin_df.sort_values('body_mass_g', ascending=False).index[:5]
penguin_df['add_label'] = [i in top5_index for i in penguin_df.index]
BasicPlotter.multi_mod_plot(penguin_df, score_cols=['flipper_length_mm', 'body_mass_g'], colour_col='bill_length_mm',
                            marker_col='island', label_dots=['add_label', 'body_mass_g'], output_path=out_dir+"doLabel",
                            title='#4 Flipper length vs body mass',
                            alpha=1, msize=25, formats=['png'])
# ---


# ***BasicPlotter.basic_venn
# Plot the overlap of lists of ingredients (incomplete) as a Venn diagram.
ingredients = {"Cookies": {'butter', 'sugar', 'flour', 'baking powder', 'chocolate'},
               'Apple pie': {'butter', 'sugar', 'flour', 'baking powder', 'apples'},
               'Bread': {'flour', 'yeast', 'oil', 'salt'}}
BasicPlotter.basic_venn(input_sets=ingredients, plot_path=out_dir+"Ingredients", formats=['png'])
# ---

# ***BasicPlotter.overlap_heatmap
# Do the overlap again but now once with the Jaccard index and once as fraction.
ingredients = {"Cookies": {'butter', 'sugar', 'flour', 'baking powder', 'chocolate'},
               'Apple pie': {'butter', 'sugar', 'flour', 'baking powder', 'apples'},
               'Bread': {'flour', 'yeast', 'oil', 'salt'}}
for mode in ['Fraction', 'Jaccard']:
    BasicPlotter.overlap_heatmap(inter_sets=ingredients, title="Ingredients overlap as "+mode, plot_path=out_dir+"Ingredients_"+mode,
                                 xsize=10, ysize=6, mode=mode, annot_type='Jaccard' if mode == 'Jaccard' else 'Absolute', formats='png')
# ---

# ***BasicPlotter.upset_plotter
# Use the same sets from the previous function, but now visualize the overlap with an UpsetPlot.
ingredients = {"Cookies": {'butter', 'sugar', 'flour', 'baking powder', 'chocolate'},
               'Apple pie': {'butter', 'sugar', 'flour', 'baking powder', 'apples'},
               'Bread': {'flour', 'yeast', 'oil', 'salt'}}
BasicPlotter.upset_plotter(inter_sets=ingredients, sort_categories_by='input', title_tag='Ingredients overlap',
                           plot_path=out_dir+"Ingredients", intersection_plot_elements=4, formats=['png'])
# ---


# ***BasicPlotter.cumulative_plot
# This is most instructive with diverging data e.g. logFC from RNA-seq. We use the RNA data from a study on the histone
# mark H379me2 (10.1038/s41467-022-35070-2) and a formatted file from their supplements.
rna_file = 'ExampleData/41467_2022_35070_MOESM4_ESM_E16Sub.txt'
rna_table = pd.read_table(rna_file, sep='\t', header=0)
# Add a column that groups genes into rough bins of how much of the gene body is covered.
rna_table['binned H3K79me2 GB Coverage'] = pd.cut(rna_table['H3K79me2 GB Coverage'], bins=2).astype('string')
BasicPlotter.cumulative_plot(rna_table, x_col='logFC', hue_col='binned H3K79me2 GB Coverage', palette='glasbey_cool', xlimit=[-1.5, 2],
                             add_all=True, output_path=out_dir, numerate=True, title=None, vertical_line=0, table_width=0.4, table_x_pos=1.2, formats=['png'])
# ---
# _________________________________________________________________________________________________________
# Various
# _________________________________________________________________________________________________________
# ***Various.fn_patternmatch
import Various
my_files = Various.fn_patternmatch('ExampleData/BirdCollection/Bird_*.txt')
print(my_files)
# ---
open("docs/gallery/src.Various.fn_patternmatch.txt", 'w').write(str(my_files))


# _________________________________________________________________________________________________________
# UniProtAPI
# _________________________________________________________________________________________________________
# ***UniProt_API.uniprot_domains
import UniProt_API
# Look up the annotated domains and regions of three examples.
protein_domains, protein_regions, missed_proteins, failed_requests = UniProt_API.uniprot_domains(protein_names=['KDM6A', 'DNMT3A', 'STAT2'], species='human', n_cores=1)
print(protein_domains)
print(protein_regions)
# ---
open("docs/gallery/src.UniProtAPI.uniprot_domains.txt", 'w').write(str(protein_domains) + '\n\n' + str(protein_regions))


# _________________________________________________________________________________________________________
# CoveragePlots
# _________________________________________________________________________________________________________
# ***CoveragePlots.plotHeatmap
from pybedtools import BedTool
import CoveragePlots
out_dir = 'docs/gallery/'
# Let's plot the signal of two bigwig files from IHEC (https://ihec-epigenomes.org/epiatlas/data/) in a small set of peaks and compare that to the signal in their shuffled locations.
peaks = BedTool("ExampleData/H3K27acPeaks_chr21.narrowPeak")
shuffled_peaks = peaks.shuffle(genome='hg38', chrom=True, seed=12)
bigwigs = ['ExampleData/IHECRE00000013_chr21.bigwig', 'ExampleData/IHECRE00000017_chr21.bigwig']
CoveragePlots.plotHeatmap(beds_to_plot=[peaks, shuffled_peaks], bed_labels=['Original', 'Shuffled'], bigwigs=bigwigs, bw_labels=['Sample1', 'Sample2'],
                          out_dir=out_dir, out_tag='ExampleCoveragePlot', mode='scale', perGroup=True, title='',
                          scaled_size=500, start_label='Peak start', end_label='Peak end')
# ---


# _________________________________________________________________________________________________________
# Heatmaps
# _________________________________________________________________________________________________________
# ***Heatmaps.base_code
# Block that has to be executed for all.
import Heatmaps
import seaborn as sns
out_dir = 'docs/gallery/'
penguin_df = sns.load_dataset('penguins')   # Example data from seaborn.
# ---

# ***Heatmaps.heatmap_cols
# Create a heatmap with several blocks that allow separate metrics to be shown side-by-side.
# For an example pick a few random rows and give them names so we can recognize them in the heatmaps.
sub_penguin_df = penguin_df.sample(5, random_state=12)
sub_penguin_df['Name'] = ['Pesto', 'Pickles', 'Popcorn', 'Pretzel', 'Pudding']
# Now we need to define which heatmap-blocks we want to show.
cmap_cols = {0: {'cols': ['bill_length_mm', 'bill_depth_mm'],
                 'cmap': 'mako',
                 'cbar_label': 'mm'},
             1: {'cols': ["body_mass_g"],
                 'cmap': 'viridis',  # Just for the sake of a different colour.
                 'cbar_label': 'g'}
             }
# We can additionally choose to add labels into cells, can be any other column.
annot_cols = {"body_mass_g": 'island'}
Heatmaps.heatmap_cols(sub_penguin_df, cmap_cols=cmap_cols, plot_out=out_dir+"SubPenguins", row_label_col='Name', class_col='species',
                 x_size=14, y_size=5, annot_cols=annot_cols, width_ratios=None, wspace=0.8,
                 annot_s=10, ticksize=14, heat_ticksize=14, square=False, x_rotation=0, formats=['png'])
# ---

# ***Heatmaps.clustermap
# Let's cluster some penguins based on their measured features, after removing entries with NAs.
non_na_penguins = penguin_df[~penguin_df.isna().values.any(axis=1)]
# Once with the original values and once with z-scoring the columns (flag takes None or the axis).
# The row names will be the index numbers, which are not meaningful here, but suffices for illustration.
for do_z in [None, 0]:
    cmap = 'mako' if do_z is None else 'bwr'
    centre = None if do_z is None else 0
    Heatmaps.clustermap(non_na_penguins, columns=['bill_length_mm', 'bill_depth_mm', 'flipper_length_mm'],
                        row_column='index', z_score=do_z, centre=centre, cbar_label='mm', class_col='species',
                        title="Clustered penguins", plot_out=out_dir+"Penguins_"+str(do_z)+'ZScore', cmap=cmap,
                        x_size=8, y_size=10, row_cluster=True, col_cluster=True, formats=['png'])
# ---


# _________________________________________________________________________________________________________
# MEME_Formatting
# _________________________________________________________________________________________________________
# ***MEME_Formatting.base_code
# Block that has to be executed for all.
import MEME_Formatting
out_dir = 'docs/gallery/'
meme_file = 'ExampleData/Jaspar_Hocomoco_Kellis_human_meme.txt'
annotation = 'ExampleData/gencode.v38.annotation_chr21Genes.gtf'  # It's not the full one, so fewer hits expected.
# ---

# ***MEME_Formatting.meme_id_map
# Get the Ensembl ID for the TFs in our motif meme-file.
tf_ids, all_tf_names, misses = MEME_Formatting.meme_id_map(meme_file=meme_file, gtf_file=annotation, species='human')
print('TBXT', tf_ids['TBXT'])
print('MAX::MYC', tf_ids['MAX::MYC'])
# ---
open("docs/gallery/src.MEME_Formatting.meme_id_map.txt", 'w').write('TBXT' + '\t' + str(tf_ids['TBXT']) + '\n' + 'MAX::MYC' + '\t' + str(tf_ids['MAX::MYC']))

# ***MEME_Formatting.meme_monomer_map
# Useful when in need of the individual monomers or removal of the motif versions.
tf_monomer_map, all_monomer_names = MEME_Formatting.meme_monomer_map(meme_file=meme_file)
print('BHLHA15(MA0607.2)', tf_monomer_map['BHLHA15(MA0607.2)'])
print('MAX::MYC', tf_monomer_map['MAX::MYC'])
print('all monomers', all_monomer_names[:4])
# ---
open("docs/gallery/src.MEME_Formatting.meme_monomer_map.txt", 'w').write('BHLHA15(MA0607.2)' + '\t' + str(tf_ids['BHLHA15(MA0607.2)']) +
                                                                         '\n' + 'MAX::MYC' + '\t' + str(tf_ids['MAX::MYC']) +
                                                                         '\n' + 'all monomers' + '\t' + str(all_monomer_names[:4]))

# ***MEME_Formatting.subset_meme
# Subset a meme-file, which is useful for example for excluding TFs that are not expressed.
MEME_Formatting.subset_meme(meme_file, motif_names=['MAX::MYC', 'TBXT'], out_file=out_dir+"Subset_meme.txt",
                            include_dimers=True, exact_match=False)
print(open(out_dir+"Subset_meme.txt").read())
# ---

# _________________________________________________________________________________________________________
# MakeTitlesGreatAgain
# _________________________________________________________________________________________________________
# ***MakeTitlesGreatAgain
# Run it with a mini example with one bib entry, here from Python, but you can also run the command
# that's sent to subprocess in your bash.
import subprocess
subprocess.call("python3 src/MakeTitlesGreatAgain.py ExampleData/ExampleBib.bib > ExampleData/ExampleBib_great.bib", shell=True)
print(open('ExampleData/ExampleBib.bib').readlines()[2])
print(open('ExampleData/ExampleBib_great.bib').readlines()[2])
# ---
open('docs/gallery/MakeTitlesGreatAgain_1.txt', 'w').write(open('ExampleData/ExampleBib.bib').readlines()[2])
open('docs/gallery/MakeTitlesGreatAgain_2.txt', 'w').write(open('ExampleData/ExampleBib_great.bib').readlines()[2])

