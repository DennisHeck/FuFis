import pandas as pd
pd.options.display.max_columns = None
pd.options.display.max_rows = None

"""Gather the code that creates plots and output to be shown in the documentation. The code blocks also serve
as template for the code in the documentation itself. This is why the imports are above the clode blocks directly."""

# _________________________________________________________________________________________________________
# GTF_Processing
# _________________________________________________________________________________________________________

# ***GTF_Processing.base_code
import src.GTF_Processing as GTF_Processing
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
import src.Bed_Analysis as Bed_Analysis
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


# _________________________________________________________________________________________________________
# GOEnrichment
# _________________________________________________________________________________________________________
# ***GOEnrichment.go_enrichment1
import src.GOEnrichment as GOEnrichment
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
import src.GenomeLifter as GenomeLifter
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
import src.BigWig_Counter as BigWig_Counter
# Take a mini bed-file and get the signal from two chr21 bigwig files.
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
import src.BasicPlotter as BasicPlotter
import pandas as pd
import seaborn as sns
out_dir = 'docs/gallery/'
penguin_df = sns.load_dataset('penguins')   # Example data from seaborn.
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


# Barplot with one bar per group.
# ***BasicPlotter.basic_bars
avg_flipper_length = pd.DataFrame(penguin_df.groupby('species')['flipper_length_mm'].mean())
BasicPlotter.basic_bars(penguin_df, x_col='species', y_col='flipper_length_mm', formats=['png'],
                        x_order=['Chinstrap', 'Adelie', 'Gentoo'], title='Group of pinguins on land is called a waddle',
                        output_path=out_dir, y_label='Flipper length [mm]', rotation=None, palette='glasbey_cool')
# ---

# ***BasicPlotter.basic_bars2
avg_bill_length = pd.DataFrame(penguin_df.groupby('species')['bill_length_mm'].mean())
BasicPlotter.basic_bars(avg_bill_length, x_col='species', y_col='bill_length_mm',
                        x_order=['Chinstrap', 'Adelie', 'Gentoo'], title='Group of pinguins on land is called a waddle',
                        output_path=out_dir, y_label='Bill length [mm]', rotation=None, palette='glasbey_cool')
# ---

