import pandas as pd
pd.options.display.max_columns = None
pd.options.display.max_rows = None

"""Gather the code that creates plots and output to be shown in the documentation. The code blocks also serve
as template for the code in the documentation itself. This is why the imports are above the clode blocks directly."""

# Block that has to be executed for all.
import src.BasicPlotter as BasicPlotter
import pandas as pd
import seaborn as sns
out_dir = 'docs/gallery/'
penguin_df = sns.load_dataset('penguins')

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

# ***Bed_Analysis.inter_heatmap
# For the example, let's create three BedTool object. The first one with three large regions, the second repeating
# two of those regions, and the last having multiple small regions inside one of those.
large_regions = BedTool('\n'.join(['chr1\t1\t1000', 'chr1\t2000\t3000', 'chr1\t4000\t5000']), from_string=True)
subset_regions = BedTool('\n'.join(['chr1\t1\t1000', 'chr1\t2000\t3000']), from_string=True)
small_regions = BedTool('\n'.join(['chr1\t1\t10', 'chr1\t11\t20', 'chr1\t21\t30']), from_string=True)
multi_bed_dict = {'Large peaks': large_regions,
                  'Subset peaks': subset_regions,
                  'Small peaks': small_regions}
Bed_Analysis.inter_heatmap(multi_bed_dict, region_label='peaks', plot_path=out_dir, annot_nums=True,  x_size=10, y_size=7,
                           wspace=1, hspace=0.6, width_ratios=[0.05, 0.05, 0.96], height_ratios=[0.05, 0.97], formats=['png'])
# ---

# ***Bed_Analysis.upset_to_reference
# Similar to the asymmetric above, this UpSet plot is meant for checking the overlap of bed-files with different sizes.
# However, this one shows the overlap only with respect to one selected bed-file. Let's use the same example bed objects
# as for the previous functions, and show the overlap only with respect to the second set of regions. One of the two
# peaks overlaps with both the large regions and the small regions, and one only with the large regions.
large_regions = BedTool('\n'.join(['chr1\t1\t1000', 'chr1\t2000\t3000', 'chr1\t4000\t5000']), from_string=True)
subset_regions = BedTool('\n'.join(['chr1\t1\t1000', 'chr1\t2000\t3000']), from_string=True)
small_regions = BedTool('\n'.join(['chr1\t1\t10', 'chr1\t11\t20', 'chr1\t21\t30']), from_string=True)
multi_bed_dict = {'Large peaks': large_regions,
                  'Subset peaks': subset_regions,
                  'Small peaks': small_regions}
Bed_Analysis.upset_to_reference(bed_files=multi_bed_dict, ref_tag='Subset peaks', y_label='Intersecting regions',
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








