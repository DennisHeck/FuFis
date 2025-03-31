import os
import gzip
import pandas as pd
import numpy as np
import Heatmaps
import BasicPlotter
import JaccardMaps


def harry_plotter(abc_folder, tag_order=None, plot_path='', same_enhancer=True, colnaming=True):
    """A collection of function plotting metrics for all ABCpp files in a folder.
    @param abc_folder: directory from which all gzipped files with the substring ABCpp_scoredInteractions will be taken.
    @param tag_order: Order of the prefixes of the interaction files. If not given take the order from os.listdir().
    @param same_enhancer: If True will create a Jaccard heatmap of the interactions, if not an asymmetric heatmap to
    @param colnaming: Whether the name tag for the files is at the column label position or the prefix to _ABCpp_.
    account for multiple overlaps.
    """
    if colnaming:
        abc_files = {x.split('_ABCpp_scoredInteractions_')[-1].split('.txt')[0]: abc_folder+'/'+x for x in os.listdir(abc_folder)
                     if 'ABCpp_scoredInteractions' in x and x.endswith('.gz')}
    else:
        abc_files = {x.split('_ABCpp_scoredInteractions')[0]: abc_folder + '/' + x for x in os.listdir(abc_folder)
                     if 'ABCpp_scoredInteractions' in x and x.endswith('.gz')}
    if not tag_order:
        tag_order = list(abc_files.keys())

    if len(tag_order) <= 10:
        palette = 'tab10'
    elif len(tag_order) <= 20:
        palette = 'tab20'
    else:
        palette = 'glasbey'

    # Plot the number of interactions.
    num_inter = {}
    for tag in tag_order:
        num_inter[tag] = len(gzip.open(abc_files[tag], 'rt').readlines()) - 1
    interaction_count = pd.DataFrame(num_inter.items(), columns=["Run", 'Number of interactions'])

    BasicPlotter.basic_bars(interaction_count, y_col='Run', x_col='Number of interactions', y_size=1+len(tag_order),
                            x_size=6, title="Number of interactions per run", output_path=plot_path, rotation=90,
                            palette=palette, font_s=16)

    # Get number of enhancers per gene and vice versa.
    enh_gene_counter = {t: {} for t in tag_order}
    gene_enh_counter = {t: {} for t in tag_order}
    for tag in tag_order:
        abc_head = {x: i for i, x in enumerate(gzip.open(abc_files[tag], 'rt').readline().strip().split('\t'))}
        with gzip.open(abc_files[tag], 'rt') as abc_in:
            abc_in.readline()
            for entry in abc_in:
                enh = '\t'.join([entry.strip().split('\t')[abc_head[c]] for c in ['#chr', 'start', 'end']])
                gene = entry.strip().split('\t')[abc_head['Ensembl ID']]
                if enh not in enh_gene_counter[tag]:
                    enh_gene_counter[tag][enh] = 0
                enh_gene_counter[tag][enh] += 1
                if gene not in gene_enh_counter[tag]:
                    gene_enh_counter[tag][gene] = 0
                gene_enh_counter[tag][gene] += 1
    genes_per_enh = pd.DataFrame.from_dict(enh_gene_counter, orient='columns').fillna(0)# columns=['Run', "Genes per enhancer"])
    enh_per_gene = pd.DataFrame.from_dict(gene_enh_counter, orient='columns').fillna(0)

    # Plot the average per gene/enhancer across all runs.
    genes_per_enh['Average genes per enhancer'] = genes_per_enh.mean(axis=1)
    genes_per_enh['type'] = 'Enhancer'  # Needed to work as long-format with seaborn.
    enh_per_gene['Average enhancers per gene'] = enh_per_gene.mean(axis=1)
    enh_per_gene['type'] = 'Gene'

    BasicPlotter.basic_violin(genes_per_enh, y_col='Average genes per enhancer', x_col='type', boxplot=False,
                              title="Number of genes per enhancer", output_path=plot_path, colour='#3a77e0',
                              font_s=15, numerate=True, palette=None, xsize=5, ysize=7, rotation=0)

    BasicPlotter.basic_violin(enh_per_gene, y_col='Average enhancers per gene', x_col='type', boxplot=False,
                              title="Number of enhancers per gene", output_path=plot_path, colour='#eda64a',
                              font_s=15, numerate=True, palette=None, xsize=5, ysize=7, rotation=0)

    # And Plots where each run has its own boxplot.  TODO proper y_col label
    BasicPlotter.basic_violin(genes_per_enh[tag_order], y_col=None, x_col=None, x_order=tag_order, boxplot=True,
                              title="Number of genes per enhancer", output_path=plot_path+"GenesPerEnh", font_s=15,
                              numerate=False, palette=palette, xsize=5+len(tag_order), ysize=8, rotation=90)
    BasicPlotter.basic_violin(enh_per_gene[tag_order], y_col=None, x_col=None, x_order=tag_order, boxplot=True,
                              title="Number of enhancers per gene", output_path=plot_path+"EnhPerGene", font_s=15,
                              numerate=False, palette=palette, xsize=5+len(tag_order), ysize=8, rotation=90)

    if not same_enhancer:
        Heatmaps.interaction_intersection_diff_enhancer(abc_folder=abc_folder, tag_order=tag_order,
                                                        plot_path=plot_path, x_size=16, y_size=10, annot_s=15)
    else:
        interaction_sets = {c: set() for c in abc_files.keys()}
        for tag, file in abc_files.items():
            inter_head = {x: i for i, x in enumerate(gzip.open(file, 'rt').readline().strip().split('\t'))}
            for line in gzip.open(file, 'rt').read().strip().split('\n')[1:]:
                this_gene = line.split('\t')[inter_head['Ensembl ID']].split('_')[0]
                this_peak = line.split('\t')[inter_head['PeakID']]
                interaction_sets[tag].add(this_gene + '#' + this_peak)

        BasicPlotter.overlap_heatmap(interaction_sets, title="ABC interaction overlap", plot_path=plot_path,
                                     mode='Fraction', annot_type='Absolute', xsize=12, ysize=8)

