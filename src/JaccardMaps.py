import gzip
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import itertools


def jc_mat(abc_files, plot_path, tag_order=None):
    """
    Plot a symmetric heatmap of the Jaccard index of interactions of the files in abc_files.
    Assumes that the candidate enhancers are identical.
    @param abc_files: Dictionary with the path to the ABC-files whose interactions to compare {tag: path}.
    """

    interaction_sets = {c: set() for c in abc_files.keys()}
    if not tag_order:
        inter_tags = list(interaction_sets.keys())
    else:
        inter_tags = tag_order

    for tag, file in abc_files.items():
        inter_head = {x: i for i, x in enumerate(gzip.open(file, 'rt').readline().strip().split('\t'))}
        for line in gzip.open(file, 'rt').read().strip().split('\n')[1:]:
            this_gene = line.split('\t')[inter_head['Ensembl ID']].split('_')[0]
            this_peak = line.split('\t')[inter_head['PeakID']]
            interaction_sets[tag].add(this_gene + '#' + this_peak)

    combos = list(itertools.combinations(inter_tags, 2))
    jc_mat = np.ones([len(interaction_sets), len(interaction_sets)])

    for comb in combos:
        jaccard = len(interaction_sets[comb[0]] & interaction_sets[comb[1]]) / len(interaction_sets[comb[0]].union(interaction_sets[comb[1]]))
        jc_mat[inter_tags.index(comb[0])][inter_tags.index(comb[1])] = jaccard
        jc_mat[inter_tags.index(comb[1])][inter_tags.index(comb[0])] = jaccard

    f, ax = plt.subplots(figsize=(jc_mat.shape[0]+3, jc_mat.shape[1]))
    shared_heat = sns.heatmap(jc_mat, cmap="mako", ax=ax, square=True, vmin=0, vmax=1, rasterized=True,
                              yticklabels=inter_tags, xticklabels=inter_tags, fmt='', annot=jc_mat.round(2),
                              annot_kws={'size': 14}, cbar_kws={'label': "Jaccard index", 'pad': 0.01})
    shared_heat.set_xticklabels(shared_heat.get_xmajorticklabels(), fontsize=13, rotation=90, fontweight='bold')
    shared_heat.set_yticklabels(shared_heat.get_ymajorticklabels(), fontsize=13, rotation=0, fontweight='bold')
    shared_heat_cbar = shared_heat.collections[0].colorbar
    shared_heat_cbar.ax.tick_params(labelsize=14)
    shared_heat_cbar.ax.yaxis.label.set_fontsize(14)
    plt.suptitle("Jaccard index of ABC interactions", y=0.93, size=20, fontweight='bold')
    plt.subplots_adjust(wspace=0.02)
    plt.savefig(plot_path + "_JaccardInteractions.pdf", bbox_inches='tight')
    plt.close()

