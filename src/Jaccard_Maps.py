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


def jc_heatmap(inter_sets, title="", plot_path='', xsize=12, ysize=8, annot=True, mode='JC', annot_type='JC',
               annot_size=13, matrix_only=False):
    """
    Plot a symmetric heatmap of the Jaccard index of the sets given in inter_sets. Combinations are formed by the keys
    of the dictionary.
    @param mode: 'JC' to show the Jaccard index of the intersection, 'Fraction' to show the percentage.
    @param annot_type: 'JC' to get the JC index written, 'abs' to get the absolute number of shared items.
    @param matrix_only: Whether only the matrix should be given without creating a plot, will return two dfs, one with
    the shared metric and the other the cell labels.
    """
    if mode == 'JC':
        cbar_label = "Jaccard index"
    elif mode == 'Fraction':
        cbar_label = "Fraction shared [ (row & column) / row ]"
    else:
        print("ERROR: Invalid mode", mode)
        return
    set_tags = list(inter_sets.keys())
    combos = list(itertools.combinations(set_tags, 2))
    shared_mat = np.zeros([len(inter_sets), len(inter_sets)])
    annot_mat = np.zeros([len(inter_sets), len(inter_sets)])

    for comb in combos:
        shared = len(inter_sets[comb[0]] & inter_sets[comb[1]])
        if shared == 0:
            continue
        if mode == 'JC':
            shared_mat[set_tags.index(comb[0])][set_tags.index(comb[1])] = shared / len(inter_sets[comb[0]].union(inter_sets[comb[1]]))
            shared_mat[set_tags.index(comb[1])][set_tags.index(comb[0])] = shared / len(inter_sets[comb[0]].union(inter_sets[comb[1]]))
        elif mode == "Fraction":
            shared_mat[set_tags.index(comb[0])][set_tags.index(comb[1])] = shared / len(inter_sets[comb[0]])
            shared_mat[set_tags.index(comb[1])][set_tags.index(comb[0])] = shared / len(inter_sets[comb[1]])
        if annot_type == "JC":
            annot_mat[set_tags.index(comb[0])][set_tags.index(comb[1])] = shared / len(inter_sets[comb[0]].union(inter_sets[comb[1]]))
            annot_mat[set_tags.index(comb[1])][set_tags.index(comb[0])] = shared / len(inter_sets[comb[0]].union(inter_sets[comb[1]]))
        else:
            annot_mat[set_tags.index(comb[0])][set_tags.index(comb[1])] = len(inter_sets[comb[0]] & inter_sets[comb[1]])
            annot_mat[set_tags.index(comb[1])][set_tags.index(comb[0])] = len(inter_sets[comb[0]] & inter_sets[comb[1]])
    if annot_type == 'JC':
        annot_mat = annot_mat.round(2)
    else:
        for n in range(len(set_tags)):  # Fill the diagnoal with the absolute number of items.
            annot_mat[n][n] = len(inter_sets[set_tags[n]])
            if len(inter_sets[set_tags[n]]) > 0:
                shared_mat[n][n] = 1
        annot_mat = annot_mat.astype(int)

    if matrix_only:
        shared_df = pd.DataFrame(shared_mat, columns=set_tags, index=set_tags)
        annot_df = pd.DataFrame(annot_mat, columns=set_tags, index=set_tags)
        return shared_df, annot_df

    f, ax = plt.subplots(figsize=(xsize, ysize))
    shared_heat = sns.heatmap(shared_mat, cmap="mako", ax=ax, square=True, vmin=0, vmax=1, rasterized=True,
                              yticklabels=set_tags, xticklabels=set_tags, fmt='',
                              annot=None if not annot else annot_mat,
                              annot_kws={'size': annot_size}, cbar_kws={'label': cbar_label, 'pad': 0.01})
    shared_heat.set_xticklabels(shared_heat.get_xmajorticklabels(), fontsize=13, rotation=90, fontweight='bold')
    shared_heat.set_yticklabels(shared_heat.get_ymajorticklabels(), fontsize=13, rotation=0, fontweight='bold')
    shared_heat_cbar = shared_heat.collections[0].colorbar
    shared_heat_cbar.ax.tick_params(labelsize=14)
    shared_heat_cbar.ax.yaxis.label.set_fontsize(14)
    plt.suptitle(title, y=0.93, size=20, fontweight='bold')
    plt.subplots_adjust(wspace=0.02)
    plt.savefig((plot_path + "_SharedHeatmap.pdf").replace(" ", "_"), bbox_inches='tight')
    plt.close()
