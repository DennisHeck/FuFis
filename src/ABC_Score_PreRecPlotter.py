import pandas as pd
from matplotlib import pyplot as plt
import sklearn.metrics as metrics
from matplotlib import cm
import numpy as np
import math
from scipy import stats
import gzip
# import ColoursAndShapes

# CARE the following shapes are not plotted for whatever reason: ['1', '2', '3', '4', '+', 'x']
marker_shapes = ['o', 'D', 'P', 's', '*', 'X', '<', 'p', 'd', '^', 'v', '>', 'H', '$O$', '$D$', '$U$', '$Y$', '$N$']
# https://davidmathlogic.com/colorblind/
two_contrasts = [['#E1BE6A', '#40B0A6'],  # orange, teal
                 ['#FFC20A', '#0C7BDC'],  # yellow blue
                 ['#994F00', '#006CD1']]  # brown blue
# https://personal.sron.nl/~pault/
tol_bright = ['#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']
tol_vibrant = ['#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988', '#BBBBBB']

palette = tol_vibrant


def pre_rec_fetcher(binning, ele_gene_pairs, true_enhancers, score_col='ABC Score', lower_limit=0, upper_limit=1,
                    mark_thresh=None, mode='pr'):

    def fetch_prerec(thresholds):
        rec_prec_list = []
        recall_hit = []
        precision_hit = []
        for thresh in thresholds:
            true_positives = true_enhancers[pd.to_numeric(true_enhancers[score_col]) >= thresh]
            # How many of the positives did we find.
            recall = true_positives.shape[0] / true_enhancers.shape[0]
            if ele_gene_pairs.dtypes['Significant'] == bool:
                false_positives = ele_gene_pairs[(ele_gene_pairs['Significant'] == 0) & (pd.to_numeric(ele_gene_pairs[score_col]) >= thresh)]
                true_negatives = len(ele_gene_pairs[ele_gene_pairs['Significant'] == 0])
            else:
                false_positives = ele_gene_pairs[(ele_gene_pairs['Significant'] == 'False') & (pd.to_numeric(ele_gene_pairs[score_col]) >= thresh)]
                true_negatives = len(ele_gene_pairs[ele_gene_pairs['Significant'] == 'False'])
            if true_positives.shape[0] != 0:
                precision = true_positives.shape[0] / (true_positives.shape[0] + false_positives.shape[0])
            else:
                precision = 0
            if true_negatives > 0:
                fpr = len(false_positives) / true_negatives
            else:
                fpr = 0
            if mode == 'pr':
                rec_prec_list.append([recall, precision, thresh])
            else:  # ROC
                rec_prec_list.append([fpr, recall, thresh])

            # if round(recall, 2) == 0.7:
            #     recall_hit.append(thresh)
            # if round(precision, 2) == 0.1:
            #     precision_hit.append(thresh)
            # if round(thresh, 4) == 0.02:
            #     print(recall, precision)
        # print(score_col, '70% recall', round(np.mean(recall_hit), 4))
        # print(score_col, '10% precision', round(np.mean(precision_hit), 4))
        return rec_prec_list

    binned_list = fetch_prerec([lower_limit + (upper_limit-lower_limit) * x / binning for x in range(0, binning+1)])
    if mark_thresh:
        singular_pos = fetch_prerec([mark_thresh])
    else:
        singular_pos = None
    return binned_list, singular_pos


def pre_rec_plotter(plot_df, sig_col, plot_cols, steps=10000, output_path='', colours=None, no_plot=False,
                    recall_start=0, zorder=None, legend_s=18, mode='pr', thresh=None, legend_out=False,
                    lines_solid=False):
    if not colours:
        if len(plot_cols) == 2:
            colours = two_contrasts[0]
        else:
            colours = palette*10
    lines = ['solid', 'dashed', 'solid', 'dashdot']*10
    pr_dict = {c: None for c in plot_cols}
    if plot_df.dtypes[sig_col] == bool:
        true_enhancers = plot_df[plot_df[sig_col] == 1]
    else:
        true_enhancers = plot_df[plot_df[sig_col].isin(['True', 'TRUE', 'true', 'T', '1'])]

    f, ax = plt.subplots(1, figsize=(8, 8))
    plt.subplots_adjust(left=0.08, right=0.98)
    ax.set_axisbelow(True)
    ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='both')
    auc_coll = []
    auc_output = []
    if not zorder:
        zorder = [i+10 for i in range(len(plot_cols))]
    for n, this_col in enumerate(plot_cols):
        lower_lim = min(plot_df[this_col])
        upper_lim = max(plot_df[this_col])
        rec_prec_list, thresh_pos = pre_rec_fetcher(steps, plot_df, true_enhancers, score_col=this_col,
                                                    lower_limit=lower_lim, upper_limit=upper_lim, mark_thresh=thresh, mode=mode)
        pr_dict[this_col] = rec_prec_list
        sorted_by_recall = sorted(rec_prec_list, key=lambda x: x[0])
        sorted_by_recall = [x for x in sorted_by_recall if round(x[0], 3) >= recall_start]
        auc = metrics.auc([x[0] for x in sorted_by_recall], [x[1] for x in sorted_by_recall])
        auc_coll.append(this_col + ' AUC: ' + str(round(auc, 3)))
        auc_output.append([this_col, auc])
        print(this_col, round(auc, 4))

        plt.plot([x[0] for x in rec_prec_list], [x[1] for x in rec_prec_list],
                 linestyle='solid' if lines_solid else lines[n], c=colours[n],
                 label=this_col, linewidth=4, zorder=zorder[n])

        # dot_pos = sorted([x for x in rec_prec_list if round(x[0], 1) == 0.4], key=lambda x: abs(x[0]-0.4))
        # dot_label = '40% recall' if n == 0 else ""
        # plt.scatter(dot_pos[0][0], dot_pos[0][1], c=tol_bright[n], s=100, label=dot_label, edgecolors='k', zorder=20)

        if thresh_pos:
            plt.scatter(thresh_pos[0][0], thresh_pos[0][1], c=colours[n], s=100, marker='*', label='cut-off ('+str(thresh)+')', edgecolors='k', zorder=21)

        # with open(output_path + this_set['label'] + "_PreRecList.txt", 'w') as output:
        #     output.write('\t'.join(["Recall", "Precision", "Threshold"]) + '\n')
        #     for entry in rec_prec_list:
        #         output.write('\t'.join([str(e) for e in entry]) + '\n')

    if legend_out:
        ax.legend(prop={'size': legend_s, 'weight': 'normal'}, loc='upper right',
                  bbox_to_anchor=(2 if type(legend_out) == bool and not legend_out else legend_out, 1))
    else:
        ax.legend(prop={'size': legend_s, 'weight': 'normal'})
    plt.xlabel('Recall' if mode == 'pr' else 'FPR', size=22)
    plt.ylabel('Precision' if mode == 'pr' else 'TPR', size=22)
    ax.tick_params(axis='both', labelsize=18)
    plt.ylim(0, 1)
    plt.xlim(recall_start, 1)
    ax.set_title(str(len(true_enhancers))+' sig/ '+str(len(plot_df)) + ' interactions' + '\n' + '\n'.join(auc_coll), y=1.05)
    # plt.title('Performance of the ABC model\n'+str(len(plotter_list[0]['true']))+' significant / '+str(len(plotter_list[0]['df'])) + ' interactions')
    ax.set_facecolor('white')
    if not no_plot:
        f.savefig(output_path + '.pdf', bbox_inches='tight')
    plt.close()
    return auc_output, pr_dict


def pr_thresh_coloured(binning, plot_df, score_col, true_set, thresh, output_path, vline=None, hline=None):
    """Plot PR curves, but as scatter and coloured by the threshold."""
    marker_shapes = ['o', '^']
    lower_lim = plot_df[score_col].min().min()
    upper_lim = plot_df[score_col].max().max()
    cmap = cm.get_cmap('viridis')
    norm = plt.Normalize(lower_lim, upper_lim)
    f, ax = plt.subplots(1, figsize=(9, 8))
    plt.subplots_adjust(left=0.08, right=0.98)
    ax.set_axisbelow(True)
    ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='both')
    for c, col in enumerate(score_col):
        rec_prec_list, thresh_pos = pre_rec_fetcher(binning, plot_df, true_set, score_col=col, lower_limit=min(plot_df[col]),
                                                    upper_limit=max(plot_df[col]), mark_thresh=thresh)
        sorted_by_recall = sorted(rec_prec_list, key=lambda x: x[0])
        auc = metrics.auc([x[0] for x in sorted_by_recall], [x[1] for x in sorted_by_recall])
        print(col, round(auc, 3))
        plt.scatter([x[0] for x in rec_prec_list], [x[1] for x in rec_prec_list], c=[cmap(norm(x[2])) for x in rec_prec_list],
                    s=16, zorder=12, marker=marker_shapes[c], label=col)
    if vline:
        ax.axvline(vline, color='grey', linestyle='dotted')
    if hline:
        ax.axhline(hline, color='grey', linestyle='dashdot')
    cax = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.8)
    cax.set_label('score cut-off', size=18)
    cax.ax.tick_params(labelsize=18)
    plt.xlabel('Recall', size=22)
    plt.ylabel('Precision', size=22)
    ax.tick_params(axis='both', labelsize=18)
    ax.legend(markerscale=3, scatterpoints=1, fontsize=16, labelspacing=0.4, prop={'weight': 'normal', 'size': 18})

    plt.ylim(0, 1)
    plt.xlim(0, 1)

    ax.set_facecolor('white')
    f.savefig(output_path + '_PRColoured.pdf', bbox_inches='tight')


def score_effect_correlation(score_df, score_cols, effect_col, output_path, out_type='pdf'):
    """Plot the correlation of the expression change to the ABC score. Also report the Spearman correlation."""

    for c, col in enumerate(score_cols):
        score_list = [(math.log(x[0]), x[1]) for x in zip(score_df[col], score_df[effect_col]) if x[0] > 0]
        rho, pval = stats.spearmanr(a=score_list, axis=0)
        print(col, rho, pval)

        f, ax = plt.subplots(figsize=(8, 8))
        ax.set_axisbelow(True)
        ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='both')
        for sig in score_list:
            ax.scatter(sig[0], sig[1], c='#ff2a00', marker='o', alpha=1, zorder=24, s=6)
        plt.ylabel(effect_col, size=16)
        plt.xlabel('log(ABC Score)', size=16)
        ax.tick_params(axis='both', labelsize=16)
        plt.title(col + '\neffect size versus score rho='+str(round(rho, 4))+'\n#'+str(len(score_list)), fontsize=18, y=1.025)
        ax.set_facecolor('white')
        f.savefig(output_path + "_" + col + '_EffectVsScore.'+out_type, bbox_inches='tight')
        plt.close()


def score_distribution(tag_list, score_files, ranks, output_path, colours=palette):
    """Plot a histogram of all the top fraction of scores in a gzipped ABCpp output file."""
    f, ax = plt.subplots(figsize=(10, 5))
    ax.set_axisbelow(True)
    ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='both', zorder=1)
    for n, (tag, file) in enumerate(zip(tag_list, score_files)):
        head = {x: i for i, x in enumerate(gzip.open(file, 'rt').readline().strip().split('\t'))}
        these_scores = [float(x.strip().split('\t')[head['ABC-Score']]) for x in gzip.open(file, 'rt').readlines()[1:ranks+1]]
        print("Hist", tag)
        print(np.mean(these_scores), np.std(these_scores))
        curr_median = np.median(these_scores)
        print(curr_median, np.median([abs(x - curr_median) for x in these_scores]))
        plt.hist(x=these_scores, bins=100, label=tag, alpha=0.7, color=colours[n], zorder=20)
    plt.ylabel("Occurrences", size=16)
    plt.xlabel('ABC Score', size=16)
    ax.tick_params(axis='both', labelsize=16)
    plt.title('Distribution of top '+str(ranks)+' scores', fontsize=18)
    ax.set_facecolor('white')
    ax.legend(prop={'size': 14, 'weight': 'bold'})
    f.savefig(output_path + '_ABCHist.pdf', bbox_inches='tight')
    plt.close()


def compare_scores_colour_metric(plot_df, score_cols, colour_col, marker_col=None, output_path=''):
    """Compares two scores and a chosen metric on interaction level. For each entry in plot_df plot one dot with
    [x,y] based on score_col and colour on colour_col."""
    if marker_col:
        main_list = list(zip(plot_df[score_cols[0]].values.tolist(), plot_df[score_cols[1]].values.tolist(),
                             plot_df[colour_col].values.tolist(), plot_df[marker_col].values.tolist()))
    else:
        main_list = list(zip(plot_df[score_cols[0]].values.tolist(), plot_df[score_cols[1]].values.tolist(),
                             plot_df[colour_col].values.tolist()))

    cmap = cm.get_cmap('viridis')
    norm = plt.Normalize(0, max(plot_df[colour_col]))
    # norm = plt.Normalize(0, 1000, clip=True)
    alpha_val = 0.6
    f, ax = plt.subplots(figsize=(8, 6))
    if marker_col:
        plt.scatter(x=[x[0] for x in main_list if x[3]], y=[x[1] for x in main_list if x[3]], c=[cmap(norm(x[2])) for x in main_list if x[3]],
                    s=40, edgecolors=None, zorder=12, marker='*', label=marker_col, alpha=alpha_val)
        plt.scatter(x=[x[0] for x in main_list if not x[3]], y=[x[1] for x in main_list if not x[3]], c=[cmap(norm(x[2])) for x in main_list if not x[3]],
                    s=20, edgecolors=None, zorder=12, marker='o', alpha=alpha_val)
        ax.legend(prop={'size': 12, 'weight': 'bold'})
    else:
        plt.scatter(x=[x[0] for x in main_list], y=[x[1] for x in main_list], c=[cmap(norm(x[2])) for x in main_list],
                    s=25, edgecolors=None, zorder=12, alpha=alpha_val)
    ax.set_xlabel(score_cols[0], fontsize=16)
    ax.set_ylabel(score_cols[1], fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    if 'ABC' in score_cols[0] and 'ABC' in score_cols[1]:
        upper_corner = max([ax.get_ylim()[1], ax.get_xlim()[1]])
        lower_corner = min([ax.get_ylim()[0], ax.get_xlim()[0]])
        ax.plot([lower_corner, upper_corner], [lower_corner, upper_corner], linestyle='dotted', color='grey')
    cax = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.8)
    cax.set_label(colour_col, size=14)
    f.savefig(output_path + score_cols[0] +'Vs' + score_cols[1] +'_' + colour_col + '_marker' + str(marker_col) + '.pdf', bbox_inches='tight')


