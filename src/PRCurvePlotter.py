import pandas as pd
from matplotlib import pyplot as plt
import sklearn.metrics as metrics
from matplotlib import cm
import numpy as np


def pre_rec_fetcher(binning, ele_gene_pairs, sig_col, score_col, lower_limit=0, upper_limit=1, mark_thresh=None):

    def fetch_prerec(thresholds):
        rec_prec_list = []
        recall_hit = []
        for thresh in thresholds:
            # How many of the positives did we find.
            if ele_gene_pairs.dtypes[sig_col] in [bool, 'int64']:
                total_positives = ele_gene_pairs[ele_gene_pairs[sig_col] == 1].shape[0]
                true_positives = ele_gene_pairs[(ele_gene_pairs[sig_col] == 1) & (pd.to_numeric(ele_gene_pairs[score_col]) >= thresh)]
                false_positives = ele_gene_pairs[(ele_gene_pairs[sig_col] == 0) & (pd.to_numeric(ele_gene_pairs[score_col]) >= thresh)]
            else:
                total_positives = ele_gene_pairs[ele_gene_pairs[sig_col] == 'True'].shape[0]
                true_positives = ele_gene_pairs[(ele_gene_pairs[sig_col] == 'True') & (pd.to_numeric(ele_gene_pairs[score_col]) >= thresh)]
                false_positives = ele_gene_pairs[(ele_gene_pairs[sig_col] == 'False') & (pd.to_numeric(ele_gene_pairs[score_col]) >= thresh)]

            recall = true_positives.shape[0] / total_positives
            if true_positives.shape[0] != 0:
                precision = true_positives.shape[0] / (true_positives.shape[0] + false_positives.shape[0])
            else:
                precision = 0
            rec_prec_list.append([recall, precision, thresh])
            if round(recall, 2) == 0.70:
                recall_hit.append(thresh)
        # print(score_col, '70% recall', round(np.mean(recall_hit), 4))
        return rec_prec_list

    binned_list = fetch_prerec([lower_limit + (upper_limit-lower_limit) * x / binning for x in range(0, binning+1)])
    if mark_thresh:
        singular_pos = fetch_prerec([mark_thresh])
    else:
        singular_pos = None
    return binned_list, singular_pos


def pr_thresh_coloured(binning, plot_df, score_cols, sig_col, thresh, output_path):
    """Plot PR curves, but as scatter and coloured by the threshold."""
    marker_shapes = ['^', 'v', 'P', '*', 'h']
    lower_lim = min([plot_df[c].min() for c in score_cols])
    upper_lim = max([plot_df[c].max() for c in score_cols])
    cmap = cm.get_cmap('viridis')
    norm = plt.Normalize(lower_lim, upper_lim)
    f, ax = plt.subplots(1, figsize=(9, 8))
    plt.subplots_adjust(left=0.08, right=0.98)
    ax.set_axisbelow(True)
    ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='both')
    auc_coll = []
    for c, col in enumerate(score_cols):
        rec_prec_list, thresh_pos = pre_rec_fetcher(binning, plot_df, sig_col, score_col=col, lower_limit=lower_lim,
                                                    upper_limit=upper_lim, mark_thresh=thresh)
        sorted_by_recall = sorted(rec_prec_list, key=lambda x: x[0])
        auc = metrics.auc([x[0] for x in sorted_by_recall], [x[1] for x in sorted_by_recall])
        auc_coll.append(col + ' AUC: ' + str(round(auc, 3)))
        print(col, round(auc, 3))
        plt.scatter([x[0] for x in rec_prec_list], [x[1] for x in rec_prec_list], c=[cmap(norm(x[2])) for x in rec_prec_list],
                    s=45, zorder=12, marker=marker_shapes[c], label=col + ' ' + str(round(auc, 3)), linewidths=0)
        plt.plot([x[0] for x in rec_prec_list], [x[1] for x in rec_prec_list], linestyle='-', linewidth=1, c='#e3e2e1')
    cax = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.8)
    cax.set_label('score cut-off', size=16)
    cax.ax.tick_params(labelsize=14)
    plt.xlabel('Recall', size=20)
    plt.ylabel('Precision', size=20)
    ax.tick_params(axis='both', labelsize=16)
    # ax.set_title('\n'.join(auc_coll))
    ax.legend(markerscale=3, scatterpoints=1, fontsize=16, labelspacing=0.4, prop={'weight': 'bold', 'size': 18})

    plt.ylim(0, 1)
    plt.xlim(0, 1)
    # plt.title('PR curve coloured by score threshold', fontsize=18)
    ax.set_facecolor('white')
    f.savefig(output_path + '_'.join(score_cols).replace(' ', '') + '_' + sig_col.replace(' ', '') + '_PRColoured.pdf', bbox_inches='tight')

