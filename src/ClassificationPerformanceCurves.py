import pandas as pd
from matplotlib import pyplot as plt
import sklearn.metrics as metrics
from matplotlib import cm
import numpy as np
import random
import copy
import ColoursAndShapes

"""Collection of functions to create precision-recall curves or receiver operator characteristic curves and
calculating the area under them."""


def performance_fetcher(binning, entries_df, true_col="Significant", score_col='ABC Score', lower_limit=0, upper_limit=1,
                    mode='pr'):
    """
    Helper function to iterate from lower_limit to upper_limit in a given number of steps and calculating the
    precision and recall for each step, or the false-positive rate and recall.

    Args:
        binning: number of steps between lower_limit and upper_limit. For each step the entries above the threshold are
            considered to be predicted as True.
        entries_df: DataFrame with the entries for which the classification was run.
        true_col: Name of the column with booleans to find the true entries.
        score_col: Column with a score for which the curves should be calculated.
        lower_limit: From where the binning should start.
        upper_limit: Where the binning should end.
        mode: 'pr' for Precision-Recall, otherwise getting FPR-TPR (ROC).

    Returns:
        list:
            - **performance_list**: List with an entry for each tested threshold with [Recall, Precision, threshold] for mode 'pr', otherwise [FPR, TPR, threshold].
        """

    thresholds = [lower_limit + (upper_limit-lower_limit) * x / binning for x in range(0, binning+1)]

    performance_list = []
    for thresh in thresholds:
        true_positives = entries_df[(entries_df[true_col]) & (pd.to_numeric(entries_df[score_col]) >= thresh)]
        recall = true_positives.shape[0] / entries_df[true_col].sum()
        false_positives = entries_df[(~entries_df[true_col]) & (pd.to_numeric(entries_df[score_col]) >= thresh)]
        negatives = (~entries_df[true_col]).sum()

        if true_positives.shape[0] != 0:
            precision = true_positives.shape[0] / (true_positives.shape[0] + false_positives.shape[0])
        else:
            precision = 0
        if negatives > 0:
            fpr = len(false_positives) / negatives
        else:
            fpr = 0
        if mode == 'pr':  # Precision andRecall
            performance_list.append([recall, precision, thresh])
        else:  # ROC: FPR and TPR
            performance_list.append([fpr, recall, thresh])

    return performance_list


def classification_plotter(df, sig_col, score_cols, add_random=False, steps=10000, mode='pr', output_path='', colours='glasbey',
                           no_plot=False, recall_start=0, zorder=None, title_tag='', legend_s=18, legend_out=False, x_size=8,
                           y_size=8, font_s=14, line_styles=None, colour_by_threshold=False, threshold_cmap='viridis',
                           formats=['pdf']):
    """
    Plots a Precision-Recall curve or a receiver operator characteristic curve based on a DataFrame and calculates the
    area under the curve.

    Args:
        df: Pandas DataFrame with each row being an entry that should be classified.
        sig_col: Name of the column that identifies the true entries.
        score_cols: List of column names in the DataFrame for which the curves should be plotted, all in the same plot. It is assumed that a high score means a higher predicted likelihood to be true.
        add_random: If a curve should be added that randomly orders the entries. For specifying the colour of 'Random', add a colour to the colour list, otherwise it will be grey.
        steps: Number of steps into which the range between the lowest and highest score will be separated, and for each the performance calculated.
        mode: 'pr' to get a Precision-Recall curve, otherwise a ROC curve.
        no_plot: To only get the list of performance values.
        recall_start: In case it is known that a certain range of the recall is not covered, limit the whole calculation and plotting to [recall_start, 1].
        zorder: List of integers defining the zorder of the score_cols.
        colour_by_threshold: If True, do the plot as scatter, and colour each dot by the threshold it was taken from. Uses the range from all score_cols.
        threshold_cmap: The colourmap which should be used for the scatter when colour_by_threshold is True.

    Returns:
        tuple:
            - **auc_output**: List of the score_cols and the respective AUPRC.
            - **performance_dict**: Dictionary with {score_col: [Recall, Precision, threshold] for mode 'pr', otherwise [FPR, TPR, threshold] for all tested thresholds}.
    """
    if colours and 'glasbey' in colours:
        colours = ColoursAndShapes.glasbey_palettes[colours]
    if not line_styles:
        line_styles = ['solid', 'dashed', 'dotted', 'dashdot'] * (int(np.ceil(len(score_cols) // 4 + 1)))

    if colour_by_threshold:
        whole_min = df[score_cols].min().min()
        whole_max = df[score_cols].max().max()
        cmap = cm.get_cmap(threshold_cmap)
        norm = plt.Normalize(whole_min, whole_max)

    plot_df = copy.deepcopy(df)
    if plot_df.dtypes[sig_col] != bool:
        plot_df[sig_col] = plot_df[sig_col].isin(['True', 'TRUE', 'true', 'T', '1'])
    if add_random:
        random.seed(1234)
        plot_df['Random'] = random.sample(range(plot_df.shape[0]), plot_df.shape[0])
        score_cols += ['Random']
        colours += ['#767B8E']
        line_styles += ['dotted']
    performance_dict = {c: None for c in score_cols}

    f, ax = plt.subplots(1, figsize=(x_size, y_size))
    plt.subplots_adjust(left=0.08, right=0.98)
    ax.set_axisbelow(True)
    ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='both')
    auc_coll = []
    auc_output = []
    if not zorder:
        zorder = [i+10 for i in range(len(score_cols))]
    for n, this_col in enumerate(score_cols):
        lower_lim = min(plot_df[this_col])
        upper_lim = max(plot_df[this_col])
        performance_list = performance_fetcher(steps, plot_df, sig_col, score_col=this_col, lower_limit=lower_lim,
                                               upper_limit=upper_lim, mode=mode)
        performance_dict[this_col] = performance_list
        sorted_by_recall = sorted(performance_list, key=lambda x: x[0])
        sorted_by_recall = [x for x in sorted_by_recall if round(x[0], 3) >= recall_start]
        auc = metrics.auc([x[0] for x in sorted_by_recall], [x[1] for x in sorted_by_recall])
        auc_coll.append(this_col + ' AUC: ' + str(round(auc, 3)))
        auc_output.append([this_col, auc])
        print(this_col, round(auc, 4))

        if colour_by_threshold:
            plt.scatter([x[0] for x in performance_list], [x[1] for x in performance_list],
                        c=[cmap(norm(x[2])) for x in performance_list],
                        s=16, zorder=12, marker=ColoursAndShapes.marker_shapes[n], label=this_col)
        else:
            plt.plot([x[0] for x in performance_list], [x[1] for x in performance_list],
                     linestyle=line_styles[n], c=colours[n],
                     label=this_col, linewidth=4, zorder=zorder[n])

    if colour_by_threshold:
        cax = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.8)
        cax.set_label('score cut-off', size=font_s+4)
        cax.ax.tick_params(labelsize=font_s+4)
        ax.legend(markerscale=3, scatterpoints=1, fontsize=font_s+2, labelspacing=0.4, prop={'weight': 'normal', 'size': font_s+4})

    if legend_out:
        ax.legend(prop={'size': legend_s, 'weight': 'normal'}, loc='upper right',
                  bbox_to_anchor=(2 if type(legend_out) == bool and not legend_out else legend_out, 1))
    else:
        ax.legend(prop={'size': legend_s, 'weight': 'normal'})
    plt.xlabel('Recall' if mode == 'pr' else 'FPR', size=font_s+8)
    plt.ylabel('Precision' if mode == 'pr' else 'TPR', size=font_s+8)
    ax.tick_params(axis='both', labelsize=18)
    plt.ylim(0, 1)
    plt.xlim(recall_start, 1)
    ax.set_title(title_tag+'\n'+str(plot_df[sig_col].sum())+' true/ '+str(len(plot_df)) + ' interactions' + '\n' + '\n'.join(auc_coll),
                 fontsize=font_s+2, y=1.05)
    ax.set_facecolor('white')
    if not no_plot:
        if type(formats) != list:
            formats = [formats]
        for form in formats:
            f.savefig(output_path + '_' + mode + ('_ColouredbyThreshold'*colour_by_threshold) + '.'+form, bbox_inches='tight', format=form)
    plt.close()

    return auc_output, performance_dict


