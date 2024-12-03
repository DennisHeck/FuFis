import seaborn as sns
from scipy import stats
import itertools
import copy
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import to_hex
import ColoursAndShapes

font_s = 16
# red_colours = ['#2e0000', '#a11b1b', '#ed6109', '#ffd86b', '#ffc400']
# blue_colours = ['#00003b', '#2c46b0', '#1231c9', '#46cdf2', '#ffc400']
# default_c = red_colours + blue_colours
qual_map = plt.get_cmap('tab10', 10)
default_c = [to_hex(qual_map(i)) for i in range(10)]


def ks_test(curr_df, grouping_col, value_col):
    """Runs a Kolmogorov smirnov test for each non-redundant pairwise comparison."""
    samples = sorted(curr_df[grouping_col].drop_duplicates().tolist())
    order_dict = {x: i for i, x in enumerate(samples)}
    comparisons = list(itertools.combinations(samples, 2))
    p_list = []
    for comp in comparisons:
        sample1_list = curr_df[curr_df[grouping_col] == comp[0]][value_col].values.tolist()
        sample2_list = curr_df[curr_df[grouping_col] == comp[1]][value_col].values.tolist()
        p_list.append([comp, stats.ks_2samp(sample1_list, sample2_list)[1]])
    p_nested = [[0 for x in range(y + 1)] for y in reversed(list(range(len(samples) - 1)))]
    p_background = [['white' for x in range(y + 1)] for y in reversed(list(range(len(samples) - 1)))]
    print(p_list)
    for i, p in enumerate(p_list):
        p_nested[order_dict[p[0][0]]][abs(order_dict[p[0][1]] - (len(samples) - 1))] = round(p[1], 5)
        if p[1] <= 0.05:
            p_background[order_dict[p[0][0]]][abs(order_dict[p[0][1]] - (len(samples) - 1))] = '#ed5739'
    for i, p in enumerate(p_nested):
        p_nested[i] += ['' for x in range((len(samples) - 1) - len(p_nested[i]))]
        p_background[i] += ['white' for x in range((len(samples) - 1) - len(p_background[i]))]

    return p_nested, p_background, samples


def cumu_plot(plot_df, fetch_col, grouping_name, output_path, table_width=0.3, table_x_pos=1.2, xlimit=None,
              vertical_line=False, hue_order=None, numerate=False, title=None):
    """Plot the cumulative distribution of all sets of grouping names in the plotting df.
    Adds a table with a Kolmogorov-Smirnov test for each non-redunant pairwise comparison next to the plot."""
    p_df = plot_df.copy()
    sns.set_style('ticks')
    sns.set_context('talk', font_scale=1.15)
    if hue_order:
        gene_set_cols = hue_order
    else:
        gene_set_cols = p_df[grouping_name].drop_duplicates().to_list()
    gene_set_cols_num = []
    different_groups = len(gene_set_cols)
    if different_groups == 3:
        linestyles = ['solid', 'dashed', 'dashdot', 'solid']
    else:
        linestyles = ['solid', 'dashed', 'solid', 'dashdot']

    if numerate:
        hue_order = []
        for i, x in enumerate(gene_set_cols):
            gene_set_cols_num.append(str(gene_set_cols[i]) + ' (#' + str(sum(p_df[grouping_name] == x)) + ')')
            hue_order.append(str(gene_set_cols[i]) + ' (#' + str(sum(p_df[grouping_name] == x)) + ')')
        for i, x in enumerate(gene_set_cols):
            p_df.loc[p_df[grouping_name] == x, grouping_name] = gene_set_cols_num[i]
    if different_groups == 2:
        palette = ColoursAndShapes.two_contrasts[0]
    elif different_groups <= len(ColoursAndShapes.tol_vibrant):
        palette = ColoursAndShapes.tol_vibrant[:different_groups]
    else:
        palette = [to_hex(c) for c in ColoursAndShapes.glasbey_palettes['glasbey']][:different_groups]
    cumu = sns.displot(p_df, x=fetch_col, hue=grouping_name, kind='ecdf', hue_order=hue_order,
                       palette=palette, height=8, aspect=1.3, zorder=12)
    # cumu._legend.set_title('')
    cumu.fig.subplots_adjust(top=0.93)
    cumu.fig.suptitle('Cumulative distribution of ' + fetch_col if not title else title, size=24, x=0.28, y=1.01)
    cumu.ax.set_xlabel(fetch_col, size=22, fontweight='bold')
    cumu.ax.set_ylabel('Proportion', size=22, fontweight='bold')
    if xlimit:
        cumu.set(xlim=(xlimit[0], xlimit[1]))
    colours_order = []
    for a in range(len(cumu.ax.lines)):
        colours_order.append(cumu.ax.get_lines()[a].get_color().lower())
        cumu.ax.lines[a].set_linestyle(linestyles[a % len(linestyles)])
        cumu.ax.lines[a].set_linewidth('4')
    colour_dict = {x: i for i, x in enumerate(colours_order)}
    print(colour_dict)
    print(different_groups)
    for a in range(len(cumu.ax.lines)):
        if to_hex(cumu.legend.get_lines()[a].get_color()) in colour_dict:  # Legend has it even if there was no data.
            cumu.legend.get_lines()[a].set_linestyle(linestyles[colour_dict[to_hex(cumu.legend.get_lines()[a].get_color())] % len(linestyles)])
            cumu.legend.get_lines()[a].set_linewidth('3')
    for spine in cumu.ax.spines.values():  # Increase spline width.
        spine.set_linewidth(3)
    p_nested, p_background, samples = ks_test(p_df, grouping_name, fetch_col)
    samples = [str(x).split('(')[0][:int(len(str(x).split('(')[0]) / 2)] + '\n' + str(x).split('(')[0][int(len(str(x).split('(')[0]) / 2):]
               for x in samples]
    if len(samples) > 1:
        ks_table = plt.table(cellText=np.asarray(p_nested).T, cellColours=np.asarray(p_background).T,
                             rowLabels=list(reversed(samples))[:-1], colLabels=samples[:-1], loc='bottom right',
                             cellLoc='center', rowLoc='center', zorder=12, bbox=[table_x_pos, -0.11, table_width, 0.35])
        ks_table.auto_set_font_size(False)
        ks_table.set_fontsize(10)

    cumu.ax.grid(True, axis='both', color='#c7c7c7', linewidth=1, which='major')
    cumu.ax.set_facecolor('white')
    if vertical_line is not None:
        cumu.ax.axvline(x=vertical_line, color='#6e6e6e', zorder=1)

    cumu.savefig(output_path+fetch_col.replace(' ', '')+'_'+grouping_name+'.pdf', bbox_inches='tight')
    plt.close()
