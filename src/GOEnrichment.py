from gprofiler import GProfiler
from matplotlib import cm
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt
import numpy as np


def go_enrichment(go_genes, title_tag, out_tag, max_terms=4, font_s=16, organism='hsapiens',
                  wanted_sources=['GO:MF', 'GO:BP', 'KEGG', 'REAC', 'HP']):
    """Requires a dictionary with sets of genes. Will run GO enrichment for each of the sets and allow a maximum of
    max_terms per dict key. There will be written one plot for each GO source.
    go_genes = {class: [class genes] for class in classes}
    out_tag: Will be appended by GO source + .pdf"""
    cmap = cm.get_cmap('summer')
    norm = plt.Normalize(0, 0.05)
    gp = GProfiler(return_dataframe=True)
    dict_keys = list(go_genes.keys())
    term_fetcher = {c: {s: [] for s in wanted_sources} for c in dict_keys}
    for cell in dict_keys:
        if len(go_genes[cell]) > 1:
            query_df = gp.profile(organism=organism, query=list(go_genes[cell]))
            for source in wanted_sources:
                term_fetcher[cell][source] = query_df[query_df['source'] == source]
        else:
            print(cell, 'empty gene set')
    for source in wanted_sources:
        print(source)
        term_collection = set()
        for cell in dict_keys:
            if len(term_fetcher[cell][source]) > 0:
                term_collection = term_collection.union(set(term_fetcher[cell][source]['name'].to_list()[:max_terms]))

        these_terms = list(term_collection)
        # First count in how many cell types a term is present and sort according to that.
        term_occs = []
        for term in these_terms:
            hits = 0
            for cell in dict_keys:
                curr_df = term_fetcher[cell][source]
                if len(curr_df) > 0:
                    if len(curr_df[curr_df['name'] == term]) == 1:
                        hits += 1
            term_occs.append([term, hits])
        sorted_terms = [x[0] for x in sorted(term_occs, key=lambda x: x[1])]

        main_list = []  # X-pos, Y-pos, gene fraction, p-value.
        for x, cell in enumerate(dict_keys):
            num_genes = len(go_genes[cell])
            for y, term in enumerate(sorted_terms):
                curr_df = term_fetcher[cell][source]
                if len(curr_df) > 0:
                    match_entry = curr_df[curr_df['name'] == term]
                    if len(match_entry) == 1:
                        main_list.append([x, y, match_entry['intersection_size'].values[0] / num_genes,
                                          match_entry['p_value'].values[0]])
        if main_list:
            size_col = [s[2] for s in main_list]
            scaled_min, scaled_max = 20, 200
            if len(main_list) > 1:
                scaled_sizes = [((scaled_max - scaled_min) * (s - min(size_col)) / (max(size_col) - min(size_col))) + scaled_min
                                for s in size_col]
            else:
                scaled_sizes = [scaled_max]
            format_terms = []
            for term in sorted_terms:
                if len(term) > 40 and ' ' in term:
                    closest_blank = \
                    sorted([[i, abs(i - len(term) // 2)] for i, x in enumerate(term) if x == ' '], key=lambda x: x[1])[0][0]
                    format_terms.append(term[:closest_blank] + '\n' + term[closest_blank + 1:])
                elif len(term) > 40:
                    format_terms.append(term[:len(term) // 2] + '-\n' + term[len(term) // 2:])
                else:
                    format_terms.append(term)

            fig_height = int(len(sorted_terms) * 0.4)
            f, ax = plt.subplots(figsize=(3 + len(dict_keys), fig_height if fig_height > 6 else 6))
            ax.yaxis.grid(True, which='major', color='#cccccc', zorder=1)
            ax.set_ylim(-0.5, len(these_terms) + 0.5)
            plt.scatter(x=[x[0] for x in main_list], y=[y[1] for y in main_list], c=[cmap(norm(c[3])) for c in main_list],
                        s=scaled_sizes, edgecolors='k', zorder=12)
            ax.set_yticks([y for y in range(len(sorted_terms))])
            ax.set_yticklabels(format_terms, fontweight='bold', size=font_s - 8)
            ax.set_xticks([x for x in range(len(dict_keys))])
            if np.mean([len(c) for c in dict_keys]) > 4:
                ax.set_xticklabels([c for c in dict_keys], size=font_s - 6, ha='right', fontweight='bold', rotation=45)
            else:
                ax.set_xticklabels([c for c in dict_keys], size=font_s - 6, ha='center', fontweight='bold')
            cax = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.5)
            cax.set_label('adj. p-value', size=font_s - 4)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            if len(scaled_sizes) > 2:
                h1 = Line2D([0], [0], marker='o', markersize=np.sqrt(scaled_min), color='black', linestyle='None')
                h2 = Line2D([0], [0], marker='o', markersize=np.sqrt((scaled_min + scaled_max) / 2), color='black',
                            linestyle='None')
                h3 = Line2D([0], [0], marker='o', markersize=np.sqrt(scaled_max), color='black', linestyle='None')
                plt.legend([h1, h2, h3], [str(round(min(size_col), 4)), str(round((min(size_col) + max(size_col)) / 2, 4)),
                                          str(round(max(size_col), 4))], loc="upper right", markerscale=1, scatterpoints=1,
                           fontsize=font_s - 4, title="gene fraction", bbox_to_anchor=(1.5, 1))
            else:
                h3 = Line2D([0], [0], marker='o', markersize=np.sqrt(scaled_max), color='black', linestyle='None')
                plt.legend([h3], [str(round(max(size_col), 4))], loc="upper right", markerscale=1, scatterpoints=1,
                           fontsize=font_s - 6, title="gene fraction", bbox_to_anchor=(1.5, 1))
            for c in range(len(dict_keys)):
                plt.axvline(c, color='#858585')
            # Simulate background stripes to find the cells easier.
            ax.set_title(source + ": " + title_tag, size=font_s)

            f.savefig((out_tag + "_" + source + "_max"+str(max_terms)+".pdf").replace(' ', ''),
                      bbox_inches='tight')
            plt.close()
