from gprofiler import GProfiler
from matplotlib import cm
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt
import numpy as np


def go_enrichment(go_genes, title_tag='', out_tag='', max_terms='all', font_s=16, organism='hsapiens', numerate=False,
                  background=None, godf_only=False, cmap='plasma', fig_width=None, fig_height=None, legend_out=None,
                  wanted_sources=['GO:MF', 'GO:BP', 'KEGG', 'REAC', 'HP', 'WP'], keywords={}, rotation=45,
                  custom_dfs=None):
    """Requires a dictionary with sets of genes. Will run GO enrichment for each of the sets. There will be written
    one plot for each GO source.
    @param max_terms: Allow a maximum of max_terms per dict key. Use 'all' to get all.
    @param go_genes: {class: [class genes] for class in classes}
    @param background: Same as go_genes but holding the set of background genes for each entry in go_genes.
    @param custom_dfs: {gene set: {go_source: DF}}; In case an external df was used or a g:Profiler was customized,
    to just use the plotting part. Must have the same format as the g:Profiler DFs.
    @param out_tag: Will be appended by GO source + .pdf"""
    keywords = {k: [t.lower() for t in vals] for k, vals in keywords.items()}  # We compare lower strings.
    cmap = cm.get_cmap(cmap)
    norm = plt.Normalize(0, 0.05)
    dict_keys = list(go_genes.keys())
    term_fetcher = {c: {s: [] for s in wanted_sources} for c in dict_keys}
    df_fetcher = {c: None for c in dict_keys}

    if not custom_dfs:
        gp = GProfiler(return_dataframe=True)
        for cell in dict_keys:
            if len(go_genes[cell]) >= 1:
                query_df = gp.profile(organism=organism, query=list(go_genes[cell]), no_evidences=False,
                                      background=None if not background else list(background[cell]))
                query_df = query_df[(~query_df['name'].str.contains("REACTOME")) & (query_df['name'] != 'WIKIPATHWAYS')
                                    & (~query_df['name'].str.contains("KEGG")) & (~query_df['name'].str.contains('HP root'))]
                query_df['Gene fraction'] = query_df['intersection_size'] / len(go_genes[cell])
                df_fetcher[cell] = query_df
                for source in wanted_sources:
                    term_fetcher[cell][source] = query_df[query_df['source'] == source]
            else:
                print(cell, 'empty gene set')
        if godf_only:
            return df_fetcher
    else:
        term_fetcher = custom_dfs

    for source in wanted_sources:
        print(source)
        term_collection = set()
        for cell in dict_keys:
            if len(term_fetcher[cell][source]) > 0:
                if str(max_terms).lower() == 'all':
                    term_collection = term_collection.union(set(term_fetcher[cell][source]['name'].to_list()))
                else:
                    term_collection = term_collection.union(set(term_fetcher[cell][source]['name'].to_list()[:max_terms]))

        these_terms = list(term_collection)
        if source in keywords:
            these_terms = [x for x in these_terms if np.any([k.lower() in x.lower() for k in keywords[source]])]
        # First count in how many cell types a term is present and sort according to that.
        term_occs = []
        for term in these_terms:
            hits = 0
            pvals = []
            for cell in dict_keys:
                curr_df = term_fetcher[cell][source]
                if len(curr_df) > 0:
                    if len(curr_df[curr_df['name'] == term]) == 1:
                        hits += 1
                        pvals.append(next(iter(curr_df[curr_df['name'] == term]['p_value'].values)))
            term_occs.append([term, hits, min(pvals)])
        sorted_terms = [x[0] for x in sorted(term_occs, key=lambda x: (x[1], -x[2]))]

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
            if len(main_list) > 1 and min(size_col) != max(size_col):  # Otherwise all get the same size.
                scaled_sizes = [((scaled_max - scaled_min) * (s - min(size_col)) / (max(size_col) - min(size_col))) + scaled_min
                                for s in size_col]
            else:
                scaled_sizes = [scaled_max] * len(size_col)
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
            # If we have multiple gene groups, do a bubble plot aligning groups vertically.
            this_height = max([int(len(sorted_terms) * 0.4), 6]) if not fig_height else fig_height
            f, ax = plt.subplots(figsize=(fig_width if fig_width else (3 + len(dict_keys)), this_height))
            ax.yaxis.grid(True, which='major', color='#f0eded', zorder=1)
            ax.set_ylim(-0.5, len(these_terms) - 0.5)
            if len(go_genes) > 1:
                ax.set_xlim(-0.2, len(dict_keys) - 0.9)
                plt.scatter(x=[x[0] for x in main_list], y=[y[1] for y in main_list], c=[cmap(norm(c[3])) for c in main_list],
                            s=scaled_sizes, edgecolors=[cmap(norm(c[3])) for c in main_list], zorder=12)

                ax.set_xticks([x for x in range(len(dict_keys))])
                ax.set_xticklabels([c for c in dict_keys], size=font_s - 6, ha='right' if rotation else 'center',
                                       fontweight='bold', rotation=rotation)
                if numerate:
                    x_counts = {c: len(vals) for c, vals in go_genes.items()}
                    ax.set_xticklabels([x._text + '\n(#' + str(x_counts[x._text]) + ')' for x in ax.get_xmajorticklabels()])
                cax = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.5)
                x_offset = cax.ax.get_position().x1
                cax.set_label('adj. p-value', size=font_s - 4)
                for c in range(len(dict_keys)):
                    plt.axvline(c, color='#c9c7c7', ymax=len(these_terms)-1, linewidth=0.5)
                ax.set_title(source + ": " + title_tag, size=font_s)

            else:  # If we only have one group, use the x-axis for the log-q value.
                ax.set_xlim(3, max([abs(np.log2(x[3])) for x in main_list])*1.1)
                plt.scatter(x=[abs(np.log2(x[3])) for x in main_list], y=[y[1] for y in main_list],
                            c=['#112791' for c in main_list],  # cmap(norm(c[3]))
                            s=scaled_sizes, edgecolors=['#112791' for c in main_list], zorder=12)
                ax.set_yticks([y for y in range(len(sorted_terms))])
                ax.set_yticklabels(format_terms, fontweight='bold', size=font_s - 8)
                ax.set_xlabel('-log2 FDR', size=font_s - 4)
                ax.set_title(
                    source + ": " + title_tag + (' (#' + str(len(next(iter(go_genes.values())))) + ')') * numerate,
                    size=font_s)
                x_offset = 1.1

            ax.set_yticks([y for y in range(len(sorted_terms))])
            ax.set_yticklabels(format_terms, fontweight='bold', size=font_s - 8)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            if len(scaled_sizes) > 2:
                h1 = Line2D([0], [0], marker='o', markersize=np.sqrt(scaled_min), color='black', linestyle='None')
                h2 = Line2D([0], [0], marker='o', markersize=np.sqrt((scaled_min + scaled_max) / 2), color='black',
                            linestyle='None')
                h3 = Line2D([0], [0], marker='o', markersize=np.sqrt(scaled_max), color='black', linestyle='None')
                leg = plt.legend([h1, h2, h3], [str(round(min(size_col), 4)), str(round((min(size_col) + max(size_col)) / 2, 4)),
                                          str(round(max(size_col), 4))], loc="upper right", markerscale=1, scatterpoints=1,
                           fontsize=font_s - 6, title="Gene fraction", bbox_to_anchor=(legend_out if legend_out else x_offset+0.5, 1))
            else:
                h3 = Line2D([0], [0], marker='o', markersize=np.sqrt(scaled_max), color='black', linestyle='None')
                leg = plt.legend([h3], [str(round(max(size_col), 4))], loc="upper right", markerscale=1, scatterpoints=1,
                           fontsize=font_s - 6, title="Gene fraction", bbox_to_anchor=(legend_out if legend_out else x_offset+0.5, 1))
            leg.get_title().set_fontsize(11)

            f.savefig((out_tag + "_" + source + "_max"+str(max_terms)+".pdf").replace(' ', ''),
                      bbox_inches='tight')
            plt.close()

        else:
            print("No enrichment", source)

    return df_fetcher
