from gprofiler import GProfiler
from matplotlib import cm
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt
import copy
import numpy as np
import pandas as pd
import gseapy as gp
import GTF_Processing
import Various


def plot_go(mode, wanted_sources, term_fetcher, keywords, cmap, fig_width, fig_height, go_genes, dict_keys, 
            numerate, font_s, rotation, title_tag, legend_out, max_terms, out_tag, formats):
    """Takes care of processing and plotting the GO enrichment results, to be able to plot results from different origins."""

    for source in wanted_sources:
        print(source)
        term_collection = set()
        for g_set in dict_keys:
            if len(term_fetcher[g_set][source]) > 0:
                if str(max_terms).lower() == 'all':
                    term_collection = term_collection.union(set(term_fetcher[g_set][source]['name'].values))
                else:
                    term_collection = term_collection.union(set(term_fetcher[g_set][source]['name'].values[:int(max_terms)]))

        these_terms = list(term_collection)
        if source in keywords:
            these_terms = [x for x in these_terms if np.any([k.lower() in x.lower() for k in keywords[source]])]
        # First count in how many gene sets a term is present and sort according to that and secondarily either by
        # minimum p-value or the NES.
        term_occs = []
        for term in these_terms:
            hits = 0
            second_sorter = []
            for g_set in dict_keys:
                curr_df = term_fetcher[g_set][source]
                curr_df = curr_df[curr_df['name'] == term]
                if not curr_df.empty:
                    hits += 1
                    if mode == 'hypergeometric':
                        second_sorter.append(-next(iter(curr_df['FDR'].values)))
                    else:
                        second_sorter.append(next(iter(curr_df['NES'].values)))
            term_occs.append([term, hits, min(second_sorter) if mode == 'hypergeometric' else np.mean(second_sorter)])
        sorted_terms = [x[0] for x in sorted(term_occs, key=lambda x: (x[1], x[2]))]
        
        main_list = []  # X-pos, Y-pos, gene fraction, p-value.
        for x, g_set in enumerate(dict_keys):
            for y, term in enumerate(sorted_terms):
                curr_df = term_fetcher[g_set][source]
                if len(curr_df) > 0:
                    match_entry = curr_df[curr_df['name'] == term]
                    if not match_entry.empty:
                        if mode == 'hypergeometric':
                            main_list.append([x, y, match_entry['Gene fraction'].values[0],
                                            abs(np.log2(match_entry['FDR'].values[0]))])
                        elif mode == 'prerank':
                            main_list.append([x, y, match_entry['Gene fraction'].values[0],
                                              match_entry['NES'].values[0]])
                        else:
                            print("ERROR: unknown mode "+mode+" only 'hypergeometric' or 'prerank' allowed.")
        if main_list:
            if mode == 'hypergeometric':
                cmap = cm.get_cmap(cmap)
                norm = plt.Normalize(0, 0.05)
            else:
                cmap = cm.get_cmap('bwr')  # Fixed since we can have positive and negative.
                boundary = max([abs(x[3]) for x in main_list])
                norm = plt.Normalize(-boundary, boundary)

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
                if mode == 'hypergeometric':
                    plt.scatter(x=[x[3] for x in main_list], y=[y[1] for y in main_list],
                            c=['#112791' for c in main_list], s=scaled_sizes, 
                            edgecolors=['#112791' for c in main_list], zorder=12)
                    ax.set_xlim(3, max([x[3] for x in main_list])*1.1)
                    ax.set_xlabel('-log2 FDR', size=font_s - 4)

                else:
                    plt.scatter(x=[x[3] for x in main_list], y=[y[1] for y in main_list],
                            c=[cmap(norm(c[3])) for c in main_list], s=scaled_sizes, 
                            edgecolors=[cmap(norm(c[3])) for c in main_list], zorder=12)
                    minx = min([x[3] for x in main_list])
                    maxx = max([x[3] for x in main_list])
                    ax.set_xlim(minx*(1.2 if minx < 0 else 0.8), maxx*(1.2 if maxx > 0 else 0.8))
                    ax.set_xlabel('NES', size=font_s - 4)
                    plt.axvline(0, color="#a7a8a7", linestyle="--")
                    
                    cax = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.5)
                    x_offset = cax.ax.get_position().x1
                    cax.set_label('NES', size=font_s - 4)
                    ax.set_title(source + ": " + title_tag, size=font_s)
                
                ax.set_yticks([y for y in range(len(sorted_terms))])
                ax.set_yticklabels(format_terms, fontweight='bold', size=font_s - 8)
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
            if type(formats) != list:
                formats = [formats]
            for form in formats:
                f.savefig((out_tag + "_" + source + "_max"+str(max_terms)+"."+form).replace(' ', '').replace(':', ''),
                            bbox_inches='tight', format=form)
            plt.close()
        else:
            print("No enrichment", source)


def go_enrichment(go_genes, title_tag='', out_tag='', max_terms='all', organism='hsapiens', background=None,
                  numerate=False, godf_only=False, wanted_sources=['GO:MF', 'GO:BP', 'KEGG', 'REAC', 'HP', 'WP'],
                  keywords={}, cmap='plasma', fig_width=None, fig_height=None, legend_out=None, rotation=45, font_s=16,
                  custom_dfs=None, formats=['pdf']):
    """Requires a dictionary with sets of genes. Will run GO enrichment for each of the sets. One plot will be written
    for each GO source. For only one gene set uses the x-axis for indicating the FDR, multiple sets will be
    separated on the x-axis and the FDR value shown as colour. Note, root terms of the databases are manually
    filtered out (e.g., HP root), but some might be missed. If there are multiple terms with the same name, the first
    one is picked. Uses the Python package of g:Profiler:
     https://biit.cs.ut.ee/gprofiler/gost.

    Args:
        go_genes: {class: [class genes] for class in classes}, or as sets. Safest identifiers are Ensembl IDs.
        max_terms: Allow a maximum of max_terms per dict key. Use 'all' to get all.
        organism: 'hsapiens' for human, 'mmusculus' for mouse. Check the g:Profiler page for all options.
        background: Same as go_genes but holding the set of background genes for each entry in go_genes. The keys of
            the two dictionaries will be matched. Leave empty to not use any background.
        numerate: Show the size of the gene sets in parentheses on the x-axis.
        godf_only: Skip the plotting step and only return the df.
        wanted_sources: Which databases to plot.
        keywords: A dictionary of {source: list of keywords}, e.g. {GO:BP: ['vascular', 'angio']}, to limit the
            plot to terms containing any of the listed strings. Capitalization doesn't matter, string comparison is done
            on the lowered strings.
        custom_dfs: Very specific use case. {gene set: {go_source: DF}}; In case an external df was used or
            a g:Profiler was customized, to just use the plotting part. Must have the same format as the g:Profiler DFs.

    Returns:
        - Returns a dict of {key df} with a df for each gene set as provided by the gprofiler package, which includes which genes matched to which terms.
    """
    keywords = {k: [t.lower() for t in vals] for k, vals in keywords.items()}  # We compare lower strings.

    dict_keys = list(go_genes.keys())
    term_fetcher = {c: {s: [] for s in wanted_sources} for c in dict_keys}
    df_fetcher = {c: None for c in dict_keys}

    if not custom_dfs:
        gp = GProfiler(return_dataframe=True)
        for g_set in dict_keys:
            if len(go_genes[g_set]) >= 1:
                query_df = gp.profile(organism=organism, query=list(go_genes[g_set]), no_evidences=False,
                                      background=None if not background else list(background[g_set]))
                query_df = query_df[(~query_df['name'].str.contains("REACTOME")) & (query_df['name'] != 'WIKIPATHWAYS')
                                    & (~query_df['name'].str.contains("KEGG")) & (~query_df['name'].str.contains('HP root'))]
                query_df['Gene fraction'] = query_df['intersection_size'] / len(go_genes[g_set])
                query_df = query_df.rename({"p_value": "FDR"}, axis=1)  # According to the docu it's corrected.
                df_fetcher[g_set] = query_df
                for source in wanted_sources:
                    term_fetcher[g_set][source] = query_df[query_df['source'] == source]
            else:
                print(g_set, 'empty gene set')
        if godf_only:
            return df_fetcher
    else:
        term_fetcher = custom_dfs

    plot_go(mode='hypergeometric', wanted_sources=wanted_sources, term_fetcher=term_fetcher, keywords=keywords, cmap=cmap, fig_width=fig_width, fig_height=fig_height,
            go_genes=go_genes, dict_keys=dict_keys, numerate=numerate, font_s=font_s, rotation=rotation,
            title_tag=title_tag, legend_out=legend_out, max_terms=max_terms, out_tag=out_tag, formats=formats)


    return df_fetcher



def gsea_prerank(go_genes, weight=0, gsea_plot_out=None, out_tag='', title_tag='', max_terms='all', numerate=False, godf_only=False, translate_ensembl=False,
                 wanted_sources=['c5.hpo', 'c2.cp.wikipathways', 'c5.go.mf', 'c2.cp.reactome', 'c5.go.bp'], nes_sign='both',
                 gmt_path_pattern="/projects/abcpp/work/base_data/GSEA_gmt/human/*.v2024.1.Hs.symbols.gmt",
                 annotation='/projects/abcpp/work/base_data/gencode.v38.annotation.gtf', keywords={}, fig_width=None,
                 fig_height=None, legend_out=None, rotation=45, font_s=16, cores=20, formats=['pdf']):
    """
    Uses GSEApy prerank function to look whether gene sets from gmt files are enriched at the top or bottom of the provided lists of genes.
    CARE: The input will be automatically sorted descendingly by the first column.

    Args:
        go_genes: {class: Series or pd.FataFrame that will be sorted descendingly by the first column. The index has to be gene names or Ensembl IDs.}
        weight: How much the values should be weighted for the running-sum statistics. Set to 0 to not consider the value, which is then equal to a Kolmogorov-Smirnov statistic.
        gsea_plot_out: Whether the typical GSEA plots should be created for the significant results. GSEA then creates a folder at out_tag.
        max_terms: Allow a maximum of max_terms per dict key. Use 'all' to get all (FDR ≤ 5%). If not 'all', sort the ones with FDR ≤ 5% by the absolute NES.
        numerate: Show the size of the gene sets in parentheses on the x-axis.
        godf_only: Skip the plotting step and only return the df from GSEApy.
        translate_ensembl: Set to True if the index in the Series/DataFrame are Ensembl IDs, they will then be translated to gene names.
        wanted_sources: Which databases to plot. Will be checked against the path given in gmt_path_pattern. With the default paths vailable options are: ['c2.cp.kegg_medicus', 'c5.hpo', 'c2.cp.wikipathways', 'c5.go.mf', 'c5.go.cc', 'c2.cp.reactome', 'c5.go.bp'].
        nes_sign: Whether the enriched terms should be should be filtered for 'positive' or 'negative' NES. By Default keep 'both'.
        gmt_path_pattern: File system pattern to the gmt files.
        annotation: Path to a gtf-file, required if Ensembl IDs should be translated to symbols.
        keywords: A dictionary of {source: list of keywords}, e.g. {'c5.go.bp': ['vascular', 'angio']}, to limit the
            plot to terms containing any of the listed strings. Capitalization doesn't matter, string comparison is done
            on the lowered strings.

    Returns:
        - Returns a dict of {key: {source: df}} with a df for each gene set and GO source as provided by the GSEA package, which includes which genes matched to which terms.
    """


    keywords = {k: [t.lower() for t in vals] for k, vals in keywords.items()}  # We compare lower strings.
    gmt_files = Various.fn_patternmatch(gmt_path_pattern)
    source_tags = {"go.bp": 'GOBP', 'hpo': 'HP', 'cp.kegg_medicus': 'KEGG_MEDICUS', 'cp.reactome': 'REACTOME', 
                   'go.cc': 'GOCC', 'cp.wikipathways': 'WP', 'go.mf': 'GOMF'}

    dict_keys = list(go_genes.keys())
    term_fetcher = {c: {s: [] for s in wanted_sources} for c in dict_keys}
    df_fetcher = {c: {} for c in dict_keys}

    for g_set in dict_keys:
        if not go_genes[g_set].empty:
            g_set_df = copy.deepcopy(go_genes[g_set])

            if translate_ensembl:
                mapped_identifiers, _ = GTF_Processing.match_gene_identifiers(g_set_df.index, gtf_file=annotation, 
                                                                            scopes="symbol,alias,uniprot", 
                                                                            fields="ensembl,symbol", ensemblonly=False)
                g_set_df.index = [mapped_identifiers[g]['symbol'] for g in g_set_df.index]

            for gmt_file, source in gmt_files.items():
                if source in wanted_sources:
                    pre_res = gp.prerank(rnk=g_set_df,  # CARE sorts again descendingly by the first column.
                                        gene_sets=gmt_file,
                                        min_size=1,
                                        max_size=len(g_set_df)+1,
                                        permutation_num=1000, 
                                        outdir=gsea_plot_out,
                                        threads=cores,
                                        weight=weight,
                                        seed=1234,
                                        no_plot=not bool(gsea_plot_out),
                                        verbose=False,
                                        graph_num=len(open(gmt_file).readlines()),
                                        ).res2d
                    
                    pre_res = pre_res[pre_res['FDR q-val'] <= 0.05]
                    if nes_sign == 'positive':
                        pre_res = pre_res[pre_res['NES'] > 0]
                    elif nes_sign == 'negative':
                        pre_res = pre_res[pre_res['NES'] < 0]

                    # Strip database name.
                    if '.'.join(source.split('.')[1:]) in source_tags:
                        pre_res['Term'] = pre_res['Term'].str.replace(source_tags['.'.join(source.split('.')[1:])]+"_", '')
                    pre_res['Term'] = pre_res['Term'].str.replace('_', ' ')
                    pre_res = pre_res.sort_values(by='NES', ascending=False, key=abs)
                    pre_res['Gene fraction'] = [float(x.replace('%', '')) / 100 for x in pre_res['Gene %'].values]
                    # Rename to use the same processing as the g:profiler parsing.
                    pre_res = pre_res.rename({"Name": "mode", 'Term': 'name', 'FDR q-val': 'FDR'}, axis=1)
                    df_fetcher[g_set][source] = pre_res
                    term_fetcher[g_set][source] = pre_res
        else:
            print(g_set, 'empty gene set')
    if godf_only:
        return df_fetcher

    plot_go(mode='prerank', wanted_sources=wanted_sources, term_fetcher=term_fetcher, keywords=keywords, cmap='', fig_width=fig_width, fig_height=fig_height,
            go_genes=go_genes, dict_keys=dict_keys, numerate=numerate, font_s=font_s, rotation=rotation,
            title_tag=title_tag, legend_out=legend_out, max_terms=max_terms, out_tag=out_tag, formats=formats)

    return df_fetcher


