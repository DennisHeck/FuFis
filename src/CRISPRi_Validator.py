import pandas as pd
from timeit import default_timer as clock
import pybedtools
import gzip
import copy
import numpy as np
import os
import json
from datetime import datetime
from multiprocessing import Pool
import GenomeLifter
import GTF_Processing


"""
Takes a list of hg19 enhancer-gene interactions with a score and intersects them with experimentally validated
interactions from three CRISPRi screens. Returns a df, which can be used to plot PR curves and the sorts.
The formatting of the validation screen result files is included here to enhance reproducibility. Especially 
Fulco has some special catches.
"""


def gtf_information():
    hg19_annotation = "/projects/abcpp/work/base_data/gencode.v19.annotation.gtf"
    hg38_annotation = '/projects/abcpp/work/base_data/gencode.v38.annotation.gtf'

    if hg19_annotation.endswith('.gz'):  # For Gschwind we already mapped the Ensembl IDs prior.
        hg19_opener = gzip.open(hg19_annotation, 'rt')
    else:
        hg19_opener = open(hg19_annotation)
    gene_name_map = {}
    for entry in hg19_opener:
        if not entry.startswith('#') and entry.split('\t')[2] == 'gene':
            gene_name_map[entry.split('\t')[8].split('gene_name "')[-1].split('"; ')[0]] = entry.split('\t')[8].split('gene_id "')[-1].split('"; ')[0].split('.')[0]

    hg19_tss_dict = GTF_Processing.gene_window_bed(gtf_file=hg19_annotation, extend=1, tss_type='5', dict_only=True)
    hg38_tss_dict = GTF_Processing.gene_window_bed(gtf_file=hg38_annotation, extend=1, tss_type='5', dict_only=True)
    return gene_name_map, hg19_tss_dict, hg38_tss_dict


def match_scores(args):
    """
    Takes an ABC-scoring file and validated interactions to fetch the matching scores. Iterate over multiple
    screens if present.
    """
    # abc_file, score_tag, validation_beds, validation_sets, match_cols=['ABC-Score'], gene_col='Ensembl ID', genome_version='hg19',
    #                      inter_f=1e-09, inter_F=1e-09, inter_e=False = args
    start_f = clock()
    score_tag, file_info, validation_beds, validation_sets, gene_name_map = args
    print(score_tag)

    fetch_cols = file_info['score']
    if type(fetch_cols) != list:
        fetch_cols = list(fetch_cols)
    if 'avg_cols' in file_info:
        fetch_cols += file_info['avg_cols']

    abc_rows = pd.read_table(file_info['file'], sep='\t', header=0)
    abc_head = {x : i for i, x in enumerate(abc_rows.columns)}
    abc_rows = abc_rows.astype({"start": int, "end": int}).astype({"#chr": str, 'start': str, 'end': str}).values
    print('file read', clock() - start_f)
    chr_prefix = "chr" if not str(abc_rows[0][0]).startswith('chr') else ''

    if file_info['genome_version'] == 'hg19':
        print('Lifting to hg38', score_tag)
        abc_rows = GenomeLifter.genome_lifter(abc_rows, input_version='hg19', output_version='hg38')
    peaks_bed = pybedtools.BedTool('\n'.join(set([chr_prefix + '\t'.join(x[:3]) for x in abc_rows])),
                                    from_string=True)
    
    matched_scores = {}
    for screen in validation_beds:
        # First find the enhancers that intersect a validated region at all.
        hit_peaks = set([str(x).strip() for x in peaks_bed.intersect(validation_beds[screen], u=True, f=file_info['inter_f'], F=file_info['inter_F'], e=file_info['inter_e'])])
        peak_preturb_map = {x: set() for x in hit_peaks}
        for inter in peaks_bed.intersect(validation_beds[screen], wo=True, f=file_info['inter_f'], F=file_info['inter_F'], e=file_info['inter_e']):
            peak_preturb_map['\t'.join(inter.fields[:3])].add('\t'.join(inter.fields[3:6]))
        
        gene_misses = set()
        matched_scores[screen] = {i: {c: [] for c in [score_tag+' '+x for x in fetch_cols]} for i in validation_sets[screen]}
        # Now go through the queried interactions again to find peaks that intersect a perturbed site.
        for inter in abc_rows:
            peak_str = chr_prefix + '\t'.join(inter[:3])
            if peak_str in hit_peaks:
                inter_gene = inter[abc_head[file_info['gene_col']]].split('.')[0]
                if not inter_gene.startswith('ENSG'):
                    try:
                        inter_gene = gene_name_map[inter_gene]
                    except KeyError:
                        gene_misses.add(inter_gene)
                        continue
                for target in peak_preturb_map[peak_str]:
                    if target + '\t' + inter_gene in matched_scores[screen]:
                        for match_col in fetch_cols:
                            matched_scores[screen][target + '\t' + inter_gene][score_tag + ' ' + match_col].append(
                                float(inter[abc_head[match_col]]))
        if gene_misses:
            print(screen, score_tag, 'missed genes', len(gene_misses), gene_misses)
    print(score_tag, 'processed', clock() - start_f)
    return matched_scores
                

def get_crispris_from_joint(to_validate, joint_file, screens=['gasperini', 'tap', 'fulco', 'gschwind'],
                            all_required=True, distance=False, meta_out='', meta_tag='', n_cores=1, tmp_dir=''):
    """
    Very similar to the other function below, and admittedly ugly, but we don't need all the screen-specific 
    catches in this one here. And it's all in hg38.
    Takes a dictionary with files of interactions, which should be bed-styled with a header holding the column names
    for the gene and a score column (e.g. ABC-Score). The genes may be symbols or IDs. Will return a df for each
    CRISPRi screen and the score of the intersecting interactions from the file. If multiple enhancer intersect a
    validated interaction the sum is taken.
    :param to_validate: {label: {
                       'file': file_path,
                       'score': column name which score to take, can be an iterable with multiple columns,
                       'gene_col': column name to find the gene,
                       'avg_cols': a list of columns to add to the df, will be averaged instead of summed.}}
    :param screens: Names screens to validate on, can be used to reduce which CRISPRi data sets to include.
    :param all_required: If all score columns must have an intersection with the validated interactions, if False
    sets the score to None of the ones missing a score.
    :param meta_out: Directory where to write the JSON file to with the metadata of the function call.
    :param meta_tag: Prefix for the metadata file. The filename will also include the distance and all_required info.
    :param n_cores: Parallelized over the files in to_validate, as read-in is usually the most time-consuming step.
    :param tmp_dir: Optional dictionary for temporary BedTool files.
    :return: A dictionary for each screen with a df with the validated interactions and the score from the given
    interaction file, also has entries for which columns were retrieved.
    """
    start = clock()
    gene_name_map, _, _ = gtf_information()
    pybedtools.helpers.set_tempdir("/home/dhecker/tmpdir/")

    # Store where to find the files and which columns to grab. The order must be the same, it relies on the indices.
    validation_info = pd.read_table(joint_file, sep='\t', header=0)

    # We make a copy to be able to streamline some entries.
    score_files = copy.deepcopy(to_validate)
    validation_dfs = {v: {'df': None, 'score_cols': None, 'avg_cols': None} for v in screens}
    score_cols = []
    average_cols = []  # The columns from 'avg_cols' will be averaged instead of summed like the 'score' column.
    # Add default values to the score_files dictionary.
    for score in score_files:
        if 'gene_col' not in score_files[score]:
            score_files[score]['gene_col'] = 'Ensembl ID'
        if 'score' not in score_files[score]:
            score_files[score]['score'] = ['ABC-Score']
        if 'genome_version' not in score_files[score]:
            score_files[score]['genome_version'] = 'hg19'
        if 'aggregate_mode' not in score_files[score]:
            score_files[score]['aggregate_mode'] = 'sum'
        if 'inter_f' not in score_files[score]:
            score_files[score]['inter_f'] = 1e-09
        if 'inter_F' not in score_files[score]:
            score_files[score]['inter_F'] = 1e-09
        if 'inter_e' not in score_files[score]:
            score_files[score]['inter_e'] = False
    for score, vals in score_files.items():
        if not type(vals['score']) == list:
            score_files[score]['score'] = [score_files[score]['score']]
        score_cols += [score + ' ' + c for c in vals['score']]
        if 'avg_cols' in vals:
            average_cols += [score + ' ' + c for c in vals['avg_cols']]

    # After filling the defaults write a metadata file to trace the function call.
    json_dict = {'Date': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                 "to_validate": to_validate,
                 'joint_file': joint_file,
                 'screens': screens,
                 'all_required': all_required,
                 'distance': distance,
                 'meta_out': meta_out,
                 'meta_tag': meta_tag}
    json.dump(json_dict, open(meta_out+'/'+meta_tag+"_reqAll"+str(all_required)+"_Dist"+str(distance)+"_Metadata.json", 'w'))

    validation_beds = {}
    validation_sets = {}
    validation_dicts = {}
    for screen in screens:
        screen_df = validation_info[validation_info['Screen'] == screen]
        if distance:
            screen_df = screen_df[screen_df['Distance'] <= distance]
        validation_scores = {
            '\t'.join([str(x[c]) for c in ['#chr', 'start', 'end', 'Ensembl ID']]): {'effect': float(x['effect']),
                                                                         'p-adj': float(x['p-adj']),
                                                                         'Significant': x['Significant']} for x in
            screen_df.to_dict(orient='records')}
        for entry in validation_scores:
            validation_scores[entry].update({c: [] for c in score_cols + average_cols})
        validation_dicts[screen] = validation_scores
        validation_beds[screen] = pybedtools.BedTool('\n'.join(validation_scores.keys()), from_string=True)
        validation_sets[screen] = set(validation_scores.keys())

    process_pool = Pool(processes=n_cores)
    pool_fetcher = process_pool.map(match_scores, [[s_tag, vals, validation_beds, validation_sets, gene_name_map] for s_tag, vals in score_files.items()])
    
    for mapped_scores in pool_fetcher:
        for screen in screens:
            for inter in validation_dicts[screen]:
                validation_dicts[screen][inter].update(mapped_scores[screen][inter])

    score_cols_agg = {}
    for s_tag, values in score_files.items():   
        for s_col in values['score']:
            score_cols_agg[s_tag + ' ' + s_col] = score_files[s_tag]['aggregate_mode']


    for screen in screens:  # This is admittedly very ugly.
        all_hits = []
        not_all = []
        for k, vals in validation_dicts[screen].items():
            if (all_required and np.all([1 if vals[score] else 0 for score in score_cols])) or \
                (not all_required and np.any([1 if vals[score] else 0 for score in score_cols])):
                inter_values = k.split('\t') + [vals['effect'], vals['p-adj'], vals['Significant']]
                inter_values += [(sum(vals[col]) if score_cols_agg[col] == 'sum' else np.mean(vals[col])) 
                                 if vals[ col] else None for col in score_cols]
                inter_values += [np.mean(vals[col]) if vals[col] else 0 for col in average_cols]
                all_hits.append(inter_values)
            else:
                not_all.append([k, vals])
        print("Validated interactions without any match", len(not_all))

        compare_cols = ['chr', 'start', 'end', 'gene'] + list(list(validation_dicts[screen].values())[0].keys())
        col_types = ['str', 'int64', 'int64', 'str', 'float64', 'float64', bool] + ['float64'] * len(
            score_cols + average_cols)
        compare_abcpp_df = pd.DataFrame(all_hits, columns=compare_cols)
        compare_abcpp_df = compare_abcpp_df.astype({c: col_types[i] for i, c in enumerate(compare_cols)})

        validation_dfs[screen]['df'] = compare_abcpp_df
        validation_dfs[screen]['score_cols'] = score_cols
        validation_dfs[screen]['avg_cols'] = average_cols

    pybedtools.helpers.cleanup(verbose=False, remove_all=False)
    
    print('Total', clock() - start)
    return validation_dfs


def get_crispris(to_validate, screen_folder, screens=['gasperini', 'tap', 'fulco', 'gschwind'],
                 all_required=True, distance=False):
    """
    Takes a dictionary with files of interactions, which should be bed-styled with a header holding the column names
    for the gene and a score column (e.g. ABC-Score). The genes may be symbols or IDs. Will return a df for each
    CRISPRi screen and the score of the intersecting interactions from the file. If multiple enhancer intersect a
    validated interaction the sum is taken.
    :param to_validate: {label: {
                       'file': file_path,
                       'score': column name which score to take, can be an iterable with multiple columns,
                       'gene': column name to find the gene,
                       'avg_cols': a list of columns to add to the df, will be averaged instead of summed.}}
    :param screen_folder: Path to the folder with the CRISPRi result files.
    :param screens: Names screens to validate on, can be used to reduce which CRISPRi data sets to include.
    :param all_required: If all score columns must have an intersection with the validated interactions, if False
    sets the score to 0 of the ones missing a score.
    :return: A dictionary for each screen with a df with the validated interactions and the score from the given
    interaction file, also has entries for which columns were retrieved.
    """

    pybedtools.helpers.set_tempdir("/home/dhecker/tmpdir/")
    
    # Store where to find the files and which columns to grab. The order must be the same, it relies on the indices.
    validation_info = {'gasperini': {'validation_file': screen_folder+'/GSE120861_all_deg_results.at_scale.txt',
                                     'columns': {'chr': 'target_site.chr', 'start': 'target_site.start', 'end': 'target_site.stop', 'gene': 'ENSG',
                                                'effect': 'fold_change.transcript_remaining', 'p-adj': 'pvalue.empirical.adjusted'},
                                     'genome_version': 'hg19'
                                     },
                       # 'review_gasp': {'validation_file': '/Users/dennis/Desktop/PaperReviews/ML_EP_BioInf23/Gasperini2019.at_scale.ABC.TF.erole.grouped.train.txt',
                       #               'columns': {'chr': 'chrEnhancer', 'start': 'startEnhancer',
                       #                           'end': 'endEnhancer', 'gene': 'EnsemblD',
                       #                           'effect': 'pValueAdjusted',
                       #                           'p-adj': 'pValueAdjusted'},
                       #               },
                       'tap': {'validation_file': screen_folder+"/TAP_5MB1Sig.txt",
                               'columns': {'chr': 'enh_chr', 'start': 'enh_start', 'end': 'enh_end', 'gene': 'gene',
                                           'effect': 'manual_lfc', 'p-adj': 'pval_adj_allTests'},
                               'genome_version': 'hg19'
                               },
                       'fulco': {'validation_file': screen_folder+'/41588_2019_538_MOESM3_ESM.xlsx',
                               'columns': {'chr': 'chr', 'start': 'start', 'end': 'end', 'gene': 'Gene',
                                           'effect': 'Fraction change in gene expr', 'p-adj': 'Adjusted p-value',
                                           'Significant': 'Significant'},
                                 'genome_version': 'hg19'
                               },
                       'gschwind': {'validation_file': screen_folder+'/EPCrisprBenchmark_ensemble_data_GRCh38_wEnsembl.tsv',
                               'columns': {'chr': 'chrom', 'start': 'chromStart', 'end': 'chromEnd', 'gene': 'Ensembl ID',
                                           'effect': 'EffectSize', 'p-adj': 'pValueAdjusted',
                                           'Significant': 'Significant'},
                                    'genome_version': 'hg38'}
                       }

    # We make a copy to be able to streamline some entries.
    score_files = copy.deepcopy(to_validate)
    validation_dfs = {v: {'df': None, 'score_cols': None, 'avg_cols': None} for v in validation_info}
    for mode in screens:
        print(mode)
        if mode == 'fulco' and not np.any(["GraphReg" in x for x in score_files]):  # Praise Excel.
            fulco_df = pd.read_excel(validation_info[mode]['validation_file'], sheet_name='Supplementary Table 6a', header=1, engine='openpyxl')
            # It's PVT1 in the gene annotation, but PVT1-TSS1 in 6a.
            fulco_df.loc[fulco_df.Gene == 'PVT1-TSS1', 'Gene'] = "PVT1"
            fulco_df['Gene'] = [gene_name_map[g] for g in fulco_df['Gene'].values]
            validation_rows = [{c: str(x[val_c]) for c, val_c in validation_info[mode]['columns'].items()} for x in fulco_df.to_dict(orient='records')]
            validation_scores = {'\t'.join([x[c] for c in ['chr', 'start', 'end', 'gene']]): {'effect': float(x['effect']), 'p-adj': float(x['p-adj']), 'Significant': x['Significant'] == 'True'} for x in validation_rows}
        else:
            head = {x: i for i, x in enumerate(open(validation_info[mode]['validation_file']).readline().strip().split('\t'))}
            validation_rows = [{c: x.split('\t')[head[val_c]] for c, val_c in validation_info[mode]['columns'].items()} for x in open(validation_info[mode]['validation_file']).read().strip().split('\n')[1:] if 'not_applicable' not in x]

            if not validation_rows[0]['gene'].startswith('ENS'):  # Already convert symbols to Ensembl IDs.
                for n, row in enumerate(validation_rows):
                    row['gene'] = gene_name_map[row['gene']]
                    validation_rows[n] = row
            validation_rows = [x for x in validation_rows if x['p-adj']]  # Bless the person who just left an empty field.
            if mode == 'gschwind':
                validation_scores = {
                    '\t'.join([x[c] for c in ['chr', 'start', 'end', 'gene']]): {'effect': float(x['effect']), 'p-adj': float(x['p-adj']), 'Significant': x['Significant'] == 'True'} for x in validation_rows}
            else:
                validation_scores = {'\t'.join([x[c] for c in ['chr', 'start', 'end', 'gene']]): {'effect': float(x['effect']), 'p-adj': float(x['p-adj']), 'Significant': False} for x in validation_rows}
        if distance:
            if validation_info[mode]['genome_version'] == 'hg19':
                this_tss_dict = hg19_tss_dict
            else:
                this_tss_dict = hg38_tss_dict
            validation_scores = {k: val for k, val in validation_scores.items() if k.split('\t')[3] in this_tss_dict and min([abs(int(k.split('\t')[1]) - next(iter(this_tss_dict[k.split('\t')[3]]['tss']))),
                                                                                        abs(int(k.split('\t')[2]) - next(iter(this_tss_dict[k.split('\t')[3]]['tss'])))]) <= distance}
        validation_bed = pybedtools.BedTool('\n'.join(validation_scoreskeys()), from_string=True)

        score_cols = []
        average_cols = []  # The columns from 'avg_cols' will be averaged instead of summed like the 'score' column.
        # Add default values to the score_files dictionary.
        for score in score_files:
            if 'gene' not in score_files[score]:
                score_files[score]['gene'] = 'Ensembl ID'
            if 'score' not in score_files[score]:
                score_files[score]['score'] = 'ABC-Score'
            if 'genome_version' not in score_files[score]:
                score_files[score]['genome_version'] = 'hg19'
            if 'aggregate_mode' not in score_files[score]:
                score_files[score]['aggregate_mode'] = 'sum'
            if 'inter_f' not in score_files[score]:
                score_files[score]['inter_f'] = 1e-09
            if 'inter_F' not in score_files[score]:
                score_files[score]['inter_F'] = 1e-09
            if 'inter_e' not in score_files[score]:
                score_files[score]['inter_e'] = False
        for score, vals in score_files.items():
            if not type(vals['score']) == list:
                score_files[score]['score'] = [score_files[score]['score']]
            score_cols += [score + ' ' + c for c in vals['score']]
            # score_cols += vals['score']
            if 'avg_cols' in vals:
                average_cols += [score + ' ' + c for c in vals['avg_cols']]
        for entry in validation_scores:
            validation_scores[entry].update({c: [] for c in score_cols + average_cols})

        def match_scores(abc_file, score_tag, score_col='ABC-Score', gene_col='Ensembl ID', genome_version='hg19',
                         inter_f=1e-09, inter_F=1e-09, inter_e=False):
            """Takes an ABC-scoring file and validated interactions to fetch the matching scores."""

            if os.path.isdir(abc_file):
                abc_head, abc_rows = folder_reader(abc_file, score_tag)
            else:
                if abc_file.endswith('.gz'):
                    with gzip.open(abc_file, 'rt') as abc_in:
                        abc_head = {x: i for i, x in enumerate(abc_in.readline().replace('#', '').strip().split('\t'))}
                        abc_rows = abc_in.readlines()
                else:
                    abc_head = {x: i for i, x in enumerate(open(abc_file).readline().replace('#', '').strip().split('\t'))}
                    abc_rows = open(abc_file).readlines()[1:]
            abc_rows = [x.strip().split('\t') for x in abc_rows]
            # print(abc_rows[:3])
            abc_rows = [[str(x[0]), str(int(float(x[1]))), str(int(float(x[2])))] + x[3:] for x in abc_rows]
            # print(abc_rows[:3])
            # First find the enhancers that intersect a validated region at all.
            chr_prefix = "chr" if not abc_rows[0][0].startswith('chr') else ''

            if genome_version == 'hg38' and validation_info[mode]['genome_version'] == 'hg19':
                print('Lifting to hg19', score_tag)
                abc_rows = GenomeLifter.genome_lifter(abc_rows, input_version='hg38', output_version='hg19')
            elif genome_version == 'hg19' and validation_info[mode]['genome_version'] == 'hg38':
                print('Lifting to hg38', score_tag)
                abc_rows = GenomeLifter.genome_lifter(abc_rows, input_version='hg19', output_version='hg38')

            peaks_bed = pybedtools.BedTool('\n'.join(set([chr_prefix+'\t'.join(x[:3]) for x in abc_rows])), from_string=True)
            hit_peaks = set([str(x).strip() for x in peaks_bed.intersect(validation_bed, u=True, f=inter_f, F=inter_F, e=inter_e)])
            peak_preturb_map = {x: set() for x in hit_peaks}
            gene_misses = set()
            for inter in str(peaks_bed.intersect(validation_bed, wo=True, f=inter_f, F=inter_F, e=inter_e)).strip().split('\n'):
                if '\t'.join(inter.split('\t')[:3]):  # For whatever reason there can be an empty intersection.
                    peak_preturb_map['\t'.join(inter.split('\t')[:3])].add('\t'.join(inter.split('\t')[3:6]))

            # Now go through the queried interactions again to find peaks that intersect a perturbed site.
            for inter in abc_rows:
                peak_str = chr_prefix + '\t'.join(inter[:3])
                if peak_str in hit_peaks:
                    inter_gene = inter[abc_head[gene_col]].split('.')[0]
                    if not inter_gene.startswith('ENSG'):
                        try:
                            inter_gene = gene_name_map[inter_gene]
                        except KeyError:
                            gene_misses.add(inter_gene)
                            continue
                    for target in peak_preturb_map[peak_str]:
                        if target + '\t' + inter_gene in validation_scores:
                            validation_scores[target + '\t' + inter_gene][score_tag + ' ' + score_col].append(
                                float(inter[abc_head[score_col]]))
            print(score_tag, score_col, 'missed genes', len(gene_misses), gene_misses)


        score_cols_agg = {}
        for s_tag, values in score_files.items():
            for s_col in values['score']:
                match_scores(values['file'], score_tag=s_tag, score_col=s_col, gene_col=values['gene'],
                             genome_version=values['genome_version'], inter_f=values['inter_f'],
                             inter_F=values['inter_F'], inter_e=values['inter_e'])
                score_cols_agg[s_tag + ' ' + s_col] = score_files[s_tag]['aggregate_mode']
            if 'avg_cols' in values:
                for col in values['avg_cols']:
                    match_scores(values['file'], score_tag=s_tag, score_col=col, gene_col=values['gene'],
                                 genome_version=values['genome_version'], inter_f=values['inter_f'],
                                 inter_F=values['inter_F'], inter_e=values['inter_e'])

        compare_cols = ['chr', 'start', 'end', 'gene'] + list(list(validation_scores.values())[0].keys())
        col_types = ['str', 'int64', 'int64', 'str', 'float64', 'float64', bool] + ['float64']*len(score_cols+average_cols)
        all_hits = []
        not_all = []
        for k, vals in validation_scores.items():
            if (all_required and np.all([1 if vals[score] else 0 for score in score_cols])) or \
                    (not all_required and np.any([1 if vals[score] else 0 for score in score_cols])):
                all_hits.append(k.split('\t') + [vals['effect'], vals['p-adj'], vals['Significant']] +
                                [(sum(vals[col]) if score_cols_agg[col]=='sum' else np.mean(vals[col])) if vals[col] else None for col in score_cols] +  # TODO
                                [np.mean(vals[col]) if vals[col] else 0 for col in average_cols])
            else:
                not_all.append([k, vals])
        print("Interactions with only partial columns", len(not_all))

        compare_abcpp_df = pd.DataFrame(all_hits, columns=compare_cols)
        compare_abcpp_df = compare_abcpp_df.astype({c: col_types[i] for i, c in enumerate(compare_cols)})
        if (mode != 'fulco' and mode != 'gschwind') or np.any(["GraphReg" in x for x in score_files]):
            compare_abcpp_df['Significant'] = [True if x <= 0.05 else False for x in compare_abcpp_df['p-adj'].values]
        compare_abcpp_df[validation_info[mode]['columns']['effect']] = compare_abcpp_df['effect']

        validation_dfs[mode]['df'] = compare_abcpp_df
        validation_dfs[mode]['score_cols'] = score_cols
        validation_dfs[mode]['avg_cols'] = average_cols

        pybedtools.helpers.cleanup(verbose=False, remove_all=False)

    return validation_dfs


