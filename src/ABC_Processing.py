import os
import numpy as np
import gzip
import pandas as pd
from timeit import default_timer as clock
from multiprocessing import Pool


def score_summer(args):
    """Get the vector of summed activities across genes, ordered by gene_idx.
    Placed on the highest level to be parallelizable."""
    mode, run, file, gene_idx = args
    run_scores = np.zeros([len(gene_idx)])
    with gzip.open(file, 'rt') as abc_in:
        abc_head = {x: i for i, x in enumerate(abc_in.readline().strip().split('\t'))}
        for entry in abc_in:
            entry = entry.strip().split('\t')
            this_gene = entry[abc_head['Ensembl ID']].split('.')[0]
            if mode == 'AdapxC':
                run_scores[gene_idx[this_gene]] += (
                            float(entry[abc_head['adaptedActivity']]) * float(entry[abc_head['Contact']]))
            elif mode == 'SumA':
                run_scores[gene_idx[this_gene]] += float(entry[abc_head['signalValue']])
            elif mode == 'AxC':
                run_scores[gene_idx[this_gene]] += (
                        float(entry[abc_head['signalValue']]) * float(entry[abc_head['Contact']]))
    return run_scores


def abc_gene_activities(interaction_files, annotation, mode='AdapxC', cores=1, id_idx=-1):
    """Build a matrix of genes * samples with an activity value as entry, which is based on the interactions
    mapped to the gene in a sample. The activity value is the sum of adaptedActivity * contact of all enhancers mapped
    to a gene.
    @param interaction_files: Path to the ABC-files. All files containing ABCpp_scoredInteractions and are
    gzipped are used. They are identified via their column string meaning *_ABCpp_scoredInteractions_IDENTIFIER.txt.gz.
    Alternatively give a dictionary of {run_tag: ABC-file}.
    @param annotation: gtf-file that was used for the ABC-scoring. Is used to get the list of genes that can possibly
    be in the matrix. Can be gzipped or uncompressed. Version suffixes are removed.
    @param mode: Can be SumA for summing the signalValue column, AdapxC for summing the product of adaptedActivity and
    contact, or AxC for summing the product of signalValue and contact.
    @param id_idx: Which entry of filename.split('_') to take as identifier.
    @return gene_scores_df: Pandas Dataframe of genes * samples including Ensembl IDs as index and the column string
    of the ABC runs as columns. Note that there is no normalization done, a direct comparison between runs is not
    advised, unless the activity and contact measurements were comparable."""
    start = clock()
    if mode not in ['AdapxC', 'SumA', 'AxC']:
        print("ERROR: unknown mode")
        return
    # First get the vector of possible genes.
    if annotation.endswith('.gz'):
        anno_opener = gzip.open(annotation, 'rt')
    else:
        anno_opener = open(annotation, 'r')
    genes = list(set([x.split('\t')[8].split('gene_id "')[-1].split('"; ')[0].split('.')[0]
                      for x in anno_opener.readlines() if not x.startswith('#') and x.split('\t')[2] == 'gene']))
    gene_idx = {g: i for i, g in enumerate(genes)}

    # Now get the samples/ABC runs for the other dimension.
    if not type(interaction_files) == dict:
        abc_runs = {x.split('.txt')[0].split('_')[id_idx]: interaction_files + '/' + x
                    for x in os.listdir(interaction_files) if "ABCpp_scoredInteractions" in x and x.endswith('.gz')}
    else:
        abc_runs = interaction_files

    process_pool = Pool(processes=cores)
    pool_summer = process_pool.map(score_summer, [[mode, run, file, gene_idx] for run, file in abc_runs.items()])
    process_pool.close()

    # Convert to a dataframe to add the gene IDs and run tags again.
    gene_scores_df = pd.DataFrame(pool_summer, columns=gene_idx.keys(), index=abc_runs.keys()).T
    print("ABC scores summed", clock() - start)

    return gene_scores_df


def interaction_fetcher_helper(args):
    """
    File reader to fetch ABC interactions. In a separate function to enable parallelization.
    """
    file, file_tag, cutoff = args
    file_interactions = set()
    with gzip.open(file, 'rt') as abc_file:
        head = {x: i for i, x in enumerate(abc_file.readline().strip().split('\t'))}
        for entry in abc_file:
            entry = entry.strip().split('\t')
            interaction = entry[head['#chr']] + '-' + entry[head['start']] + '-' + entry[head['end']] + '#' + \
                          entry[head['Ensembl ID']].split('.')[0]
            score = float(entry[head['ABC-Score']])
            if score >= cutoff:
                file_interactions.add(interaction)
    return [file_tag, file_interactions]


def interaction_fetcher(pattern, pooled=True, cutoff=0.02, n_cores=1):
    """
    Takes the path to a ABC-scoring file or a pattern to multiple ABC-scoring files to get the gABC interactions surpassing the threshold in the
    format chr-start-end#EnsemblID. The version from the Ensembl ID is cut.

    Args:
        pattern: File path pattern e.g., ABC/ABCpp_scoredInteractions_*.bed, or the path to just one individual file.
        pooled: If True return one set of interactions pooled across all files. Otherwise create a dictionary with an entry for each wildcard matching the pattern.
        cutoff: Threshold when to fetch a gABC interaction.
    """
    start = clock()
    if '*' in pattern:
        abc_files = Various.fn_patternmatch(pattern)
    else:  # Assume we have just one file.
        abc_files = {pattern: pattern.split('/')[-1]}

    process_pool = Pool(processes=n_cores)
    pool_fetcher = process_pool.map(interaction_fetcher_helper,
                                    [[file, file_tag, cutoff] for file, file_tag in abc_files.items()])
    process_pool.close()

    if pooled:
        inter_dict = set.union(*[x[1] for x in pool_fetcher])
    else:
        inter_dict = {x[0]: x[1] for x in pool_fetcher}

    print(clock() - start, 'interactions fetched')
    return inter_dict

