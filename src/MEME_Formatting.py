from itertools import chain
import GTF_Processing

"""Collection of functions related to formatting or subsetting TF motif files in meme-format."""


def meme_monomer_map(meme_file):
    """
    Takes a motif file in meme format and creates a map of the motif name to a list of the constituent monomers,
    removing any motif versions in the list values. E.g. {'FOS(MA1951.1)': ['FOS'], 'FOXJ2::ELF1': ['FOXJ2', 'ELF1']}.
    """
    tfs = [x.split('\n\n')[0].split(' ')[0] for x in open(meme_file).read().split('MOTIF ')[1:]]
    tf_monomer_map = {t: [x.split('(')[0] for x in t.split('::')] for t in tfs}
    all_monomer_names = list(set(chain(*tf_monomer_map.values())))
    return tf_monomer_map, all_monomer_names


def meme_id_map(meme_file, gtf_file, species='human'):
    """
    Takes a meme-file and returns the list of TF gene ids belonging to that TF. Note that it assumes a certain
    syntax for the different motif versions. Species is required for looking up names missing in the
    gtf file with the MyGene.info API.
    @param species: usually 'mouse' or 'human'
    """
    tf_monomer_map, all_tf_names = meme_monomer_map(meme_file)
    tf_ids = {t: [] for t in tf_monomer_map.keys()}
    misses = []
    mapped_names, missed_names = GTF_Processing.match_gene_identifiers(all_tf_names, gtf_file=gtf_file, species=species,
                                                                       scopes="symbol", fields="ensembl")
    for tf in tf_ids:
        sub_tfs = [x.split('(')[0] for x in tf.split('::')]
        for sub in sub_tfs:
            if sub in mapped_names:
                tf_ids[tf].append(mapped_names[sub]['ensembl'])
            else:
                misses.append(tf)
                print(sub, 'name not mappable')
                if tf in tf_ids:  # Might have been already deleted as part of an earlier dimer.
                    del tf_ids[tf]

    return tf_ids, all_tf_names, misses


def subset_meme(meme_file, motif_names, out_file, include_dimers=True, exact_match=False):
    """Takes a meme file and writes a new one containing only the ones present in motif_names.
    Different motif versions are included.
    @param include_dimers: Also adds dimers containing one of the motif_names.
    @param exact_match: If the given list of names is really the original motifs' names and not just TF names."""

    meme = open(meme_file).read()
    header_block = meme.split("\nMOTIF")[0]
    tf_blocks = meme.split('MOTIF ')[1:]

    with open(out_file, 'w') as output:
        output.write(header_block + '\n')
        for block in tf_blocks:
            if not exact_match:
                this_tf = [x.split(' ')[0].split('(')[0] for x in block.split('\n\n')[0].strip().split('::')]
            else:
                this_tf = [block.split('\n\n')[0].strip().split(" ")[0]]
            if include_dimers:
                if sum([sub_t in motif_names for sub_t in this_tf]) > 0:  # If any matches.
                    output.write("MOTIF " + block)
            else:
                if sum([sub_t in motif_names for sub_t in this_tf]) == len(this_tf):  # If all match.
                    output.write("MOTIF " + block)
