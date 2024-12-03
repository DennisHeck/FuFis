import os
import gzip
import argparse

"""Takes a folder with output files from the ABC scoring of STARE, and converts them into a bed-file, with the
target gene of the enhancer as csv in the last column. That means that only enhancers are written to the new file
 that also got linked to a gene in the respective scoring file. Expects the STARE files to contain 
'ABCpp_scoredInteractions', and to end with .txt.gz. Optionally another substring can be given that has to be 
contained."""

parser = argparse.ArgumentParser()
parser.add_argument("--abc_folder", required=True)
parser.add_argument("--substring", required=False, default="")
parser.add_argument("--out_folder", required=True)
args = parser.parse_args()

if not os.path.isdir(args.out_folder):
    os.mkdir(args.out_folder)

abc_files = [x for x in os.listdir(args.abc_folder) if 'ABCpp_scoredInteractions' in x and args.substring in x
             and x.endswith(".gz")]

for file in abc_files:
    print(file)
    region_gene_map = {}
    with gzip.open(args.abc_folder + '/' + file, 'rt', encoding='utf-8') as data:
        abc_head = {x: i for i, x in enumerate(data.readline().strip().replace('#', '').split('\t'))}
        for row in data.readlines():
            row = row.strip().split('\t')
            peak_str = '\t'.join([row[abc_head[c]] for c in ['chr', 'start', 'end']])
            if peak_str not in region_gene_map:
                region_gene_map[peak_str] = set()
            region_gene_map[peak_str].add(row[abc_head['Ensembl ID']])

    # Write to new file with genes as an additional csv column.
    with open(args.out_folder + '/' + file.replace('scoredInteractions', 'perEnhancer').replace('.gz', ''), "w") as output:
        for r, vals in region_gene_map.items():
            if vals:
                output.write(r + "\t" + ",".join(vals) + "\n")
            else:
                output.write(r + "\t" + "-" + "\n")  # Don't leave empty so that the column number is continuous.





