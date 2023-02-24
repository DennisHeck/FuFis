import gzip
from pybedtools import BedTool
from collections import Counter
import subprocess
import os
import argparse
import TSS_Fetcher

"""For each gene call TFBS in a promoter window at the 5'TSS with Fimo. Then convert Fimo's output into a matrix
of genes x TFs with the counted binding sites."""

parser = argparse.ArgumentParser()
parser.add_argument('--organism', required=True)
parser.add_argument('--version', required=True)
parser.add_argument("--gtf", required=True)
parser.add_argument("--PWMs", required=True)
parser.add_argument("--fasta", required=True)
parser.add_argument("--fimo_src", required=True)
parser.add_argument("--out_dir", required=True)
parser.add_argument("--extend", default=1000, type=int, help="Number of bp to extend the promoter in each direction.")
parser.add_argument("--open_regions", default=False, type=str, help="Optional bed-file to restrict the promoter regions to.")
parser.add_argument("--write_sequence", default=False, type=bool)
args = parser.parse_args()

if not args.out_dir.endswith('/'):
    args.out_dir += '/'

if not os.path.isdir(args.out_dir):
    os.mkdir(args.out_dir)

if not os.path.isfile(args.fasta+'.fai'):
    print("Genome file missing, creating with samtools faidx")
    subprocess.call('samtools faidx ' + args.fasta + ' -o', shell=True)

seq_out = args.out_dir + args.organism + '_' + args.version + '_PromoterSequences_'+str(args.extend*2//1000)+'kb.fa'
new_meme_file = args.out_dir + args.PWMs.split('/')[-1].split('.')[0] + "_" + str(args.extend*2//1000)+'kb_bg.txt'
fimo_out = args.out_dir + args.organism + '_' + args.version + '_Promoter_'+str(args.extend*2//1000)+'kb_fimo.tsv.gz'

# ---------------------------------------------------------------------------------------------------
# Get Promoter and run Fimo
# ---------------------------------------------------------------------------------------------------
if args.gtf.endswith('.gz'):
    gtf_rows = [x.strip().split('\t') for x in gzip.open(args.gtf, 'rt').readlines() if not x.startswith('#')]
else:
    gtf_rows = [x.strip().split('\t') for x in open(args.gtf).readlines() if not x.startswith('#')]
genes = set([x[8].split('gene_id "')[-1].split('";')[0].split('.')[0] for x in gtf_rows])

promoter_bed = TSS_Fetcher.gene_window_bed(args.gtf, args.extend, gene_set=set(), dict_only=False, open_regions=args.open_regions)
# Force the regions to be within the chromosome boundaries.
promoter_bed = promoter_bed.slop(g=args.fasta+'.fai', b=0)
promoter_bed.sequence(fi=args.fasta, fo=seq_out, name=True)
print('sequence fasta stored at', seq_out)

# Now get the base content of the written fasta-file and create an adapted meme file with the background.
base_occs = Counter(''.join([x.strip().lower().replace('n', '') for x in open(seq_out).readlines() if not x.startswith('>')]))
cg_content = (base_occs['c'] + base_occs['g']) / sum(base_occs.values()) / 2
at_content = (base_occs['a'] + base_occs['t']) / sum(base_occs.values()) / 2

meme = open(args.PWMs).read()
meme = [x for x in meme.split('\n\n') if not x.startswith('Background')]
new_background = 'Background letter frequencies:\nA ' + str(round(at_content, 5)) + ' C ' + str(round(
    cg_content, 5)) + ' G ' + str(round(cg_content, 5)) + ' T ' + str(round(at_content, 5))
new_meme = meme[:3] + [new_background] + meme[3:]
open(new_meme_file, 'w').write('\n\n'.join(new_meme))
print('new meme file with matching background at', new_meme_file)

print("Running Fimo")
write_seq = '' if args.write_sequence else '--skip-matched-sequence'
out_suffix = '' if args.write_sequence else " > " + fimo_out.replace('.gz', '')  # Fimo has its own names in this case.

bashCommand = args.fimo_src + " --thresh 0.0001 " + write_seq + " --verbosity 1 --bfile --motif-- " + new_meme_file \
              + " " + seq_out + " > " + fimo_out.replace('.gz', '')
subprocess.call(bashCommand, shell=True)

if args.write_sequence:
    os.rename(args.out_dir + '/' + 'fimo.tsv', fimo_out.replace('.gz', ''))
subprocess.call('gzip ' + fimo_out.replace('.gz', ''), shell=True)
print('fimo output at', fimo_out)


# ---------------------------------------------------------------------------------------------------
# Process Fimo output
# ---------------------------------------------------------------------------------------------------
print("Processing Fimo output")
fimo_head = {x: i for i, x in enumerate(gzip.open(fimo_out, 'rt').readline().strip().split('\t'))}
fimo_fetcher = ''  # Create a bed object to merge hits, with the chromosome as TF#gene#strand
for entry in [x.strip('\n').split('\t') for x in gzip.open(fimo_out, 'rt').readlines()[1:] if not x.startswith('#')]:
    if len(entry) == len(fimo_head):
        fimo_fetcher += entry[fimo_head['motif_id']] + '#' + entry[fimo_head['sequence_name']].split('::')[0] + '#' + \
                        entry[fimo_head['strand']] + '\t' + entry[fimo_head['start']] + '\t' + entry[fimo_head['stop']] + '\n'
fimo_bed = BedTool(fimo_fetcher, from_string=True).sort().merge()

tfs = [x.split('\n\n')[0].split(' ')[0] for x in open(new_meme_file).read().split('MOTIF ')[1:]]
gene_hits = {g: {tf: 0 for tf in tfs} for g in genes}

for hit in fimo_bed:
    gene_hits[hit.fields[0].split('#')[1]][hit.fields[0].split('#')[0]] += 1

with open(fimo_out.replace('fimo.tsv.gz', 'Fimo_TFBS.txt'), 'w') as output:
    output.write('Ensembl ID\t' + '\t'.join(tfs) + '\n')
    for gene, vals in gene_hits.items():
        output.write(gene + '\t' + '\t'.join([str(vals[tf]) for tf in tfs]) + '\n')

subprocess.call('gzip ' + fimo_out.replace('fimo.tsv.gz', 'Fimo_TFBS.txt'), shell=True)


