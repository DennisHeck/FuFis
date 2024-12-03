import gzip
from pybedtools import BedTool
from collections import Counter
import subprocess
import os
import argparse

"""Takes a bed file with regions and runs Fimo on them. Afterwards convert Fimo's output into a matrix of regions x TFs
 with the counted binding sites."""

parser = argparse.ArgumentParser()
parser.add_argument("--bed_file", required=True)
parser.add_argument("--PWMs", required=True)
parser.add_argument("--fasta", required=True)
parser.add_argument("--fimo_src", required=True)
parser.add_argument("--out_dir", required=True)
parser.add_argument("--write_sequence", default=False, type=bool)
args = parser.parse_args()

if not args.out_dir.endswith('/'):
    args.out_dir += '/'

if not os.path.isdir(args.out_dir):
    os.mkdir(args.out_dir)

if not os.path.isfile(args.fasta+'.fai'):
    print("Genome file missing, creating with samtools faidx")
    subprocess.call('samtools faidx ' + args.fasta, shell=True)

seq_out = args.out_dir + 'Sequences.fa'
new_meme_file = args.out_dir + args.PWMs.split('/')[-1].split('.')[0] + "_bg.txt"
fimo_out = args.out_dir + 'Out_fimo.tsv.gz'

# ---------------------------------------------------------------------------------------------------
# Get regions and run Fimo
# ---------------------------------------------------------------------------------------------------

# Force the regions to be within the chromosome boundaries.
inter_bed = BedTool(args.bed_file).slop(g=args.fasta + '.fai', b=0)
inter_bed.sequence(fi=args.fasta, fo=seq_out, name=True)
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
fimo_fetcher = ''  # Create a bed object to merge hits, with the chromosome as TF#region#strand
all_regions = set()
peak_name_map = {}
if len(open(peak_file).readline().strip().split('\t')) >= 4:
    peak_name_map = {x.strip().split('\t')[3]: x.split('\t')[0]+':'+x.split('\t')[1]+'-'+x.split('\t')[2] for x in open(peak_file).readlines() if not x.startswith('#')}
with gzip.open(fimo_out, 'rt') as fimo_in:
    fimo_head = {x: i for i, x in enumerate(fimo_in.readline().strip().split('\t'))}
    for row in fimo_in:
        if not row.startswith('#'):
            entry = row.strip().split('\t')
            # if len(entry) == len(fimo_head):
            # The naming scheme of the sequences in the fasta file might differ.
            if '::' in entry[fimo_head['sequence_name']]:
                seq_region = entry[fimo_head['sequence_name']].split('::')[1]
            else:
                seq_region = peak_name_map[entry[fimo_head['sequence_name']]]
            all_regions.add(seq_region)
            fimo_fetcher += entry[fimo_head['motif_id']] + '#' + seq_region + '#' + \
                            entry[fimo_head['strand']] + '\t' + entry[fimo_head['start']] + '\t' + entry[fimo_head['stop']] + '\n'
# NOTE CARE without the s flag the merge call will ignore any strand information.
fimo_bed = BedTool(fimo_fetcher, from_string=True).sort().merge()

tfs = [x.split('\n\n')[0].split(' ')[0] for x in open(new_meme_file).read().split('MOTIF ')[1:]]
region_hits = {e: {tf: 0 for tf in tfs} for e in all_regions}

for hit in fimo_bed:
    region_hits[hit.fields[0].split('#')[1]][hit.fields[0].split('#')[0]] += 1

with open(fimo_out.replace('fimo.tsv.gz', 'Fimo_TFBS.txt'), 'w') as output:
    output.write('region\t' + '\t'.join(tfs) + '\n')
    for region, vals in region_hits.items():
        output.write(region + '\t' + '\t'.join([str(vals[tf]) for tf in tfs]) + '\n')

fimo_count_out = fimo_out.replace('fimo.tsv.gz', 'Fimo_TFBS.txt')
subprocess.call('gzip ' + fimo_count_out, shell=True)
print("FIMO hits counted")
print(fimo_count_out + '.gz')
