import numpy as np
from pybedtools import BedTool
import subprocess
import argparse

"""
CARE: Will likely get killed on a local machine, needs considerable amount of memory.

For a scenario where one wants to have the average read counts across replicates in a unified set of peaks.
For all bam_files the read counts located fully in the bed peaks are counted, normalised for the read depth
in the bam-file (CPM), and then averaged across the bam-files. In case of one replicate, len(bam_files) is 1.
The bam-files should already be filtered in the same way as they were for constructing the peak file. E.g. if
MACS2 removed the duplicates, they should also be removed from the bam-files.
Writes a bed_file at out_path, with an appended column with the average read counts."""

parser = argparse.ArgumentParser()
parser.add_argument('--bed', required=True, help='Path to the bed-file in which regions to count.')
parser.add_argument('--bams', required=True, help='Bam-files separated by comma without whitespaces.')
parser.add_argument('--out', required=True, help='Full path to the output file which will be created.')
parser.add_argument('--cpm', default=True, type=bool, help='Whether to average for total read counts')
args = parser.parse_args()

bed_file = args.bed
bam_files = args.bams.split(',')
out_path = args.out
print('Bed:', bed_file)
print('Bams:', bam_files)
print('Out:', out_path)

bed = BedTool(''.join([x for x in open(bed_file).readlines() if not x.startswith('#') and (x[0].lower() == 'c' or x[0].isdigit())]), from_string=True)
avg_counts = {str(x).strip(): [] for x in bed}
for bam in bam_files:
    if args.cpm:
        total_reads = int(subprocess.check_output("samtools view -c " + bam, shell=True).decode("utf-8").strip())
    else:
        total_reads = 1000000  # So that we stick to the raw read counts.
    print('total reads', total_reads, bam)
    these_counts = bed.coverage(bam, F=1, counts=True)
    for count in these_counts:
        avg_counts[str('\t'.join(count.fields[:-1]))].append(int(count.fields[-1]) / total_reads * 1000000)

with open(out_path, 'w') as output:
    if open(bed_file).readline().startswith('#'):
        output.write(open(bed_file).readline())
    for peak, count in avg_counts.items():
        output.write(peak + '\t' + str(np.mean(count)) + '\n')


