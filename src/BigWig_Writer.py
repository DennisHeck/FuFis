import subprocess
import os
import argparse

"""Writes a CPM-normalised bigwig-file for each bam in a folder."""

parser = argparse.ArgumentParser()
parser.add_argument('--bam_folder', required=True, help="Folder with the bam-files")
parser.add_argument('--out_folder', required=True, help="Output folder, will be created if not existent")
parser.add_argument("--effectiveGenomeSize", required=True, help="Genome size, see https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html")
args = parser.parse_args()

bam_files = [args.bam_folder + '/' + x for x in os.listdir(args.bam_folder) if x.endswith('.bam')]

if not os.path.isdir(args.out_folder):
    os.mkdir(args.out_folder)

for file in bam_files:
    print(file)
    # If the index file is missing we copy the bam file first, to avoid permission issues.
    if not os.path.isfile(file + '.bai'):
        print('indexing')
        subprocess.call('samtools index -@ 20 ' + file, shell=True)
    print('bamCoverage -b ' + file + ' -p 20 --binSize 20 --normalizeUsing CPM --effectiveGenomeSize '
                    + args.effectiveGenomeSize + ' -o ' + args.out_folder + '/' +
                    file.split('/')[-1].replace('.bam', '.bw'))#, shell=True)

