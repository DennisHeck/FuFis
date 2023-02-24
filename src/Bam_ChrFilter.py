import subprocess
import argparse
"""Calls different samtool commands to create a bam file with only the specified chromosomes. Also creates the index."""

parser = argparse.ArgumentParser()
parser.add_argument('--bam', required=True, help='Path to the bam-file to filter from.')
parser.add_argument('--chr', required=True, help='Which chromosome to filter for.')
args = parser.parse_args()

filter_out = args.bam.replace('.bam', '_'+args.chr+'.sam')
subprocess.call('samtools view -@ 30 -h ' + args.bam + ' ' + args.chr + ' > ' + filter_out, shell=True)
subprocess.call('samtools view -@ 30 -S -b ' + filter_out + ' > ' + filter_out.replace('.sam', '.bam'), shell=True)
subprocess.call('samtools index ' + filter_out.replace('.sam', '.bam'), shell=True)
subprocess.call('rm ' + filter_out, shell=True)
