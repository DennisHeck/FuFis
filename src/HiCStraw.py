import hicstraw
import argparse
import os
import subprocess
from pathlib import Path

"""With hicstraw to write the coo-format with one file per chromosome of an hic-file or URL."""

parser = argparse.ArgumentParser()
parser.add_argument("--hic", required=True, help='Path to .hic file or an URL.')
parser.add_argument("--normalization", required=True, help="Which normalization to apply, has to be present in the hic file, will otherwise fail.")
parser.add_argument("--resolution", default=5000, type=int, help="Resolution with which to write the output.")
parser.add_argument("--out_folder", required=True, help="Folder where the files (one per chr) will be written to.")

args = parser.parse_args()
hic_source = args.hic
hic = hicstraw.HiCFile(hic_source)
print("Genome version", hic.getGenomeID())
print("Available resolutsions", hic.getResolutions())

if not os.path.isdir(args.out_folder):
  os.mkdir(args.out_folder)

for chrom in hic.getChromosomes():
    if chrom.name.replace("chr", '') not in ['M', 'All']:
        print(chrom.name, chrom.length)
        chrom_out_path = Path(args.out_folder, chrom.name + "_Contacts.txt")
        result = hicstraw.straw("observed", args.normalization, hic_source, chrom.name, chrom.name, 'BP', args.resolution)
        with open(chrom_out_path, 'w') as chr_out:
            for i in range(len(result)):
                chr_out.write("{0}\t{1}\t{2}".format(result[i].binX, result[i].binY, result[i].counts) + '\n')
        subprocess.call("gzip -f " + str(chrom_out_path), shell=True)

