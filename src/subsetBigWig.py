#!/usr/bin/env python
# Source https://gist.github.com/dpryan79/740f2d00ce6b509ab9644fc43418996c
import pyBigWig
import argparse
import os.path
import sys

parser = argparse.ArgumentParser(description="Subset a single chromosome from a bigWig file.")
parser.add_argument("--headerToo", action="store_true", help="Subset the header too, which can expedite some downstream programs but possibly cause problems for others.")
parser.add_argument("input", help="Input bigWig file.")
parser.add_argument("output", help="Output bigWig file.")
parser.add_argument("chromosome", help="The name of the chromosome that should be placed in the output file")
args = parser.parse_args()

if not os.path.exists(args.input):
    sys.exit("{} does not exist!\n".format(args.input))

bw = pyBigWig.open(args.input)

# sanity checking, ensure the chromosome exists
if args.chromosome not in bw.chroms():
    bw.close()
    sys.exit("{} isn't a chromosome in {}!\n".format(args.chromosome, args.input))

bwO = pyBigWig.open(args.output, "w")
if args.headerToo:
    hdr = [(args.chromosome, bw.chroms(args.chromosome))]
else:
    hdr = [(x, y) for x, y in bw.chroms().items()]
bwO.addHeader(hdr)
starts = []
ends = []
vals = []
for interval in bw.intervals(args.chromosome):
    starts.append(interval[0])
    ends.append(interval[1])
    vals.append(interval[2])
bwO.addEntries([args.chromosome] * len(starts), starts, ends, vals)
bwO.close()
bw.close()


