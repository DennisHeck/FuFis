import subprocess
import os
from timeit import default_timer as clock

"""Function to merge bigwig files with wiggletools write that creates a wig-file that is subsequently
converted to a bigwig again. Written as callable function and not commandline-callable script to
have the list of bigwig files more flexible."""

# CARE the UCSC merger behaves weirdly, which is why it wiggletools is used here


def merge_bws(bw_files, wiggle_exe, wigToBigWig_exe, chromsize_file, bw_out):
    """
    Merge the list of bw_files. Yes, it's only a tiny wrapper, but still a wrapper.
    
    Args:
        bw_files: List of bw-files.
        wiggle_exe: Full path to the wiggletools executable. If it's on path just write 'wiggletools'.
        wigToBigWig_exe: Path to UCSC's wigToBigWig executable (e.g. '/home/dhecker/UCSC_tools/wigToBigWig').
        chromsize_file: Full path to the file with the chromosome sizes (e.g. '/projects/abcpp/work/base_data/hg38.chrom.sizes').
        bw_out: Full path where the merged bigwig should be stored (e.g. 'Folder/Merged.bw').
    """
    start = clock()
    wig_out = '.'.join(bw_out.split('.')[:-1])+'.wig'
    cleaned_wig_out = wig_out.replace(".wig", '_cleaned.wig')

    if not os.path.isfile(wig_out):
        subprocess.call("time {} write {} mean {}".format(wiggle_exe, wig_out, ' '.join(bw_files)), shell=True)
    else:
        print("WARNING: wig file already exists, trying conversion to bigwig:", wig_out)
    print("Converting the merged wig to bigwig")
    
    # Check if the wig file has the chr-prefix, if not, add it, and remove odd scaffolds.
    allowed_chr = ['X', 'Y']  # We also allow ints after chr-removal.
    with open(wig_out) as wig_in, open(cleaned_wig_out, 'w') as chr_wig:
        for line in wig_in:
            # This is a rather ugly fix, but there can be many weird lines in wig files.
            line_chr = line.split('\t')[0].replace('chr', '')
            if line_chr in allowed_chr or line_chr.isdigit():
                chr_wig.write('chr'+line.replace('chr', ''))
    os.remove(wig_out)

    subprocess.call(f"time {wigToBigWig_exe} {cleaned_wig_out} {chromsize_file} {bw_out}", shell=True)
    if os.path.isfile(bw_out):
        os.remove(cleaned_wig_out)
    print(clock() - start, "bigwigs merged")

