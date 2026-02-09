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

    subprocess.call("time {} write {} mean {}".format(wiggle_exe, wig_out, ' '.join(bw_files)), shell=True)
    print("Converting the merged wig to bigwig")
    subprocess.call("time {} {} {} {}".format(wigToBigWig_exe, wig_out, chromsize_file, wig_out.replace('.wig', '.bw')), shell=True)
    if os.path.isfile(wig_out.replace('.wig', '.bw')):
        os.remove(wig_out)
    print(clock() - start, "bigwigs merged")
