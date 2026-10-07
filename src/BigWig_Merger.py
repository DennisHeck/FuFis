import subprocess
import os
from timeit import default_timer as clock

"""Function to merge bigwig files with wiggletools write_bg that creates a bedGraph-file that is subsequently
converted to a bigwig again. Written as callable function and not commandline-callable script to
have the list of bigwig files more flexible."""

# CARE using wiggletools write can create 'holes' in the track, write_bg circumvents this for whatever reason.
# CARE the UCSC merger behaves weirdly, which is why wiggletools is used here.


def merge_bws(bw_files, wiggle_exe, bedGraphToBigWig_exe, chromsize_file, bw_out):
    """
    Merge the list of bw_files. Yes, it's only a tiny wrapper, but still a wrapper.
    
    Args:
        bw_files: List of bw-files.
        wiggle_exe: Full path to the wiggletools executable (e.g. /home/dhecker/miniconda3/envs/wiggletools/bin/wiggletools). If it's on path just write 'wiggletools'.
        bedGraphToBigWig_exe: Path to UCSC's bedGraphToBigWig executable (e.g. '/home/dhecker/UCSC_tools/bedGraphToBigWig').
        chromsize_file: Full path to the file with the chromosome sizes (e.g. '/projects/abcpp/work/base_data/hg38.chrom.sizes').
        bw_out: Full path where the merged bigwig should be stored (e.g. 'Folder/Merged.bw').
    """
    start = clock()
    bg_out = '.'.join(bw_out.split('.')[:-1])+'.bg'
    cleaned_bg_out = bg_out.replace(".bg", '_cleaned.bg')

    if not os.path.isfile(bg_out):
        subprocess.call("time {} write_bg {} mean default 0 {}".format(wiggle_exe, bg_out, ' default 0 '.join(bw_files)), shell=True)
    else:
        print("WARNING: bg file already exists, trying conversion to bigwig:", bg_out)
    print("Converting the merged bg to bigwig")

    # Check if the bg file has the chr-prefix, if not, add it, and remove odd scaffolds.
    allowed_chr = ['X', 'Y']  # We also allow ints after chr-removal.
    with open(bg_out) as wig_in, open(cleaned_bg_out, 'w') as chr_wig:
        for line in wig_in:
            # This is a rather ugly fix, but there can be many weird lines in wig files.
            line_chr = line.split('\t')[0].replace('chr', '')
            if line_chr in allowed_chr or line_chr.isdigit():
                chr_wig.write('chr'+line.replace('chr', ''))
    os.remove(bg_out)

    subprocess.call(f"time {bedGraphToBigWig_exe} {cleaned_bg_out} {chromsize_file} {bw_out}", shell=True)
    if os.path.isfile(bw_out):
        os.remove(cleaned_bg_out)
    print(clock() - start, "bigwigs merged")

