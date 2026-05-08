#########################
Hi-C Processing
#########################

Hi-C files are fun and very variable in their format. Here is wrapper that uses hicstraw to bring a .hic-file into coordinate format,
with one file per chromosome and each file with bin | bin | contact.
Due to the size of Hi-C files, no example is provided here. There is an alternative bash wrapper in the `STARE repository <https://github.com/SchulzLab/STARE/blob/main/Code/Juicebox_normalization.sh>`_ that uses juicer's dump.

***************************
HiCStraw
***************************

With hicstraw to write the coo-format from a .hic-file with one file per chromosome. Works with a local hic-file or remote URL.

.. argparse::
   :ref: HiCStraw.parser
   :prog: HiCStraw.py

