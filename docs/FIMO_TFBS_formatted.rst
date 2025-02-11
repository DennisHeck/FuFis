#########################
FIMO TFBSs
#########################

`FIMO <https://meme-suite.org/meme/doc/fimo.html>`_ is a popular tool for predicting TFBS in regions. Here are scripts that make calling FIMO easier, as they
take over some preprocessing step and also process the output into a matrix of regions x TFs.
The preprocessing includes
  - writing the sequence of the regions
  - forcing the regions to be within the genome boundaries
  - adjusting the TF motif meme-file so that the base content of the bed-regions is used as background frequencies.

In the output matrix, overlapping TFBS of the same TF on the same strand are merged and counted as 1.
Take care with palindromic motifs, they are still counted on both strands.Only the regions with at least
one binding site are written into the matrix.
The first script starts with a bed-file, the second with a set of genes to get the TFBS in their promoter regions.
Both are called via the command line.
If the genome sequence file doesn't have an index, `samtools <https://www.htslib.org/>`_ has to be installed and on the
PATH for the function to create it.

CARE: Tested for meme5.4.1, later versions handle sequence names differently.

***************************
TFBS in regions
***************************
Start with a bed-file and get the matrix of regions x TFs with the TFBS.

.. argparse::
   :ref: FIMO_TFBS_inRegions.parser
   :prog: FIMO_TFBS_inRegions.py

.. code-block:: python

    import subprocess
    import pandas as pd
    
    # We use subprocess here to run the bash command for easier documentation. It's also possible to run it directly
    # in the terminal.
    bed_file = 'ExampleData/chr21_MockPeaks.bed'
    PWMs = 'ExampleData/Jaspar_Hocomoco_Kellis_human_meme.txt'
    fasta = 'ExampleData/chr21_MiniMock.fa'  # For the example it's a random sequence from chr21.
    fimo_src = "fimo"  # If it's on the PATH, otherwise full path to the executable.
    out_dir = 'docs/gallery/FIMO_inRegions'
    
    fimo_cmd = 'python3 src/FIMO_TFBS_inRegions.py --bed_file {} --PWMs {} --fasta {} --fimo_src {} --out_dir {} --write_sequence False'.format(bed_file, PWMs, fasta, fimo_src, out_dir)
    print(fimo_cmd)
    subprocess.call(fimo_cmd, shell=True)
    

.. include:: gallery/src.FIMO_TFBS_inRegions.cmd.txt
    :literal:

.. code-block:: python

    # The output matrix will be stored in the following path.
    fimo_matrix_file = 'docs/gallery/FIMO_inRegions/Out_Fimo_TFBSMatrix.txt.gz'
    fimo_matrix = pd.read_table(fimo_matrix_file, sep='\t', header=0, index_col=0)
    # With the tiny example most TF have no binding site, so let's pick a few that have some.
    print(fimo_matrix[['FOXH1', 'STAT1', 'STAT3']])
    

.. include:: gallery/src.FIMO_TFBS_inRegions.matrix.txt
    :literal:


***************************
TFBS in promoter
***************************
Start with a gtf-file and get the TFBS in their promoter regions as gene x TF matrix.

.. argparse::
   :ref: FIMO_TFBS_inPromoter.parser
   :prog: FIMO_TFBS_inPromoter.py

.. code-block:: python

    import subprocess
    import pandas as pd
    
    # The run is similar to the previous, but we use a gtf-file instead of a bed-file.
    gtf_file = 'ExampleData/gencode.v38.annotation_chr21Genes.gtf'  # With tow mock genes covered by the mock fasta.
    PWMs = 'ExampleData/Jaspar_Hocomoco_Kellis_human_meme.txt'
    fasta = 'ExampleData/chr21_MiniMock.fa'  # For the example it's a random sequence from chr21.
    fimo_src = "fimo"  # If it's on the PATH, otherwise full path to the executable.
    out_dir = 'docs/gallery/FIMO_inPromoter'
    
    fimo_cmd = 'python3 src/FIMO_TFBS_inPromoter.py --gtf {} --PWMs {} --fasta {} --fimo_src {} --out_dir {} --write_sequence False'.format(gtf_file, PWMs, fasta, fimo_src, out_dir)
    subprocess.call(fimo_cmd, shell=True)
    
    fimo_matrix_file = 'docs/gallery/FIMO_inPromoter/Out_Fimo_TFBSMatrix.txt.gz'
    fimo_matrix = pd.read_table(fimo_matrix_file, sep='\t', header=0, index_col=0)
    # With the tiny example most TF have no binding site, so let's pick a few that have some.
    print(fimo_matrix[['FOXH1', 'REST', 'PLAG1']])
    

.. include:: gallery/src.FIMO_TFBS_inPromoter.txt
    :literal:








