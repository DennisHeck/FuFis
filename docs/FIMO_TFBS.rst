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
The first script starts with a bed-file, the second with a set of genes to get the TFBS in their promoter regions.
Both are called via the command line.
If the genome sequence file doesn't have an index, `samtools <https://www.htslib.org/>`_ has to be installed and on the
PATH for the function to create it.

***************************
TFBS in regions
***************************
Start with a bed-file and get the matrix of regions x TFs with the TFBS.

.. argparse::
   :ref: FIMO_TFBS_inRegions.parser
   :prog: FIMO_TFBS_inRegions.py

.. code-block:: python

    *FIMO_TFBS_inRegions.cmd*

.. include:: gallery/src.FIMO_TFBS_inRegions.cmd.txt
    :literal:

.. code-block:: python

    *FIMO_TFBS_inRegions.matrix*

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

    *FIMO_TFBS_inPromoter*

.. include:: gallery/src.FIMO_TFBS_inPromoter.txt
    :literal:








