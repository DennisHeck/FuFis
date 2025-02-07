============
BigWig Counter
============

Function that uses pyBigWig to get the average bigwig signal in regions from a bed-file from one or multiple bigwig files.

.. autofunction:: BigWig_Counter.bigwig_counts

.. code-block:: python

    import BigWig_Counter
    # Take a mini bed-file and get the signal from two chr21 bigwig files from IHEC (https://ihec-epigenomes.org/epiatlas/data/).
    bed_file = "ExampleData/H3K27acPeaks_chr21.narrowPeak"
    bigwigs = ['ExampleData/IHECRE00000013_chr21.bigwig', 'ExampleData/IHECRE00000017_chr21.bigwig']
    bed_counts, errors = BigWig_Counter.bigwig_counts(bed_file, bigwigs, n_cores=1)
    print(bed_counts.head())
    

.. include:: gallery/src.Bigwig_Counter.txt
    :literal:


