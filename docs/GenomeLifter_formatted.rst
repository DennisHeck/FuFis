============
GenomeLifter
============

A handy function to lift genomic coordinates between genome versions.

.. autofunction:: GenomeLifter.genome_lifter

.. code-block:: python

    import GenomeLifter
    from pybedtools import BedTool
    
    # Lift an example list of hg38 regions to hg19.
    hg38_bed_file = "ExampleData/H3K27acPeaks_chr21.narrowPeak"
    hg38_regions = BedTool(hg38_bed_file)
    print('hg38 coordinates:')
    print(''.join([str(x) for x in hg38_regions[:3]]))
    hg19_regions, unliftable = GenomeLifter.genome_lifter(hg38_regions, input_version='hg38', output_version='hg19')
    print('hg19 coordinates:')  # Note the output is now a list.
    print(hg19_regions[:3])
    

.. include:: gallery/src.GenomeLifter.genome_lifter_hg38.txt
    :literal:

.. include:: gallery/src.GenomeLifter.genome_lifter_hg19.txt
    :literal:









