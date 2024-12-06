============
GTF processing
============


.. code-block:: python

    import src.GTF_Processing as GTF_Processing
    annotation = 'ExampleData/gencode.v38.annotation_chr21Genes.gtf'
    gene_list = ['ENSG00000160294', 'ENSG00000279493', 'ENSG00000279720']
    promoter_regions = GTF_Processing.gene_window_bed(gtf_file=annotation, extend=200, gene_set=gene_list, tss_type='5')
    print(promoter_regions)
    

.. include:: docs/gallery/src.GTF_Processing.txt
   :literal:

.. automodule:: GTF_Processing
   :members:
   :undoc-members:
   :show-inheritance:





