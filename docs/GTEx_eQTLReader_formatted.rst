============
GTEx eQTL Reader
============

Function that reads the files of fine-mapped eQTLs from GTEx to ease working with the data with pybedtools.

.. autofunction:: GTEx_eQTLReader.get_eqtls

.. code-block:: python

    # Get the fine-mapped eQTLs from files subsetted for the first 1k lines from chr21.
    import GTEx_eQTLReader
    annotation = 'ExampleData/gencode.v38.annotation_chr21Genes.gtf'
    gtex_folder = 'ExampleData/GTEx_phs000424.v8.p2_chr21/'
    eqtl_beds, unique_eqtl_beds, tissue_genes = GTEx_eQTLReader.get_eqtls(hg38_annotation=annotation, gtex_folder=gtex_folder, gtex_tissues=None, max_distance=None)
    print(eqtl_beds['Adipose_Subcutaneous']['CaVEMaN'][0])
    print(unique_eqtl_beds['Adipose_Subcutaneous']['CaVEMaN'][0])
    print(list(tissue_genes['Adipose_Subcutaneous']['CaVEMaN'])[:2])
    

.. include:: gallery/src.GTEx_eQTLReader.get_eqtls_1.txt
    :literal:

.. include:: gallery/src.GTEx_eQTLReader.get_eqtls_2.txt
    :literal:

.. include:: gallery/src.GTEx_eQTLReader.get_eqtls_3.txt
    :literal:


