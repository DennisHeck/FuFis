============
MEME formatting
============

Functions for parsing TF motif files in meme-format (https://meme-suite.org/meme/doc/meme-format.html). Note, even if
it's a defined format, there might be variations when using files that are not used in the examples here.

.. code-block:: python

    # Block that has to be executed for all.
    import MEME_Formatting
    out_dir = 'docs/gallery/'
    meme_file = 'ExampleData/Jaspar_Hocomoco_Kellis_human_meme.txt'
    annotation = 'ExampleData/gencode.v38.annotation_chr21Genes.gtf'  # It's not the full one, so fewer hits expected.
    

.. .--------------------------------------------------------------------------------------------------------------------
.. meme_id_map
.. .--------------------------------------------------------------------------------------------------------------------
.. autofunction:: MEME_Formatting.meme_id_map

.. code-block:: python

    # Get the Ensembl ID for the TFs in our motif meme-file.
    tf_ids, all_tf_names, misses = MEME_Formatting.meme_id_map(meme_file=meme_file, gtf_file=annotation, species='human')
    print('TBXT', tf_ids['TBXT'])
    print('MAX::MYC', tf_ids['MAX::MYC'])
    

.. include:: gallery/src.MEME_Formatting.meme_id_map.txt
    :literal:


.. .--------------------------------------------------------------------------------------------------------------------
.. meme_monomer_map
.. .--------------------------------------------------------------------------------------------------------------------
.. autofunction:: MEME_Formatting.meme_monomer_map

.. code-block:: python

    # Useful when in need of the individual monomers or removal of the motif versions.
    tf_monomer_map, all_monomer_names = MEME_Formatting.meme_monomer_map(meme_file=meme_file)
    print('BHLHA15(MA0607.2)', tf_monomer_map['BHLHA15(MA0607.2)'])
    print('MAX::MYC', tf_monomer_map['MAX::MYC'])
    print('all monomers', all_monomer_names[:4])
    

.. include:: gallery/src.MEME_Formatting.meme_monomer_map.txt
    :literal:


.. .--------------------------------------------------------------------------------------------------------------------
.. meme_id_map
.. .--------------------------------------------------------------------------------------------------------------------
.. autofunction:: MEME_Formatting.subset_meme

.. code-block:: python

    # Subset a meme-file, which is useful for example for excluding TFs that are not expressed.
    MEME_Formatting.subset_meme(meme_file, motif_names=['MAX::MYC', 'TBXT'], out_file=out_dir+"Subset_meme.txt",
                                include_dimers=True, exact_match=False)
    print(open(out_dir+"Subset_meme.txt").read())
    

.. include:: gallery/Subset_meme.txt
    :literal:

