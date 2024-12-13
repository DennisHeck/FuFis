============
GTF processing
============

Functions based on gtf-files (tested for GENCODE's format) and the resulting outputs.
Centred on the usage of `pybedtools <https://daler.github.io/pybedtools/>`_, but
the output can also be easily written to bed files.
To run the example code for this module, we will always start with this block of code:

.. code-block:: python

    import src.GTF_Processing as GTF_Processing
    annotation = 'ExampleData/gencode.v38.annotation_Mini.gtf'
    gene_list = ['ENSG00000160294', 'ENSG00000279493', 'ENSG00000279720']
    

.. autofunction:: GTF_Processing.gene_window_bed

.. code-block:: python

    promoter_regions = GTF_Processing.gene_window_bed(gtf_file=annotation, extend=200, gene_set=gene_list, tss_type='5')
    print(promoter_regions)
    

.. include:: gallery/src.GTF_Processing.gene_window_bed.txt
    :literal:


.. autofunction:: GTF_Processing.gene_feature_table

.. code-block:: python

    # Gets the biotypes as annotated in a gtf-file, and also provides a high-level annotation.
    feature_table = GTF_Processing.gene_feature_table(gtf_file=annotation, gene_set=gene_list)
    print(feature_table)
    

.. include:: gallery/src.GTF_Processing.gene_feature_table.txt
    :literal:


.. autofunction:: GTF_Processing.gene_body_bed

.. code-block:: python

    gene_bodies = GTF_Processing.gene_body_bed(gtf_file=annotation, gene_set=gene_list, dict_only=False)
    print(gene_bodies)
    

.. include:: gallery/src.GTF_Processing.gene_body_bed.txt
    :literal:


.. autofunction:: GTF_Processing.gene_feature_bed

.. code-block:: python

    gene_exons = GTF_Processing.gene_feature_bed(gtf_file=annotation, feature='exon', gene_set=gene_list, dict_only=False,
                                                 merge=False, keep_strand=False)
    # This is list has 118 entries, so let's only print the first 5.
    print(''.join([str(x) for x in gene_exons[:5]]))
    

.. include:: gallery/src.GTF_Processing.gene_feature_bed.txt
    :literal:


.. autofunction:: GTF_Processing.gene_introns_bed

.. code-block:: python

    # Getting the introns is a bit more tricky, as they are not explicitly annotated. We get them by subtracting all
    # other annotations from the gene bodies.
    gene_introns = GTF_Processing.gene_introns_bed(gtf_file=annotation, gene_set=gene_list)
    print(''.join([str(x) for x in gene_introns[:5]]))
    

.. include:: gallery/src.GTF_Processing.gene_introns_bed.txt
    :literal:


.. autofunction:: GTF_Processing.match_gene_identifiers

.. code-block:: python

    # The best task of all, map different identifiers. Let's start with Ensembl IDs to gene symbols, which can often
    # be successfully done by a gtf-file alone.
    gene_ids = ['ENSG00000160294', 'ENSG00000279493', 'ENSG00000279720']
    mapped_ids, missed_ids = GTF_Processing.match_gene_identifiers(gene_ids, gtf_file=annotation, species='human',
                                                                   fields="symbol")
    print(mapped_ids)
    

.. include:: gallery/src.GTF_Processing.match_gene_identifiers.ensembl.txt
    :literal:

.. code-block:: python

    # Next we start from symbols, which is more prone to failing. We can also query for the Entrez ID.
    gene_symbols = ['AXL', 'MYO10', 'ATP5SL']
    mapped_symbols, missed_symbols = GTF_Processing.match_gene_identifiers(gene_symbols, gtf_file=annotation, species='human',
                                                                           fields="ensembl,entrezgene")
    print(mapped_symbols)
    

.. include:: gallery/src.GTF_Processing.match_gene_identifiers.symbols.txt
    :literal:

.. code-block:: python

    # If we want to lookup Entrez IDs, we have to add entrezgene to the scope in which mygene looks.
    gene_entrez = [4651, 558]
    mapped_entrez, missed_entrez = GTF_Processing.match_gene_identifiers(gene_entrez, gtf_file=annotation, species='human',
                                                                         scopes='symbol,entrezgene', fields="ensembl")
    print(mapped_entrez)
    

.. include:: gallery/src.GTF_Processing.match_gene_identifiers.entrez.txt
    :literal:


.. autofunction:: GTF_Processing.gene_biotypes

.. code-block:: python

    # Gets the biotypes as annotated in a gtf-file, and also provides a high-level annotation.
    biotype_dict = GTF_Processing.gene_biotypes(gtf_file=annotation, gene_set=gene_list)
    print(biotype_dict)
    

.. include:: gallery/src.GTF_Processing.gene_biotypes.txt
    :literal:


.. autofunction:: GTF_Processing.bed_to_length_dict

.. code-block:: python

    # Convenient, for example, to get the total exon length.
    gene_exons = GTF_Processing.gene_feature_bed(gtf_file=annotation, feature='exon', gene_set=gene_list, dict_only=False,
                                                 merge=False, keep_strand=False)
    exon_lengths = GTF_Processing.bed_to_length_dict(gene_exons)
    print(exon_lengths)
    

.. include:: gallery/src.GTF_Processing.bed_to_length_dict.txt
    :literal:


.. autofunction:: GTF_Processing.bed_to_feature_dict

.. code-block:: python

    # Admittedly, has quite specific use cases. It can be handy to collect all locations on gene level or on the level of
    # any other identifier.
    gene_exons = GTF_Processing.gene_feature_bed(gtf_file=annotation, feature='exon', gene_set=gene_list, dict_only=False,
                                                 merge=False, keep_strand=False)
    exon_dict = GTF_Processing.bed_to_feature_dict(gene_exons)
    print(exon_dict['ENSG00000279493'])
    

.. include:: gallery/src.GTF_Processing.bed_to_feature_dict.txt
    :literal:









