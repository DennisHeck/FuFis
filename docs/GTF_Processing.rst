============
GTF processing
============

Functions based on gtf-files (tested for GENCODE's format) and the resulting outputs.
Centred on the usage of `pybedtools <https://daler.github.io/pybedtools/>`_, but
the output can also be easily written to bed files.
To run the example code for this module, we will always start with this block of code:

.. code-block:: python

    *GTF_Processing.base_code*

.. autofunction:: GTF_Processing.gene_window_bed

.. code-block:: python

    *GTF_Processing.gene_window_bed*

.. include:: gallery/src.GTF_Processing.gene_window_bed.txt
    :literal:


.. autofunction:: GTF_Processing.gene_feature_table

.. code-block:: python

    *GTF_Processing.gene_feature_table*

.. include:: gallery/src.GTF_Processing.gene_feature_table.txt
    :literal:


.. autofunction:: GTF_Processing.gene_body_bed

.. code-block:: python

    *GTF_Processing.gene_body_bed*

.. include:: gallery/src.GTF_Processing.gene_body_bed.txt
    :literal:


.. autofunction:: GTF_Processing.gene_feature_bed

.. code-block:: python

    *GTF_Processing.gene_feature_bed*

.. include:: gallery/src.GTF_Processing.gene_feature_bed.txt
    :literal:


.. autofunction:: GTF_Processing.gene_introns_bed

.. code-block:: python

    *GTF_Processing.gene_introns_bed*

.. include:: gallery/src.GTF_Processing.gene_introns_bed.txt
    :literal:


.. autofunction:: GTF_Processing.match_gene_identifiers

.. code-block:: python

    *GTF_Processing.match_gene_identifiers.ensembl*

.. include:: gallery/src.GTF_Processing.match_gene_identifiers.ensembl.txt
    :literal:

.. code-block:: python

    *GTF_Processing.match_gene_identifiers.symbols*

.. include:: gallery/src.GTF_Processing.match_gene_identifiers.symbols.txt
    :literal:

.. code-block:: python

    *GTF_Processing.match_gene_identifiers.entrez*

.. include:: gallery/src.GTF_Processing.match_gene_identifiers.entrez.txt
    :literal:


.. autofunction:: GTF_Processing.gene_biotypes

.. code-block:: python

    *GTF_Processing.gene_biotypes*

.. include:: gallery/src.GTF_Processing.gene_biotypes.txt
    :literal:


.. autofunction:: GTF_Processing.bed_to_length_dict

.. code-block:: python

    *GTF_Processing.bed_to_length_dict*

.. include:: gallery/src.GTF_Processing.bed_to_length_dict.txt
    :literal:


.. autofunction:: GTF_Processing.bed_to_feature_dict

.. code-block:: python

    *GTF_Processing.bed_to_feature_dict*

.. include:: gallery/src.GTF_Processing.bed_to_feature_dict.txt
    :literal:









