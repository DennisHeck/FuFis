============
UniProt API
============

Wrapper function for the UniProt API.


.. autofunction:: UniProt_API.uniprot_domains

.. code-block:: python

    import UniProt_API
    # Look up the annotated domains and regions of three examples.
    protein_domains, protein_regions, missed_proteins, failed_requests = UniProt_API.uniprot_domains(protein_names=['KDM6A', 'DNMT3A', 'STAT2'], species='human', n_cores=1)
    print(protein_domains)
    print(protein_regions)
    


.. include:: gallery/src.UniProtAPI.uniprot_domains.txt
    :literal:
