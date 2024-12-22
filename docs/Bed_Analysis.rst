============
Analysis of bed-files
============

Analyses related to bed-files, where the peaks are located with respect to genes, or how multiple bed-files overlap.
Centred on the usage of `pybedtools <https://daler.github.io/pybedtools/>`_, but also works with stored bed-files.
To run the example code for this module, we will always start with this block of code:

.. code-block:: python

    *Bed_Analysis.base_code*


.. autofunction:: Bed_Analysis.gene_location_bpwise

.. code-block:: python

    *Bed_Analysis.gene_location_bpwise*

|pic1| |pic2|

.. |pic1| image:: gallery/Example_peaks_GeneFeatureLocation_bpwiseOverlap_PieChart.png
   :width: 45%

.. |pic2| image:: gallery/InclExternalExample_peaks_GeneFeatureLocation_bpwiseOverlap_PieChart.png
   :width: 45%


.. autofunction:: Bed_Analysis.intersection_heatmap

.. image:: gallery/IGV_MultiBed.png
   :width: 90%

.. code-block:: python

    *Bed_Analysis.intersection_heatmap*

.. image:: gallery/_MultiIntersectHeat.png
   :width: 90%


.. autofunction:: Bed_Analysis.upset_to_reference

.. code-block:: python

    *Bed_Analysis.upset_to_reference*

.. image:: gallery/_UpSet_Subset_peaks.png
   :width: 80%




