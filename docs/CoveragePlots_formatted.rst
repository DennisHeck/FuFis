============
Coverage Plots
============

Ease the use of the neat coverage plots from deeptools.

.. autofunction:: CoveragePlots.plotHeatmap

.. code-block:: python

    from pybedtools import BedTool
    import CoveragePlots
    out_dir = 'docs/gallery/'
    # Let's plot the signal of two bigwig files from IHEC (https://ihec-epigenomes.org/epiatlas/data/) in a small set of peaks and compare that to the signal in their shuffled locations.
    peaks = BedTool("ExampleData/H3K27acPeaks_chr21.narrowPeak")
    shuffled_peaks = peaks.shuffle(genome='hg38', chrom=True, seed=12)
    bigwigs = ['ExampleData/IHECRE00000013_chr21.bigwig', 'ExampleData/IHECRE00000017_chr21.bigwig']
    CoveragePlots.plotHeatmap(beds_to_plot=[peaks, shuffled_peaks], bed_labels=['Original', 'Shuffled'], bigwigs=bigwigs, bw_labels=['Sample1', 'Sample2'],
                              out_dir=out_dir, out_tag='ExampleCoveragePlot', mode='scale', perGroup=True, title='',
                              scaled_size=500, start_label='Peak start', end_label='Peak end')
    

.. image:: gallery/ExampleCoveragePlot_scale_perGroup_Heatmap.png
   :width: 60%
