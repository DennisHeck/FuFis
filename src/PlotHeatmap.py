import pandas as pd
import os
import subprocess
import pybedtools
import gzip
from pybedtools import BedTool

"""Wrapper for deeptool's computeMatrix followed by plotHeatmap."""


def plotHeatmap(beds_to_plot, bed_labels, bigwigs, bw_labels, out_dir, out_tag, mode, perGroup=True, title='',
                scaled_size=1000, start_label='TSS', end_label='TES', referencePoint='center', vmin='auto', vmax='auto', matrix_min=None,
                upstream=200, downstream=200, n_cores=1, cmap='plasma', legend_loc='best', show='plot, heatmap and colorbar'):
    """
    @param beds_to_plot: Either a BedTool object or the path to a bed-file. Note that if strand is a relevant
    information it should be in the 6th column (1-based).
    @param bed_labels: List of labels for the bed files, have to have the same order.
    @param bigwigs: List of bigwig files.
    @param bw_labels: Labels for the bigwig files, have to have the same order.
    @param out_dir: Path to write the intermediary files, including bedtool's objects and plots to.
    @param out_tag: Tag to label the output plots.
    @param mode: Either 'scale' or 'reference', 'scale' rescales different regions to the same size, while
    'reference' extends from one position of the regions.
    @param perGroup: If True plot the coverage of the bigwig files underneath another in one column instead
    @param scaled_size: For 'scale' mode to what size the regions should be scaled to.
    @param start_label: For 'scale' mode how to label the start and end position.
    @param end_label: Same as above.
    of making a column for each bigwig file.
    @param referencePoint: From where to extend the regions in the 'reference' mode:
    'TSS' for region start, 'TES' for region end or 'center'.
    @param vmin: Minimum value for the heatmap. The curves are unaffected by this.
    @param vmax: Maximum value for the heatmap.
    @param matrix_min: CURRENTLY NOT WORKING. Values below this value are set to this value in the matrix file, so that it's reflected in the curves.
    @param upstream: How far to extend the regions upstream in 'reference' mode.
    @param downstream: How far to extend the regions downstream in 'reference' mode.
    @param n_core: Number of cores to compute the matrix that's then used for plotting.
    @param cmap: Colormap used for the visualising the signal.
    @param legend_loc: best, upper-right, upper-left, upper-center, lower-left, lower-right, lower-center, center,
    center-left, center-right, none.
    @param show: Which parts of the plot to generate: plot, heatmap and colorbar; plot and heatmap; heatmap only; heatmap and colorbar
    """

    pybedtools.helpers.set_tempdir(out_dir)

    start_label = start_label.replace('-', '\-')  # Otherwise it crashes the process call.
    end_label = end_label.replace('-', '\-')

    if type(beds_to_plot) != list:
        beds_to_plot = list(beds_to_plot)
    if type(bed_labels) != list:
        bed_labels = list(bed_labels)

    # A few checks of the input.
    for bw in bigwigs:
        if not os.path.isfile(bw):
            print("ERROR bw not existing:", bw)
            return

    if mode != 'scale' and mode != 'reference':
        print("ERROR: mode can only be 'scale' or 'reference'")
        return

    bed_paths = []
    filtered_labels = []
    for b, bed in enumerate(beds_to_plot):
        if type(bed) == BedTool:
            if len(bed) > 0:
                temp_bed_file = out_dir + "/temp_bed"+str(b)+".bed"
                open(temp_bed_file, 'w').write(str(bed))
                bed_paths.append(temp_bed_file)
                filtered_labels.append(bed_labels[b])
        else:
            if not os.path.isfile(bed):
                print("ERROR can't find bed-file:", bed)
                return
            bed_paths.append(beds_to_plot)
            filtered_labels.append(bed_labels[b])

    if not bed_paths:
        print("ERROR: no valid bed-files remaining")
        return

    # First, we need to call computeMatrix that gets the coverage into a matrix format that is read
    # in the next step from the plotHeatmap function.
    matrix_out = out_dir + out_tag + "_" + mode + "_perGroup"*perGroup + '.gz'
    print('Computing matrix')
    if mode == 'scale':
        matrix_cmd = "computeMatrix scale-regions -S "+' '.join(bigwigs)+" -R "+ ' '.join(["'"+b+"'" for b in bed_paths])+\
                        " -R " + ' '.join(["'"+b+"'" for b in bed_paths])+\
                        " --regionBodyLength "+str(scaled_size)+" --outFileName "+matrix_out+\
                        " --startLabel "+start_label+" --endLabel "+end_label+" --binSize 10 -p " + str(n_cores)
    elif mode == 'reference':
        matrix_cmd = "computeMatrix reference-point -S "+' '.join(bigwigs)+" -R "+ ' '.join(["'"+b+"'" for b in bed_paths])+\
                        ' --referencePoint ' + referencePoint + " --upstream "+str(upstream) +\
                        " --downstream "+str(downstream)+\
                        " --outFileName "+matrix_out+" --binSize 10 -p " + str(n_cores)
    print(matrix_cmd)
    subprocess.call(matrix_cmd, shell=True)

    if matrix_min is not None:  # NOTE Currently not working, no idea why, the matrix is not accepted afterwards.
        # Open the matrix file and set all values <forced_min to forced_min.
        header_row = gzip.open(matrix_out, 'rt').readline()
        matrix_df = pd.read_table(matrix_out, sep='\t', header=None, skiprows=1)
        val_cols = matrix_df.iloc[:, 6:]  # Bit too many lines here, but pandas was refusing to do it in one.
        val_cols[val_cols < matrix_min] = matrix_min
        minned_matrix = pd.concat([matrix_df.iloc[:, :6], val_cols], axis=1)
        print(matrix_df.shape, matrix_df.shape)
        open(matrix_out.replace('.gz', ''), 'w').write(header_row)
        matrix_df.to_csv(matrix_out.replace('.gz', ''), sep='\t', header=False, index=False, mode='a')
        subprocess.call('gzip -f '+ matrix_out.replace('.gz', ''), shell=True)


    # Call the plotting function.
    print('Plotting heatmap')
    heatmap_cmd = "plotHeatmap -m " + matrix_out + ' --perGroup'*perGroup + (" --plotTitle '" + title+"'")*bool(title) +\
                    ' --startLabel ' + start_label + ' --endLabel ' + end_label + ' --colorMap ' + cmap + " --zMin "+str(vmin)+" --zMax "+str(vmax)+\
                    " --yAxisLabel 'coverage'" + " --regionsLabel "+ ' '.join(["'"+b+"'" for b in filtered_labels])+\
                    ' --samplesLabel '+' '.join(bw_labels)+' --outFileName ' + matrix_out.replace('.gz', '_Heatmap.pdf')+ \
                    " --legendLocation " + legend_loc.lower() + " --whatToShow "+"'"+show+"'"
    print(heatmap_cmd)
    subprocess.call(heatmap_cmd, shell=True)


