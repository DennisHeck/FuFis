============
FuFis (Frequently used Functions)
============

A collection of functions for bioinformatics analyses and plotting functions. All documented functions come with example
data and example code.

Installation and Introduction
*********************

To get the required packages, clone the repository and create a conda environment with the
yml-file src/condaenv.yml (and cross fingers that conda is in your favour this day). If you have `mamba <https://github.com/mamba-org/mamba>`_
installed (conda but fast), just replace the 'conda' in the commands with 'mamba'.

.. code-block:: Python

    git clone https://github.com/DennisHeck/FuFis.git
    cd FuFis/src  # Or skip the cd and give the full path to the yml-file.
    conda env create -f condaenv.yml

There might be other software that needs to be installed depending on the function, for example deeptools, which can't be easily added to the conda environment file.

Disclaimers:
 - No warranties. There are no tests to guarantee the correctness of functions. Be sceptical.
 - Double-check whether the packages and APIs of the functions originate from work that needs to be cited. Usually, this is mentioned in the docstring, but this is not kept up-to-date and does not claim to be complete.

Add your own functions:
    If you also have functions which could be useful to others, it would be great if you open a pull request to add
    it to the collection! If you submit anything, please:
     - Provide code how to run it.
     - Have the code properly annotated, ideally with a docstring that is readable by Sphinx.
     - Add example data if necessary and how the output should look like.
     - List the package requirements.

There's also an explanation of how this documentation is set up `here <file:///Users/dennis/GitHub/FuFis_Sphinx/html/DocuDocu.html.>`_.


***************************
Frequently used Flags
***************************
Some flags, especially those of plotting scripts, are used in many functions. To not bloat up the documentation
(and avoid lots of copy-pasting), those rather universal flags are listed here. Other flags that are self-explanatory,
such as *title* or *xlimit*, won't be explained.

 - **plot_path**: The path to which any plots should be saved. No format-suffix needed, that will be handled by the functions. Also, other information will be added to the final path, e.g., the label of the axes, which makes it easier to produce multiple plots without having to adjust the path to prevent them being overwritten.
 - **formats**: A list of formats in which a plot should be saved, defaults to ['pdf']. All common formats should be supported. So if you need png and svg, set it to ['png', 'svg'].
 - **palette**: This takes the usual palettes such as 'tabs20', but also has the option to name one of the categorical palettes defined by Glasbey (10.1002/col.20327), and accessed by the `colorcet package <https://colorcet.holoviz.org/user_guide/Categorical.html>`_. They are defined in the ColoursAndShapes.py script, and allow the options 'glasbey', 'glasbey_warm', 'glasbey_cool', 'glasbey_dark', 'glasbey_light'.
 - **x_size and y_size**: Dimensions of the plot (rather intuitive).
 - **font_s**: Defines a base-font size for the plot. The font size of the elements will be changed in reference to it, e.g., when the x-axis ticks are font_s in size, the title will be set to font_s+4.
 - **grid**: Whether to display a grid underneath the plot elements.
 - **vmin and vmax**: Minimum and maximum values for plots with a non-categorical colourmap.







