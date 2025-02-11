============
Documentation for the documentation
============

General setup
*************

To understand how to add functions to the documentation or to adjust it, here a description of the general setup:

- The documentation is compiled by `Sphinx <https://www.sphinx-doc.org/en/master/index.html#>`_, a Python application to generate documentations.
- `Read the Docs <https://about.readthedocs.com/>`_ can host documentations via GitHub and uses Sphinx to build them.
- If a GitHub repository is hooked on Read the Docs, each push will trigger a new build. Its status can be checked on the Read the Docs page.
- During testing, it is more practical to use a Sphinx installation on your local machine and do the building of the html files there as it's faster.


Folder structure
*************

There are currently three main folders in the repository:

| ├── **src**: Holds all the scripts
| ├── **docs**: Holds the files required for the documentation, including the markdown .rst files.
|     └── **gallery**: All the files loaded in the html files, such as plots or text files, most of them written by docs/GalleryGenerator.py
| ├── **ExampleData**: Folder with the files required to run the examples in docs/GalleryGenerator.py

The autodoc function
*************

A big advantage of Sphinx are the extensions, specifically the autodoc functions. That means that when the docstring of a function in Python is
written in a format that Sphinx understands, it only requires one line to fully display the docstring in a neat format. This also means, that whenever a commit
changes the docstring, the documentation will also be updated. For example, we have
a function in our Python script written such as:

.. code-block:: Python

    def myfunction1(arg1, arg2):
        """
        This is a toy example function.

        Args:
            arg1: First argument, very important.
            arg2: Second argument, very required.
        """

In the markdown files (ending with .rst), we'll only need a one-liner to show the docu:

.. code-block:: Python

    .. autofunction:: DocuExamples.myfunction1

The line above will be rendered in the html-file as follows:

.. autofunction:: DocuExamples.myfunction1

This also works for scripts that use argparse to be called via the command line, although the syntax is different, see https://github.com/DennisHeck/FuFis/blob/main/docs/FIMO_TFBS.rst which is the markdown file for
the script https://github.com/DennisHeck/FuFis/blob/main/src/FIMO_TFBS_inRegions.py.


Specifics of this documentation
*************

Usually, you would have to copy-paste code for running examples, so that you have it once in a Python script where you can actually
run it to get the output, and once in the rst files as code-block so that people can read it. Apparently, there's extensions
for running Python code from rst files, but installation hasn't worked yet, and managing the packages for the code execution is what I imagine the
second level of hell to be.
To make the documentation easier and avoid copy-pasting code, there's a specific setup at play here. First, there's the script docs/GalleryGenerator.py, out of
which chunks of code can be executed to create the files in docs/gallery for the documentation to load. The specialty is, that the code blocks in said script
are labelled with a made-up syntax, e.g.

.. code-block:: Python

    # ***Script1.function1
    # Generate examples to show.
    output = Script1.function1(args)
    # ---

Instead of copy-pasting this code block into an rst file, you only need to write the following into the rst file:

.. code-block:: Python

    .. code-block:: python

       *Script1.function1*

Because when running the html-build the script docs/DocGenerator.py will automatically run and fill in the code blocks at the
respective places. It creates a new rst file with the suffix _formatted.rst, which is the file you should be adding to
the docs/index.rst.


Steps for adding a documentation page
*************

Given the specifics and structure of this documentation, here the steps if you want to add a documentation page for a script:

1. Add code to docs/GalleryGenerator.py that showcases the functions you want to document and generates the files you'd like to be shown.
2. Create a .rst file in docs/, e.g. MyScript.rst. See other example .rst-files for the syntax, look at those without the _formatted suffix (see section above).
3. Add MyScript_formatted.rst (see section above) to the docs/index.rst.
4. Build the documentation locally or push to GitHub to trigger a build by Read the Docs.

A few words on automock
*************

When using autodock, it will try to import the whole scripts, including all the packages. Because package management is great, this will
fail. It's also not necessary, since you can autodock your functions without getting all the imports to work. That's where mock imports come into play.
Sphinx allows to list the packages that should be ignored and not tried to be installed. However,
listing them all by hand for all scripts is a nuisance (I tried). Currently, there are a few lines of code in docs/conf.py which fetch
all imports from Python scripts that have matching rst files and add them to the list of mock imports. If your documentation is not showing,
it might be because of import issues. It is quite unstable, and there are very weird instances where combinations of packages crash but the
individual ones work.




