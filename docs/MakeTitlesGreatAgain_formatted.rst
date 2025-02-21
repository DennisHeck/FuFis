#########################
Make Titles Great Again
#########################

Script to fix wrong capitalisation in latex reference bib-files. It prevents that originally fully upper case
words in the title are compiled into words with only the first letter capitalized, e.g. QUOKKA to Quokka.

Originates from https://github.com/SchulzLab/MakeTitlesGreatAgain .

.. code-block:: python

    # Run it with a mini example with one bib entry, here from Python, but you can also run the command
    # that's sent to subprocess in your bash.
    import subprocess
    subprocess.call("python3 src/MakeTitlesGreatAgain.py ExampleData/ExampleBib.bib > ExampleData/ExampleBib_great.bib", shell=True)
    print(open('ExampleData/ExampleBib.bib').readlines()[2])
    print(open('ExampleData/ExampleBib_great.bib').readlines()[2])
    

.. include:: gallery/MakeTitlesGreatAgain_1.txt
    :literal:

.. include:: gallery/MakeTitlesGreatAgain_2.txt
    :literal:




