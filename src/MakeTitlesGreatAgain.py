import sys
import re

"""This script outputs all title lines of a bibtex file with two {{ }} brackets around it to fix capitalization of 
most latex bibliography styles"""

if len(sys.argv) != 2:
    print("Usage python3 FixBibTexTitles.py file.bib > newFile.bib")
else:
    file = open(sys.argv[1], 'r')
    for line in file:
        line = line.strip()
        # If this is a title line with {} brackets around the title and make sure we didn't already put two brackets.
        if re.search(r"\s*title\s*=", line) and len(re.findall(r'\{.*?\}', line)) == 1 and line.count("{{") == 0:
            # Split on the { and } brackets and output with double {{ }}
            elems = re.split(r"[{}]", line)
            print("".join([elems[0], "{{", elems[1].replace('{', '').replace('}', ''), "}}", elems[2]]))
        else:
            print(line)

