import os

"""Call the GalleryGenerator script to update all plots and text output, and then edit the doc-files to
insert the code blocks from GalleryGenerator as well."""

gallery_script = 'GalleryGenerator.py'
gallery_blocks = open(gallery_script).read().split("# ***")[1:]
gallery_code = {x.split('\n')[0]: '\n'.join(x.split('\n')[1:]).split('# ---')[0] for x in gallery_blocks}

rst_files = [x for x in os.listdir('.') if x.endswith('.rst') and 'formatted' not in x and x != 'index.rst'
             and x != 'Main.rst' and x != 'DocuDocu.rst']
for rst_file in rst_files:
    rst_text = open(rst_file).read()
    for gall in gallery_code:
        rst_text = rst_text.replace('*'+gall+'*', '\n    '.join([x for x in gallery_code[gall].split('\n')]))

    open(rst_file.replace('.rst', '_formatted.rst'), 'w').write(rst_text)

