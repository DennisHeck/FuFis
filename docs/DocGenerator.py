
"""Call the GalleryGenerator script to update all plots and text output, and then edit the doc-files to
insert the code blocks from GalleryGenerator as well."""

rst_file = 'docs/Main.rst'
gallery_script = 'docs/GalleryGenerator.py'

gallery_blocks = open(gallery_script).read().split("# ***")[1:]
gallery_code = {x.split('\n')[0]: '\n'.join(x.split('\n')[1:]).split('# ---')[0] for x in gallery_blocks}

rst_text = open(rst_file).read()
for gall in gallery_code:
    rst_text = rst_text.replace('*'+gall+'*', gallery_code[gall])

open(rst_file.replace('.rst', '_replaced.rst'), 'w').write(rst_text)

