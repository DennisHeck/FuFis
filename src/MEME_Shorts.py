import subprocess
from pybedtools import BedTool

"""Collection of slightly easier access to the MEME suite tools."""


def sea(foreground, background, sequence, meme_file, file_tag, out_dir, sea_src='sea'):
    """Run MEME's SEA motif enrichment for known motifs that are relatively enriched in the foreground
    compred to the background sequences. If background is 'random', then use the shuffled foreground as background."""
    if type(foreground) == str:
        foreground = BedTool(foreground)
    if type(background) == str and background != 'random':
        background = BedTool(background)
    foreground_path = out_dir + file_tag + "_FG.fa"
    background_path = out_dir + file_tag + "_BG.fa"
    foreground.sequence(fi=sequence, fo=foreground_path, name=False)
    if type(background) == str and background.lower() == 'random':
        subprocess.call(sea_src+" --p " + foreground_path + " --m " + meme_file +
                        " --o " + out_dir + file_tag, shell=True)
    else:
        background.sequence(fi=sequence, fo=background_path, name=False)
        subprocess.call(sea_src+" --p " + foreground_path + " --n " + background_path + " --m " + meme_file +
                        " --o " + out_dir + file_tag, shell=True)


def ame(bed_obj, sequence, meme_file, file_tag, out_dir, ame_src='ame'):
    """Run MEME's AME motif enrichment for different sequence comparisons against their shuffled background."""
    if type(bed_obj) == str:
        bed_obj = BedTool(bed_obj)
    bed_obj.sequence(fi=sequence, fo=out_dir + file_tag + ".fa", name=False)
    subprocess.call(ame_src+"  --control --shuffle-- --o " + out_dir + file_tag + " " + out_dir + file_tag + ".fa "
                    + meme_file, shell=True)


def streme(foreground, background='', sequence='', file_tag='', out_dir='', nmotifs=None, streme_src='streme'):
    """Run MEME's STREME to discover motifs that are enriched in the foreground compared to background."""
    if type(foreground) == str:
        foreground = BedTool(foreground)
    if background and type(background) == str:
        background = BedTool(background)
        background_path = out_dir + file_tag + "_BG.fa"
        background.sequence(fi=sequence, fo=background_path, name=False)
        background_flag = ' --n ' + background_path
    else:
        background_flag = ''
    if nmotifs:
        nmotifs_flag = ' --nmotifs ' + str(nmotifs)
    else:
        nmotifs_flag = ''
    foreground_path = out_dir + file_tag + "_FG.fa"
    foreground.sequence(fi=sequence, fo=foreground_path, name=False)
    streme_cmd = streme_src+" --p " + foreground_path + background_flag + " --o " + out_dir + file_tag + " --minw 6" + nmotifs_flag
    print(streme_cmd)
    subprocess.call(streme_cmd, shell=True)

