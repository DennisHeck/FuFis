import subprocess
from pybedtools import BedTool

"""Collection of slightly easier access to the MEME suite tools."""


def sea(foreground, background, sequence, meme_file, file_tag, out_dir):
    """Run MEME's SEA motif enrichment for known motifs that are relatively enriched in the foreground
    compred to the background sequences. If background is 'random', then use the shuffled foreground as background."""
    if type(foreground) == str:
        foreground = BedTool(foreground)
    if type(background) == str and background != 'random':
        background = BedTool(background)
    foreground_path = out_dir + file_tag + "_FG.fa"
    background_path = out_dir + file_tag + "_BG.fa"
    foreground.sequence(fi=sequence, fo=foreground_path, name=True)
    if type(background) == str and background.lower() == 'random':
        subprocess.call("sea --p " + foreground_path + " --m " + meme_file +
                        " --o " + out_dir + file_tag, shell=True)
    else:
        background.sequence(fi=sequence, fo=background_path, name=True)
        subprocess.call("sea --p " + foreground_path + " --n " + background_path + " --m " + meme_file +
                        " --o " + out_dir + file_tag, shell=True)


def ame(bed_obj, sequence, meme_file, file_tag, out_dir):
    """Run MEME's AME motif enrichment for different sequence comparisons against their shuffled background."""
    if type(bed_obj) == str:
        bed_obj = BedTool(bed_obj)
    bed_obj.sequence(fi=sequence, fo=out_dir + file_tag + ".fa", name=True)
    subprocess.call("ame  --control --shuffle-- --o " + out_dir + file_tag + " " + out_dir + file_tag + ".fa "
                    + meme_file, shell=True)


def streme(foreground, background, sequence, file_tag, out_dir):
    """Run MEME's STREME to discover motifs that are enriched in the foreground compared to background."""
    if type(foreground) == str:
        foreground = BedTool(foreground)
    if type(background) == str:
        background = BedTool(background)
    foreground_path = out_dir + file_tag + "_FG.fa"
    background_path = out_dir + file_tag + "_BG.fa"
    foreground.sequence(fi=sequence, fo=foreground_path, name=True)
    background.sequence(fi=sequence, fo=background_path, name=True)

    subprocess.call("streme --p " + foreground_path + " --n " + background_path + " --o " + out_dir + file_tag + " --minw 6", shell=True)