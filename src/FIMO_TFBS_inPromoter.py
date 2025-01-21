
import argparse

"""For each gene call TFBS in a promoter window at the 5'TSS with Fimo. Then convert Fimo's output into a matrix
of genes x TFs with the counted binding sites. Only genes with at least one binding site are written to output."""

parser = argparse.ArgumentParser()
parser.add_argument("--gtf", required=True, help='GENCODE gtf-file to take the genes and locations from.')
parser.add_argument("--extend", default=200, type=int, help="Number of bp to extend the promoter in each direction.")
parser.add_argument("--open_regions", default=False, type=str, help="Optional bed-file to restrict the promoter regions to. "
                                                                    "That means first the promoter regions are constructed and "
                                                                    "then only those promoter parts that overlap with open_regions are kept.")
parser.add_argument("--PWMs", required=True, help='TF motif file in meme format.')
parser.add_argument("--fasta", required=True, help='Genome sequence file in fasta format.')
parser.add_argument("--fimo_src", required=True, help='Path to the Fimo executable. If fimo is on path, just type "fimo".')
parser.add_argument("--out_dir", required=True, help='Path to which to write the output to.')
parser.add_argument("--write_sequence", default='False', help='If Fimo should write the sequence matched '
                                                                       'to the motif in its output file [True, False].')

if __name__ == '__main__':
    # This is ugly and a style violation, but otherwise can't convince the Sphinx markdown to document this function.
    from pybedtools import BedTool
    import GTF_Processing
    import FIMO_TFBS_Helper

    args, seq_out, fimo_out = FIMO_TFBS_Helper.process_args(parser.parse_args())

    # Get Promoter and run FIMO.
    promoter_bed = GTF_Processing.gene_window_bed(args.gtf, args.extend, gene_set=set(), dict_only=False, open_regions=args.open_regions)
    # Force the regions to be within the chromosome boundaries.
    promoter_bed = promoter_bed.slop(g=args.fasta+'.fai', b=0)
    # Regions outside will get one bp at the border.
    promoter_bed = BedTool(''.join([str(x) for x in promoter_bed if x.length > 1]), from_string=True)
    if not promoter_bed:
        print("No remaining region after capping at locations covered by the fasta file.")
        exit()
    promoter_bed.sequence(fi=args.fasta, fo=seq_out, name=True)
    print('sequence fasta stored at', seq_out)

    new_meme_file = FIMO_TFBS_Helper.meme_fitbackground(args.PWMs, seq_out, args.out_dir)

    FIMO_TFBS_Helper.fimo_runner(args, seq_out, fimo_out, new_meme_file)

    FIMO_TFBS_Helper.fimo_processor(args, fimo_out)

