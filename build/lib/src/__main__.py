################################################################################
# SCiMS: Sex-determination for metagenomic sequences 
# 
# Author: Kobie Kirven
# Davenport Lab - Penn State University
# Date: 9-2-2021
################################################################################

# imports
import argparse
from .scims import *
from .scims.determine_sex import *
import os.path
from scipy import stats as st
import pysam

def scims():
    """
    Main function for the SCiMS program
    """

    parser = argparse.ArgumentParser(
        description="Sex Calling from Metagenomic Sequences"
    )

    parser.add_argument("-v", "--version", action="version", version="scims 0.0.1")

    ######################################
    # determine-sex
    ######################################
    parser.add_argument(
        "--x",
        dest="homogametic",
        help="ID of homogametic sex chromosome (ex. X)",
    )
    parser.add_argument(
        "--y",
        dest="heterogametic",
        help="ID of heterogametic sex chromosome (ex. Y)"
    )

    parser.add_argument(
        "--bowtie-index", dest="index", help="Path to Bowtie2 index"
    )

    parser.add_argument(
        "--ref-genome", dest="ref", help="Reference genome in FASTA or FASTA.gz format"
    )

    parser.add_argument(
        "--scaffold-names", dest="scaffold", help="Path to file with scaffold names"
    )

    parser.add_argument(
        "-1",
        dest="forward",
        help="Forward reads in FASTA or FASTQ format (PE mode only)")

    parser.add_argument("-2", dest="reverse", help="Reverse reads in FASTA or FASTQ format (PE mode only)")

    parser.add_argument("-s", dest="single", help="Single-end reads in FASTA or FASTQ format (SE mode only)" )

    parser.add_argument(
        "--o",
        dest="output",
        help="Output plot prefix",
        required="True"
    )
    parser.add_argument(
        "--t",
       dest="threads", help="Number of threads to use"
    )

    parser.add_argument(
        "--from-sam", dest="from_sam", help="Use sam file instead of preforming the alignment"
    )

    parser.add_argument(
        "--get-scaffold-lengths", dest="lengths", help="Determine scaffold lengths"
    )

    parser.add_argument(
        "--scaffold-lengths", dest="scaf_lengths", help="Path to file with scaffold lengths"
    )


    args = parser.parse_args()

    ############################################################################
    # Verify inputs
    ############################################################################

    if args.ref:
        if os.path.exists(args.ref):
            print("Using reference genome")
        else:
            print("Reference genome does not exist")
            exit()


    # Check to see if the user wants to jump in from the BAM file 
    if args.lengths:
        print("Getting scaffold lengths")
        scaffolds = read_scaffold_names(args.scaffold)
        lengths = determine_chrom_lengths(args.ref, scaffolds)
        with open(args.output, "w") as f:
            for chrom, length in lengths.items():
                f.write(f"{chrom}\t{length}\n")
        exit()

    if args.from_sam:
        alignment = args.from_sam
        if os.path.exists(args.from_sam):
            print("Using SAM file")
        else:
            print("SAM file does not exist")
            exit()

    # Preform the alignment if the "--from-bam" flag is not specified 
    else:
        # Make sure we have all the flags we need
        if not args.index:
            print("Please specify a Bowtie2 index and reference genome")
            exit()

        if not args.ref:
            if not args.scaf_lengths:
                print("Please specify a reference genome")
                exit()

        if args.forward and args.reverse:
            #Paired-end mode
            print("Aligning sequences to {}:".format(args.index))
            print("Forward reads: {}".format(args.forward))
            print("Reverse reads: {}".format(args.reverse))
            alignment = paired_reads_with_bowtie2(args.index, args.forward, args.reverse, "local", args.threads)

        elif args.single:
            #Single-end mode
            print("Aligning {} to {}:".format(args.single, args.index))
            alignment = single_reads_with_bowtie2(args.single, args.index, "local", args.threads)

        else:
            print("Please specify a FASTQ file for alignment")
            exit()

    if args.scaffold:
            
        # If the user supplies a file with the scaffold IDs to be used
        scaffolds = read_scaffold_names(args.scaffold)

        if args.scaf_lengths:
            # If the user supplies a file with the scaffold lengths
            lengths = read_scaffold_lengths(args.scaf_lengths)

        else:
            # Calculate the lengths of the scaffolds 
            print("Getting scaffold lengths for {}".format(args.ref))
            lengths = determine_chrom_lengths(args.ref, scaffolds)

        # # Extract the proper alignments from the BAM file 
        # print("Extracting propper alignments:")
        # counts = count_chrom_alignments(alignment, scaffolds)

    else:
        pass
        # # Extract the proper alignments from the BAM file 
        # print("Extracting propper alignments:")
        # counts = count_chrom_alignments(alignment)

    # Print the results 
    # print("\tA total of {} alignments met the criteria".format(str(count_seqs(counts))))

    # Determine the coverages 
    print("Determining coverages")
    window = lengths[[length for length in lengths][0]] // 100
    chrom_window_dict = get_chrom_window_coverage(alignment, lengths, window, scaffolds)

    big_list = []
    for chrom in chrom_window_dict:
        if chrom != args.heterogametic:
            big_list = big_list + list(chrom_window_dict[chrom])

    plt.style.context('seaborn')
    plt.hist([f for f in list(chrom_window_dict[args.homogametic])], bins=10, label=args.homogametic, density=True)
    plt.hist([f for f in big_list], bins=10, label=args.heterogametic, density=True)
    plt.xlabel("Coverage")
    plt.ylabel("Density")
    plt.title(args.output)
    plt.savefig(args.output + ".png", dpi=300, bbox_inches="tight")

    stat = (st.kstest([f for f in list(chrom_window_dict[args.homogametic])], [f for f in big_list]))
    print(f"P-Value:{stat[1]}")
    print(f"Distribution output to: {args.output}.png")
    print("Thank you for using SCiMS!")