# Author: Kobie Kirven
# Davenport Lab - Penn State University
# Date: 9-2-2021

# imports
import argparse
from .scims import *
from .scims.determine_sex import *
import time
import gzip


def scims():
    """
    Main function for the scims program
    """

    parser = argparse.ArgumentParser(
        description="Sex Calling for Metagenomic Sequences"
    )

    parser.add_argument("-v", "--version", action="version", version="scims 0.0.1")

    ######################################
    # determine-sex
    ######################################

    parser.add_argument(
        "-i", dest="index", help="Path to Bowtie2 index",
        required=True
    )

    parser.add_argument(
        "-r", dest="ref", help="Reference genome in FASTA or FASTA.gz format",
        required=True
    )

    parser.add_argument(
        "-sca", dest="scaffold", help="Path to file with scaffold names",
        required=True
    )

    parser.add_argument(
        "-1",
        dest="forward",
        help="Forward reads in FASTA or FASTQ format (PE mode only)")

    parser.add_argument("-2", dest="reverse", help="Reverse reads in FASTA or FASTQ format (PE mode only)")

    parser.add_argument("-s", dest="single", help="Single-end reads in FASTA or FASTQ format (SE mode only)" )

    parser.add_argument(
        "-t",
       dest="threads", help="Number of threads to use"
    )

    parser.add_argument(
        "-hom",
        dest="homogametic",
        help="ID of homogametic sex chromosome (ex. X)",
        required="True"
    )
    parser.add_argument(
        "-het",
        dest="heterogametic",
        help="ID of heterogametic sex chromesome (ex. Y)"
    )

    parser.add_argument(
        "-o",
        dest="output",
        help="Output plot prefix",
        required="True"
    )

    args = parser.parse_args()

    if type(args.forward) == str:
        #Paired-end mode
        print("Aligning sequences to {}:".format(args.index))
        alignment = paired_reads_with_bowtie2(args.index, args.forward, args.reverse, "local", args.threads)


    elif args.single:
        #Single-end mode
        print("Aligning {} to {}:".format(args.single, args.index))
        alignment = single_reads_with_bowtie2(args.single, args.index, "local", args.threads)
        print(alignment)

    print("Getting scaffold lengths for {}".format(args.ref))
    lengths = determine_chrom_lengths(args.ref, args.scaffold)

    if args.scaffold:
        scaffolds = read_scaffold_names(args.scaffold)
        print("Extracting propper alignments:")
        counts = count_chrom_alignments(alignment, scaffolds)

    else:
        print("Extracting propper alignments:")
        counts = count_chrom_alignments(alignment)

    print("\tA total of {} alignments met the criteria".format(str(count_seqs(counts))))
    normalized_length_counts = normalize_by_chrom_lengths(counts, lengths)
    with open(args.output + ".txt", "w") as fn:
        for count in  normalized_length_counts:
            fn.write("{}\t{}\n".format(count, normalized_length_counts[count]))
    compare = compare_to_homogametic(normalized_length_counts, args.homogametic)
    generate_plot(compare, args.output, args.homogametic, args.heterogametic)
    stats = stats_test(compare, args.homogametic, args.heterogametic)
    print("The average {}:{} ratio for the file was {}".format(args.heterogametic, args.homogametic, stats))
    delete_temp_file_list()



