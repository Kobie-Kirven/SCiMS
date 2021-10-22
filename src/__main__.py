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
        "-i", dest="index", help="Name of bowtie2 and bwa index",
        required=True
    )

    parser.add_argument(
        "-r", dest="ref", help="Reference genome in FASTA or FASTA.gz format",
        required=True
    )

    parser.add_argument(
        "-1",
        dest="forward",
        help="Forward reads in fasta or fastq format")

    parser.add_argument("-2", dest="reverse", help="Reverse reads in fasta or fastq format")

    parser.add_argument("-s", dest="single", help="Single-end reads in fasta or fastq format" )

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
        help="Output plot name",
        required="True"
    )

    args = parser.parse_args()

    if args.forward is not None and args.reverse is not None:
        #Paired-end mode
        print("Aligning sequences to {}:".format(args.index))
        alignment = paired_reads_with_bowtie2(args.index, args.forward, args.reverse, "local", args.threads)

        print("Getting scaffold lengths for {}".format(args.reference))
        lengths = determine_chrom_lengths(args.ref)

    elif args.single:
        #Single-end mode
        print("Aligning {} to {}:".format(args.single, args.index))
        alignment = paired_reads_with_bowtie2(args.index, args.forward, args.reverse, "local", args.threads)


    print("Extracting propper alignments:")
    counts = count_chrom_alignments(alignment)
    print("\tA total of {} alignments met the criteria")
    print(count_seqs(counts))
    normalized_length_counts = normalize_by_chrom_lengths(counts, lengths)
    compare = compare_to_homogametic(normalized_length_counts, args.homogametic)
    generate_plot(compare, args.output, args.homogametic, args.heterogametic)
    stats = stats_test(compare, args.homogametic, args.heterogametic)
    print("The average {}:{} ratio for the file was {}".format(args.heterogametic, args.homogametic, stats))
    delete_temp_file_list()



