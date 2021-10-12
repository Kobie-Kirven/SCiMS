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

    subparser = parser.add_subparsers(dest="command")

    parser.add_argument("-v", "--version", action="version", version="scims 0.0.1")

    #######################################
    # Arguments for build-index
    #######################################
    build_index_parser = subparser.add_parser(
        "build-index",
        help="Build indexes for BWA and Bowtie2",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    build_index_parser.add_argument(
        "-r",
        "--reference",
        dest="reference",
        help="Host reference genome in FASTA or FASTQ format",
        required=True
    )

    build_index_parser.add_argument(
        "-o", "--output", dest="output", help="Name of output index",
        required=True
    )

    ######################################
    # determine-sex
    ######################################
    determine_sex_parser = subparser.add_parser(
        "determine-sex",
        help="determine sex of sample",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    determine_sex_parser.add_argument(
        "-i", "--index-name", dest="index", help="Name of bowtie2 and bwa index",
        required=True
    )

    determine_sex_parser.add_argument(
        "-r", "--reference-genome", dest="ref", help="Reference genome in FASTA or FASTA.gz format",
        required=True
    )

    determine_sex_parser.add_argument(
        "-1",
        "--forward-reads",
        dest="forward",
        help="Forward reads in fasta or fastq format",
        required=True
    )

    determine_sex_parser.add_argument(
        "-2",
        "--reverse-reads",
        dest="reverse",
        help="Reverse reads in fasta or fastq format",
        required=True
    )

    determine_sex_parser.add_argument(
        "-t", "--number-of-threads", dest="threads", help="Number of threads to use",
        required=True
    )

    determine_sex_parser.add_argument(
        "-hom",
        "--homogametic",
        dest="homogametic",
        help="ID of homogametic sex chromosome (ex. X)",
        required="True"
    )
    determine_sex_parser.add_argument(
        "-het",
        "--heterogametic",
        dest="heterogametic",
        help="ID of heterogametic sex chromesome (ex. Y)"
    )

    determine_sex_parser.add_argument(
        "-o",
        "--output",
        dest="output",
        help="Output plot name",
        required="True"
    )

    args = parser.parse_args()

    if args.command == "build-index":
        buildBwaIndex(args.ref, args.output)
        buildBowtie2Index(args.ref, args.output)

    elif args.command == "determine-sex":
        alignment = paired_reads_with_bowtie2(args.index, args.forward, args.reverse, "local", threads=1)
        lengths = determine_chrom_lengths(args.ref)
        counts = count_chrom_alignments(alignment)
        normalized_length_counts = normalize_by_chrom_lengths(counts, lengths)
        compare = compare_to_homogametic(normalized_length_counts, args.homogametic)
        generate_plot(compare, args.output, args.homogametic, args.heterogametic)
        stats = stats_test(compare, args.homogametic,args.heterogametic)
        if stats[0].pvalue > stats[1].pvalue:
            print("Male (p=" + str(stats[0].pvalue) + ") two-tailed t-test")
        elif stats[1].pvalue > stats[0].pvalue:
            print("Female (p=" + str(stats[1].pvalue) + ") two-tailed t-test")
        delete_temp_file_list()



