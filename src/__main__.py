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
        help="Forward reads in fasta or fastq format",
        required=True)
    parser.add_argument("-2", dest="reverse", help="Reverse reads in fasta or fastq format", required=True)


    parser.add_argument(
        "-t",
       dest="threads", help="Number of threads to use",
        required=True
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

    fast_align = paired_reads_with_bowtie2(args.index, args.forward, args.reverse, "very-fast-local", args.threads)
    human_reads = get_human_sequences(fast_align)
    ids = get_unique_ids_in_sam(human_reads)
    reads = get_fastq_reads_in_sam(ids, args.forward, args.reverse)
    alignment = paired_reads_with_bowtie2(args.index, reads[0], reads[1], "local", args.threads)
    lengths = determine_chrom_lengths(args.ref)
    counts = count_chrom_alignments(alignment)
    print(count_seqs(counts))
    normalized_length_counts = normalize_by_chrom_lengths(counts, lengths)
    compare = compare_to_homogametic(normalized_length_counts, args.homogametic)
    generate_plot(compare, args.output, args.homogametic, args.heterogametic)
    stats = stats_test(compare, args.homogametic,args.heterogametic)
    if stats[0].pvalue > stats[1].pvalue:
        print("Male")
    elif stats[1].pvalue > stats[0].pvalue:
        print("Female")
    delete_temp_file_list()



