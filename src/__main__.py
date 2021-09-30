# Author: Kobie Kirven
# Davenport Lab - Penn State University
# Date: 9-2-2021


# imports
import argparse
from .scims import *
from .scims.determine_sex import *
import time


def scims():
    """Main function for the scims program"""

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
        help="ID of homogametic sex chromesome (ex. X)",
        required="True"
    )
    determine_sex_parser.add_argument(
        "-het",
        "--heterogametic",
        dest="heterogametic",
        help="ID of heterogametic sex chromesome (ex. Y)",
        required="True"
    )

    args = parser.parse_args()

    if args.command == "build-index":
        buildBwaIndex(args.reference, args.output)
        buildBowtie2Index(args.reference, args.output)

    elif args.command == "determine-sex":
        print("Now aligning reads with BWA..")
        timeStart = time.time()
        bwa = align_with_bwa(args.index, args.forward, args.reverse, args.threads)
        timeStop = time.time()
        print("BWA finished in {:.2f} seconds\n".format(timeStop - timeStart))

        human = get_human_sequences(bwa)
        sam = get_sex_sequences(human, args.homogametic, args.heterogametic)
        uniqueIds = get_unique_ids_in_sam(sam)
        files = get_fastq_reads_in_sam(uniqueIds, args.forward, args.reverse)

        print("Now merging files with Flash...")
        timeStart = time.time()
        merged = merge_with_flash(files[0], files[1])
        timeStop = time.time()
        print("Flash finished in {:.2f} seconds\n".format(timeStop - timeStart))

        print("Now aligning merged reads with Bowtie2...")
        timeStart = time.time()
        mergeAligned = align_merged_with_bowtie2(args.index, merged[0])
        timeStop = time.time()
        print("Bowtie2 finished in {:.2f} seconds\n".format(timeStop - timeStart))

        print("Now aligning unmerged reads with Bowtie2...")
        timeStart = time.time()
        unmergeAligned = paired_reads_with_bowtie2(args.index, merged[1], merged[2])
        timeStop = time.time()
        print("Bowtie2 finished in {:.2f} seconds\n".format(timeStop - timeStart))

        print("Preforming read rescue...")
        mergeDf = filter_merged_sam(mergeAligned)

        unmergedDf = read_rescue_unmerged(unmergeAligned)
        combDf = combine_df(mergeDf, unmergedDf)
        filtSam = get_sex_sequences(
            human, args.homogametic, args.heterogametic, is_filter=False
        )
        update = read_rescue_update(combDf, filtSam, 50)
        chromCounts = count_chrom(update, args.homogametic, args.heterogametic)
        stats = calculate_stats(chromCounts[0], chromCounts[1])
        print(
            "The proportion of {} reads to {} reads is:".format(
                args.heterogametic, args.homogametic
            )
        )
        print("{} \u00B1 {} (95% CI)".format(stats[0], stats[1]))
        print("Thank you for using SCiMS!")
        delete_temp_file_list()
