# Author: Kobie Kirven
# Davenport Lab - Penn State University
# Date: 9-2-2021


# imports
import argparse
from .scims import BuildIndex
from .scims.determine_sex import *


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
        help="Build indexes for bwa and bowtie2",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    build_index_parser.add_argument(
        "-r",
        "--reference",
        dest="reference",
        help="Host reference genome in FASTA or FASTQ format",
    )
    build_index_parser.add_argument(
        "-hom",
        "--homogametic",
        dest="homogametic",
        help="ID of homogametic sex chromesome (ex. X)",
    )
    build_index_parser.add_argument(
        "-het",
        "--heterogametic",
        dest="heterogametic",
        help="ID of heterogametic sex chromesome (ex. Y)",
    )
    build_index_parser.add_argument(
        "-o", "--output", dest="output", help="Name of output index"
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
        "-i",
        "--index-name",
        dest="index",
        help="Name of bowtie2 and bwa index",
    )

    determine_sex_parser.add_argument(
        "-1",
        "--forward-reads",
        dest="forward",
        help="Forward reads in fasta or fastq format",
    )

    determine_sex_parser.add_argument(
        "-2",
        "--reverse-reads",
        dest="reverse",
        help="Reverse reads in fasta or fastq format",
    )

    determine_sex_parser.add_argument(
        "-t",
        "--number-of-threads",
        dest="threads",
        help="Number of threads to use",
    )

    determine_sex_parser.add_argument(
        "-hom",
        "--homogametic",
        dest="homogametic",
        help="ID of homogametic sex chromesome (ex. X)",
    )
    determine_sex_parser.add_argument(
        "-het",
        "--heterogametic",
        dest="heterogametic",
        help="ID of heterogametic sex chromesome (ex. Y)",
    )

    args = parser.parse_args()

    if args.command == "build-index":
        index = BuildIndex(
            args.reference, args.homogametic, args.heterogametic, args.output
        )
        index.buildBwaIndex()
        index.buildBowtie2Index()

    elif args.command == "determine-sex":
        bwa = alignWithBwa(args.index, args.forward, args.reverse, args.threads)
        human = getHumanSequences(bwa)
        sam = getSexSequences(human, args.homogametic, args.heterogametic)
        uniqueIds = getUniqueIdsInSam(sam)
        files = getFastqReadsInSam(uniqueIds, args.forward, args.reverse)
        merged = mergeWithFlash(files[0], files[1])
        mergeAligned = alignMergedWithBowtie2(args.index, merged[0])
        unmergeAligned = pairedReadsWithBowtie2(args.index, merged[1], merged[2])
        mergeDf = readRescueMerged(mergeAligned)
        unmergedDf = readRescueUnmerged(unmergeAligned)
        combDf = combineDf(mergeDf,unmergedDf)
        filtSam = getSexSequences(human, args.homogametic, args.heterogametic, isFilter=False)
        update = readRescueUpdate(combDf,filtSam, 50)
        chromCounts = countChrom(update, args.homogametic,args.heterogametic)
        stats = calculateStats(chromCounts[0], chromCounts[1])
        print("The proportion of {} reads to {} reads is:".format(args.heterogametic, args.homogametic))
        print('{} u"\u00B1" {}'.format(stats[0], stats[1]))
        print("Thank you for using SCiMS!")
        deleteTempFileList()


