# Author: Kobie Kirven
# Davenport Lab - Penn State University
# Date: 9-2-2021


# imports
import argparse
from .scims import (BuildIndex,TestFile)
from .scims.determine_sex import GetProportion


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
        help="ID of heterogametic sex chromesome (ex. Y)",)


    args = parser.parse_args()

    if args.command == "build-index":
        index = BuildIndex(
            args.reference, args.homogametic, args.heterogametic, args.output
        )
        index.buildBwaIndex()
        index.buildBowtie2Index()

    elif args.command == "determine-sex":
        sample = GetProportion(args.index, args.forward, args.reverse, args.threads)
        sample.getHumanSequences()
        sample.getSexSequences(args.homogametic,args.heterogametic,"reads.sam")
        sample.getFastqReadsInSam()
        sample.mergeWithFlash()
        sample.alignWithBowtie2(input1="flash.extendedFrags.fastq", merged=True)
        sample.alignWithBowtie2(input1="flash.notCombined_1.fastq", input2="flash.notCombined_2.fastq",
         merged=False)

