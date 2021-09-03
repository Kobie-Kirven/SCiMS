# Author: Kobie Kirven
# Davenport Lab - Penn State University
# Date: 9-2-2021


# imports
import argparse
from .scims import BuildIndex


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
        help="get sequences of assembled chromosomes",
        description="Get the FASTA file of assembled chromosomes.",
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
    args = parser.parse_args()

    if args.command == "build-index":
        index = BuildIndex(
            args.reference, args.homogametic, args.heterogametic, args.output
        )
        index.buildBwaIndex()
        index.buildBowtie2Index()
