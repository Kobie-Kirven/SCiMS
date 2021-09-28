# Author: Kobie Kirven, Kyle McGovern
# Davenport Lab - Pennsylvania State University
# Date: 9-2-2021

# Imports
import os
from Bio import SeqIO
from .sam_files import *
import subprocess
import pandas as pd
import numpy as np
import tempfile
import math

filesList = []


def tempFilesList(fileName):
    """
    Adds the temporary file to the list of temporary files
    
    Args: fileName(str): path to the temporary file
    """
    if type(fileName) in [int, list]:
        raise TypeError
    elif fileName in filesList:
        pass
    else:
        filesList.append(fileName)


def deleteTempFileList():
    """
    Delete the list of temporary files
    """
    try:
        for file in filesList:
            os.unlink(file)
    finally:
        pass


def alignWithBwa(index, forwardReads, reverseReads, threads):
    """
        Align paired-end reads with BWA

        Parameters:
            index (str): prefix for BWA index
            forwardReads (str): path to file with forward FASTQ reads
            reverseReads (str): path to file with reverse FASTQ reads
            threads (int): number of threads to use

        Returns:
            f.name (str): Path to output SAM file
    """
    with tempfile.NamedTemporaryFile(delete=False) as f:
        subprocess.run(
            ["bwa", "mem", "-t", threads, index, forwardReads, reverseReads],
            stdout=f,
            stderr=subprocess.DEVNULL,
        )
        tempFilesList(f.name)
        return f.name


def getHumanSequences(inputSam):
    """
    Get the sequences from BWA alignment output
    that map to the human genome

    Parameters:
        inputSam (str): Path to input SAM file

    Returns:
        f.name (str): Path to output BAM file
    """
    with tempfile.NamedTemporaryFile(delete=False) as f:
        subprocess.run(
            ["samtools", "view", "-F", "4", "-F", "8", "-b", inputSam],
            stdout=f,
            stderr=subprocess.DEVNULL,
        )
        tempFilesList(f.name)
    return f.name


def getSexSequences(
    inputBam, homogamete, heterogamete, isFilter=True, mapqMin=50, tlen=75
):
    """
    Get the sequences that map to the sex chromosomes

    Parameters:
        inputBam (str): Path to input BAM file
        homogamete (str): FASTA ID for homogametic element
        heterogamete (str): FASTA ID for heterogametic element
        isFilter (bool): Filter reads based on quality thresholds
        mapqMin (int): Minimum mapq score
        tlen (int): Template length

    Returns:
        g.name (str): Path to SAM file with reads that mapped to sex chromosomes

    """
    with tempfile.NamedTemporaryFile(delete=False) as f:
        smallSam = subprocess.run(
            ["samtools", "view", "-F", "2048", "-F", "256", inputBam], stdout=f
        )
        name = f.name
        tempFilesList(f.name)

    if isFilter:
        with tempfile.NamedTemporaryFile(delete=False) as g:
            for line in open(f.name):
                rec = ReadSamLine(line)
                rec.getFeatures()
                if rec.mapq < mapqMin:
                    if rec.rnam == homogamete or rec.rnam == heterogamete:
                        if rec.tlen > tlen or rec.tlen < (tlen * -1):
                            g.write(line.encode())
            tempFilesList(g.name)
            return g.name
    else:
        with tempfile.NamedTemporaryFile(delete=False) as g:
            for line in open(f.name):
                rec = ReadSamLine(line)
                rec.getFeatures()
                if rec.rnam == homogamete or rec.rnam == heterogamete:
                    if rec.tlen > tlen or rec.tlen < (tlen * -1):
                        g.write(line.encode())
            tempFilesList(g.name)
            return g.name


def getUniqueIdsInSam(inputSam):
    """
    Get the unique query names from the SAM file

    Parameters:
        inputSam (str): Path to input SAM file

    returns:
        (list): Set of unique FASTA IDs
    """
    readIds = []
    with open(inputSam) as fn:
        for rec in fn:
            line = ReadSamLine(rec)
            line.getFeatures()
            readIds.append(line.query)
    return set(readIds)


def getFastqReadsInSam(readIds, forwardReads, reverseReads):
    """
    Get FASTQ reads for sequences in SAM file

    Parameters:
        readIds (list): List of FASTA IDs
        forwardReads (str): Path to forward FASTQ reads
        reverseReads (str): Path to reverse FASTQ reads

    Returns:
        outFastq1.name (str): Path to output forward FASTQ file
        outFastq2.name (str): Path to output reverse FASTQ file
    """

    # Open output files
    outFastq1, outFastq2 = (
        tempfile.NamedTemporaryFile(delete=False),
        tempfile.NamedTemporaryFile(delete=False),
    )

    # Get whether the file is FASTA or FASTQ
    fileType = fastaOrFastq(forwardReads)

    # Open the input read files and get those FASTQ sequences from the SAM file
    if isFileGzip(forwardReads):
        forwardReads = gzip.open(forwardReads, "rt")
        reverseReads = gzip.open(reverseReads, "rt")

    for rec1 in SeqIO.parse(forwardReads, fileType):
        if str(rec1.id)[:-2] in readIds:
            outFastq1.write(("@" + str(rec1.description)[:-2] + "\n").encode())
            outFastq1.write((str(rec1.seq) + "\n").encode())
            outFastq1.write("+\n".encode())
            quality = ""
            for qual in rec1.letter_annotations["phred_quality"]:
                quality = quality + str(chr(qual + 33))
            outFastq1.write((str(quality) + "\n").encode())

    for rec2 in SeqIO.parse(reverseReads, fileType):
        if str(rec2.id)[:-2] in readIds:
            outFastq2.write(("@" + str(rec2.description)[:-2] + "\n").encode())
            outFastq2.write((str(rec2.seq) + "\n").encode())
            outFastq2.write("+\n".encode())
            quality = ""
            for qual in rec2.letter_annotations["phred_quality"]:
                quality = quality + str(chr(qual + 33))
            outFastq2.write((str(quality) + "\n").encode())

    outFastq1.close()
    outFastq2.close()
    tempFilesList(outFastq1.name)
    tempFilesList(outFastq2.name)
    return outFastq1.name, outFastq2.name


def mergeWithFlash(pair1, pair2, mismatchRatio=0.1, maxOverlap=150, minOverlap=40):
    """
    Merge paired-end reads with FLASH

    Parameters:
        pair1 (str): Path to forward FASTQ file
        pair2 (str): Path to reverse FASTQ file
        mismatchRatio (float): Maximum mismatch ratio
        maxOverlap (int): Maximum overlap for read merger
        minOverlap (int): Minimum overlap for read merger

    Returns:
        merged (str): Path to merged FASTQ file
        unmerged1 (str): Path to forward, unmerged FASTQ file
        unmerged2 (str): Path to reverse, unmerged FASTQ file
    """
    tempdir = tempfile.TemporaryDirectory()
    subprocess.run(
        [
            "flash",
            "-x",
            str(mismatchRatio),
            "-M",
            str(maxOverlap),
            "-m",
            str(minOverlap),
            pair1,
            pair2,
            "-d",
            tempdir.name,
        ],
        stdout=subprocess.DEVNULL,
    )
    merged = convertToTemp(tempdir.name + "/out.extendedFrags.fastq")
    unmerged1 = convertToTemp(tempdir.name + "/out.notCombined_1.fastq")
    unmerged2 = convertToTemp(tempdir.name + "/out.notCombined_2.fastq")
    return merged, unmerged1, unmerged2


def convertToTemp(oldFile):
    """
    Convert a file into a temporary file

    Parameters:
        oldFile (str): Path to input file

    Returns:
        f.name (str): Path to output temporary file
    """
    with open(oldFile) as fn:
        with tempfile.NamedTemporaryFile(delete=False) as f:
            for line in fn:
                f.write(line.encode())
            return f.name


def alignMergedWithBowtie2(index, merged, threads=1, noOfAlignPerRead=50):
    """
    Align merged reads with bowtie2

    Parameters:
        index (str): Prefix for Bowite2 index
        merged (str): Path to FLASH merged, FASTQ reads
        threads (int): Number of threads to use
        noOfAlignPerRead (int): Maximum number of alignments per-read

    Returns:
        f.name (str): Path to output SAM file
    """
    with tempfile.NamedTemporaryFile(delete=False) as f:
        subprocess.run(
            [
                "bowtie2",
                "-p",
                str(threads),
                "-k",
                str(noOfAlignPerRead),
                "--local",
                "-x",
                index,
                "-U",
                merged,
            ],
            stdout=f,
            stderr=subprocess.DEVNULL,
        )
        tempFilesList(f.name)
        return f.name


def pairedReadsWithBowtie2(
    index, forwardReads, reverseReads, threads=1, noOfAlignPerRead=50
):
    """
    Align paired reads with bowtie2

    Parameters:
        index (str): Prefix for Bowite2 index
        forwardReads (str): Path to forward FASTQ reads
        reverseReads (str): Path to reverse FASTq reads
        threads (int): Number of threads to use
        noOfAlignPerRead (int): Maximum number of alignments per-read

    Returns:
        f.name (str): Path to output SAM file
    """
    with tempfile.NamedTemporaryFile(delete=False) as f:
        subprocess.run(
            [
                "bowtie2",
                "-p",
                str(threads),
                "-k",
                str(noOfAlignPerRead),
                "--local",
                "-x",
                index,
                "-1",
                forwardReads,
                "-2",
                reverseReads,
            ],
            stdout=f,
            stderr=subprocess.DEVNULL,
        )
        return f.name


def getReadDataFromMergedSam(samFile):
    """
    Get data from the SAM file of merged reads

    Parameters:
        samFile (str): Path to input SAM file

    Returns:
        fastqData (list): A list of dictionaries where each element
                            corresponds to each line in the SAM file
    """
    fastqData = []
    for rec in ParseSam(samFile):
        if rec.cigar != "*":
            numMatches = extractFromCigar("M", rec.cigar)
            numDeletions = extractFromCigar("D", rec.cigar)

            if rec.flag != 4 and rec.flag != 8:
                readData = {
                    "readName": rec.query,
                    "chrom": rec.rnam,
                    "alignmentScore": rec.align_score,
                    "numMismatches": rec.mismatches,
                    "insertSize": None,
                    "pos_1": rec.pos,
                    "pos_2": -1,
                    "numMatches": numMatches,
                    "numDeletions": numDeletions,
                }

                fastqData.append(readData)
    return fastqData


def getReadDataFromUnmergedSam(unmergedSam):
    """
    Get data from SAM file of unmerged reads

    Parameters:
        unmergedSam (str): Path to input SAM file

    Returns:
        fastqData (list): List of dictionaries where each element
                            in the list corresponds to a line in
                            the SAM file
    """
    fastqData = []
    for rec in ParseSam(unmergedSam):
        if rec.cigar != "*":
            numMatches = extractFromCigar("M", rec.cigar)
            numDeletions = extractFromCigar("D", rec.cigar)
            readData = {
                "readName": rec.query,
                "samFlag": rec.flag,
                "chrom": rec.rnam,
                "alignmentScore": rec.align_score,
                "numMismatches": rec.mismatches,
                "mateStart": rec.pnext,
                "insertSize": rec.tlen,
                "pos": rec.pos,
                "numMatches": numMatches,
                "numDeletions": numDeletions,
            }
            fastqData.append(readData)
    return fastqData


def readRescueMerged(
    mergedSam,
    diffFromMaxPercent=0.025,
    minMatch=74,
    mismatchRatio=0.025,
    numDeletions=3,
):
    """
    Preform 'read rescue' on the merged reads

    Parameters:
         mergedSam (str): Path to SAM file of merged reads
         diffFromMaxPercent (float): Maximum difference in alignment quality between
                                    a given read and the best alignment for that read
        minMatch (int): Minimum number of bases matching between the query and reference
        mismatchRatio (float): Maximum ratio of matches to mismatches between query
                                and reference
        numDeletions (float): Maximum number of deletions

    Returns:
        df (pandas data frame): Data frame with info for sequences that survived read rescue
    """
    fastqData = getReadDataFromMergedSam(mergedSam)

    if fastqData:
        # Transform the fastqData into a pandas dataframe
        df = pd.DataFrame(fastqData)

        # Filter reads that arent within a certain quality percentage
        # from the maximum scoring alignment for each read
        df["maxScoreReadName"] = df.groupby("readName")["alignmentScore"].transform(max)
        df["diffPercentMaxRound"] = np.round(
            df["maxScoreReadName"] * diffFromMaxPercent
        )
        df["scoreMaxDiff"] = df["maxScoreReadName"] - df["alignmentScore"]
        df = df[df["scoreMaxDiff"] <= df["diffPercentMaxRound"]]

        # Filter out reads with less than a certain number of bases matching.
        df = df[df["numMatches"] > minMatch]

        # Filter out reads with too many mismatches (more than 1 mismatch every 40 bp).
        df["mismatchRatio"] = df["numMismatches"].astype(float) / df[
            "numMatches"
        ].astype(float)
        df = df[df["mismatchRatio"] < mismatchRatio]

        # Filter out reads with too many deletions.
        df = df[df["numDeletions"] < numDeletions]

        # Finally keep reads that have only one annotated chromosome.
        mapper = df.groupby("readName")["chrom"].nunique().to_dict()
        df["uniqueChroms"] = df["readName"].map(mapper)
        df = df[df["uniqueChroms"] == 1]
        return df
    else:
        return False


def readRescueUnmerged(
    unmergedSam, diffFromMaxPercent=0.025, minMatch=100, numDeletions=3
):
    """"
    Preform 'read rescue' for unmerged (paired-end) sequences

    Parameters:
        unmergedSam (str): Path to SAM file of unmerged reads
        diffFromMaxPercent (float): Maximum difference in alignment quality between
                                    a given read and the best alignment for that read
        minMatch (int): Minimum number of bases matching between the query and reference
        mismatchRatio (float): Maximum ratio of matches to mismatches between query
                                and reference
        numDeletions (float): Maximum number of deletions

    Returns:
        df_pair (pandas data frame): Data frame with data for reads that survived
                                    read rescue

    """
    read_pair = []
    fastq_data = []
    fastqData = getReadDataFromUnmergedSam(unmergedSam)

    for readData in fastqData:
        read_pair.append(readData)
        if read_pair:
            if (
                abs(readData["insertSize"]) != abs(read_pair[0]["insertSize"])
                or read_pair[0]["chrom"] != readData["chrom"]
                or readData["readName"] != read_pair[0]["readName"]
            ):
                continue
            read_pair.append(readData)
            total_num_matches = read_pair[0]["numMatches"] + read_pair[1]["numMatches"]
            total_num_mismatches = (
                read_pair[0]["numMismatches"] + read_pair[1]["numMismatches"]
            )
            total_num_deletions = (
                read_pair[0]["numDeletions"] + read_pair[1]["numDeletions"]
            )
            total_alignment_score = (
                read_pair[0]["alignmentScore"] + read_pair[1]["alignmentScore"]
            )
            fastq_data.append(
                {
                    "readName": readData["readName"],
                    "numMismatches": total_num_mismatches,
                    "chrom": readData["chrom"],
                    "numMatches": total_num_matches,
                    "numDeletions": total_num_deletions,
                    "insertSize": abs(readData["insertSize"]),
                    "alignmentScore": total_alignment_score,
                    "pos_1": read_pair[0]["pos"],
                    "pos_2": read_pair[1]["pos"],
                }
            )
            read_pair = []

        else:
            if readData["samFlag"] & 2 == 2 and readData["samFlag"] & 64 == 64:
                read_pair.append(readData)

    df_pair = pd.DataFrame(fastq_data)

    # Filter reads with big insert size (gt 999)
    df_pair = df_pair[df_pair["insertSize"] < 1000]

    # Filter out reads less than 2.5% max score.
    df_pair["maxScoreReadName"] = df_pair.groupby("readName")[
        "alignmentScore"
    ].transform(max)
    df_pair["diffPercentMaxRound"] = np.round(
        df_pair["maxScoreReadName"] * diffFromMaxPercent
    )
    df_pair["scoreMaxDiff"] = df_pair["maxScoreReadName"] - df_pair["alignmentScore"]
    df_pair = df_pair[df_pair["scoreMaxDiff"] <= df_pair["diffPercentMaxRound"]]

    # Filter out reads with less than 100 bp matching.
    df_pair = df_pair[df_pair["numMatches"] >= minMatch]

    # Filter out reads with too many mismatches (more than 1 mismatch every 40 bp).
    df_pair["mismatchRatio"] = df_pair["numMismatches"].astype(float) / df_pair[
        "numMatches"
    ].astype(float)
    df_pair = df_pair[df_pair["mismatchRatio"] < 0.025]

    # Filter out reads with too many deletions.
    df_pair = df_pair[df_pair["numDeletions"] <= numDeletions]

    # Finally keep reads that have only one annotated chromosome.
    mapper = df_pair.groupby("readName")["chrom"].nunique().to_dict()
    df_pair["uniqueChroms"] = df_pair["readName"].map(mapper)
    df_pair = df_pair[df_pair["uniqueChroms"] == 1]
    return df_pair


def combineDf(df1, df2):
    """
    Combine 2 pandas data frames

    Parameters:
        df1 (pandas data frame): Data frame 1
        df2 (pandas data frame): Data frame 2

    Returns:
        (pandas df): Concatenated data frame

    """
    return pd.concat([df1, df2])


def bam2sam(bamFile):
    """
    Convert BAM files to SAM files

    Parameters:
        bamFile (str): path to BAM file

    Returns:
        f.name (str): path to SAM file
    """
    with tempfile.NamedTemporaryFile(delete=False) as f:
        subprocess.run(["samtools", "view", bamFile], stdout=f)
        tempFilesList(f.name)
        return f.name


def readRescueUpdate(df, inBam, minMapq):
    """
    Preform 'read rescue update'

    Parameters:
        df (pandas data frame): data frame from combineDf()
        inBam (str): path to input BAM file
        minMapq (int): Minimum mapq score

    Returns:
        out_sam.name (str): Path to output SAM
    """
    ok_reads = set()
    seen_starts = set()
    seen_ends = set()
    sam = inBam
    with tempfile.NamedTemporaryFile(delete=False) as out_sam:
        for line in open(sam, "r"):
            if line.startswith("@"):
                continue
            sam_fields = line.strip().split("\t")
            header = sam_fields[0]
            chrom = sam_fields[2]
            start = int(sam_fields[3])
            tlen = int(sam_fields[8])
            mapq_score = int(sam_fields[4])
            sam_flag = int(sam_fields[1])
            if (sam_flag & 256 == 256) or (sam_flag & 2048 == 2048):
                continue
            if sam_flag & 64 != 64 or tlen < 0:
                continue
            start_coord = "{}:{}".format(chrom, start)
            end_coord = "{}:{}".format(chrom, start + tlen - 1)
            if mapq_score >= minMapq:
                seen_starts.add(start_coord)
                seen_ends.add(end_coord)
        group_df_dict = {g: gdf for g, gdf in df.groupby("readName")}
        for line in open(sam, "r"):
            if line.startswith("@"):
                out_sam.write(line.encode())
                continue
            sam_fields = line.strip().split("\t")
            header = sam_fields[0]
            chrom = sam_fields[2]
            mapq_score = int(sam_fields[4])
            # Check if read has been added already. if it has write and
            if header in ok_reads:
                out_sam.write(line.encode())
                continue
            elif mapq_score >= minMapq:
                out_sam.write(line.encode())
                ok_reads.add(header)
                continue
            else:
                try:
                    sub_df = group_df_dict[header]
                except:
                    continue
                do_not_add = False
                if not sub_df.empty:
                    if not sub_df["chrom"].iloc[0] == chrom:
                        continue
                    to_add_start = []
                    to_add_end = []
                    for i, row in sub_df.iterrows():
                        coord_start = "{}:{}".format(row["chrom"], row["pos_1"])
                        if int(row["pos_2"]) != -1:
                            coord_end = "{}:{}".format(row["chrom"], row["pos_2"])
                        else:
                            coord_end = False
                        if coord_start in seen_starts:
                            do_not_add = True
                        else:
                            if coord_end:
                                if coord_end in seen_ends:
                                    do_not_add = True
                        if not do_not_add:
                            to_add_start.append(coord_start)
                            if coord_end:
                                to_add_end.append(coord_end)
                else:
                    do_not_add = True
                if not do_not_add:
                    for coord in to_add_start:
                        seen_starts.add(coord)
                    for coord in to_add_end:
                        seen_ends.add(coord)
                    ok_reads.add(header)
                    out_sam.write(line.encode())
        return out_sam.name


def countChrom(sam, hom, het):
    """
    Count the number of reads that map to the homogametic and
    heterogametic elements

    Parameters:
        sam (str): path to SAM File
        hom (str): FASTA identifier for homogametic element
        het (str): FASTA identifier for heterogametic element

    Returns:
        homCounts (int): Count of reads that mapped to homogametic element
        hetCounts (int): Count of reads that mapped to heterogametic element
    """
    hetCounts, homCounts = 0, 0
    for rec in ParseSam(sam):
        if rec.rnam == hom:
            homCounts += 1
        elif rec.rnam == het:
            hetCounts += 1
    return homCounts, hetCounts


def calculateStats(homogameticCounts, heterogameticCounts):
    """
    Calculate statistics for sequence counts

    Parameters:
        homogameticCounts (int): Number of reads that mapped to homogametic element
        heterogameticCounts (int): Number of reads that mapped to heterogametic element

    Returns:
        round(prop, 3) (float): Proportion of heterogametic to homogametic reads
                                rounded to 3 decimal places
        round(ci, 3) (float): Confidence interval for read proportions
    """
    prop = heterogameticCounts / (homogameticCounts + heterogameticCounts)
    ci = 1.96 * math.sqrt(
        (prop * (1 - prop)) / (heterogameticCounts + homogameticCounts)
    )
    return round(prop, 3), round(ci, 3)
