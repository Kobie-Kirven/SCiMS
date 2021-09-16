# Author: Kobie Kirven, Kyle McGovern
# Davenport Lab - Penn State University
# Date: 9-2-2021

import os
from Bio import SeqIO
import gzip
from .sam_files import *
import subprocess
import pandas as pd
import numpy as np
import tempfile
import os
global filesList


filesList = []

def tempFilesList(fileName):
    '''Add the tempoary file name to the list'''
    if type(fileName) in [int, list]:
        raise TypeError
    elif fileName in filesList:
        pass
    else:
        filesList.append(fileName)

def deleteTempFileList():
    '''Delete the list of temporary files'''
    try:
        for file in filesList:
            os.unlink(file)
    except:
        pass

def alignWithBwa(index, forwardReads, reverseReads, threads):
    '''
        Align paired reads with BWA and direct output to
        pipe and return the SAM file
    '''
    with tempfile.NamedTemporaryFile(delete=False) as f:
        subprocess.run(["bwa","mem","-t", threads, index, 
            forwardReads, reverseReads], stdout=f)
        tempFilesList(f.name)
        return f.name


def getHumanSequences(inputSam):
    '''
    Get the sequences from BWA alignment output
    that map to the human genome
    '''
    with tempfile.NamedTemporaryFile(delete=False) as f:
        subprocess.run(["samtools", "view", "-F", "4",
            "-F", "8", "-b", inputSam], stdout=f)
        tempFilesList(f.name)
    return f.name


def getSexSequences(inputBam, homogamete, heterogamete, mapqMin=50, tlen=75):
    '''
    Input two strings for the homogametic element and one for the
    heterogametic element
    '''
    with tempfile.NamedTemporaryFile(delete=False) as f:
        smallSam = subprocess.run(["samtools","view","-F", "2048", "-F", "256", inputBam],
            stdout=f)
        name = f.name
        tempFilesList(f.name)

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

def getUniqueIdsInSam(inputSam):
    '''Get the unique query names from the SAM file'''
    readIds = []
    with open(inputSam) as fn:
        for rec in fn:
            line = ReadSamLine(rec)
            line.getFeatures()
            readIds.append(line.query)
    return set(readIds)


def getFastqReadsInSam(readIds, forwardReads, reverseReads):

    # Open output files
    outFastq1, outFastq2 = tempfile.NamedTemporaryFile(delete=False), tempfile.NamedTemporaryFile(delete=False)

    # Get whether the file is fasta or fastq
    fileType = fastaOrFastq(forwardReads)

    # Open the input read files and get those FASTQ sequences from the SAM file
    if isFileGzip(forwardReads) == True:
        forwardReads = gzip.open(forwardReads, "rt")
        reverseReads = gzip.open(reverseReads, "rt")

    for rec1 in SeqIO.parse(forwardReads, fileType):
        if str(rec1.id)[:-2] in readIds:
            outFastq1.write(("@" + str(rec1.description) + "\n").encode())
            outFastq1.write((str(rec1.seq) + "\n").encode())
            outFastq1.write("+\n".encode())
            quality = ""
            for qual in rec1.letter_annotations["phred_quality"]:
                quality = quality + str(chr(qual + 33))
            outFastq1.write((str(quality) + "\n").encode())

    for rec2 in SeqIO.parse(reverseReads, fileType):
        if str(rec2.id)[:-2] in readIds:
            outFastq2.write(("@" + str(rec2.description) + "\n").encode())
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
    # Merge paired-end reads with flash
    tempdir = tempfile.TemporaryDirectory()
    subprocess.run(["flash", "-x",str(mismatchRatio), "-M", str(maxOverlap),"-m", str(minOverlap),
        pair1, pair2, "-d", tempdir.name])
    merged = convertToTemp(tempdir.name + "out.extendedFrags.fastq")
    unmerged1 = convertToTemp(tempdir.name + "out.notCombined_1.fastq")
    unmerged2 = convertToTemp(tempdir.name + "out.notCombined_2.fastq")
    return merged, unmerged1, unmerged2

def convertToTemp(oldFile):
    with open(oldFile) as fn:
        with tempfile.NamedTemporaryFile(delete=False) as f:
            for line in fn:
                f.wirte(line.encode)
            return f.name

def alignMergedWithBowtie2(index, merged, threads=1, noOfAlignPerRead=50):
    """Align merged reads with bowtie2"""
    with tempfile.NamedTemporaryFile(delete=False) as f:
        subprocess.run(["bowtie2","-p",str(threads),"-k", str(noOfAlignPerRead),"--local",
            "-x", index, "-U", merged], stdout=f)
        tempFilesList(f.name)
        return f.name

def pairedReadsWithBowtie2(index, forwardReads, reverseReads, threads=1, noOfAlignPerRead=50):
    """Align paired reads with bowtie2"""
    with tempfile.NamedTemporaryFile(delete=False) as f:
        subprocess.run(["bowtie2","-p",str(threads),"-k", str(noOfAlignPerRead),"--local",
            "-x", index, "-1", forwardReads, "-2", reverseReads], stdout=f)
        return f.name


def getReadDataFromMergedSam(samFile):
    '''Get data from the SAM file of merged reads'''
    fastqData = []
    for rec in ParseSam(samFile):
        if rec.cigar != "*":
            numMatches = extractFromCigar("M", rec.cigar)
            numDeletions = extractFromCigar("D", rec.cigar)

            if rec.flag != 4 and rec.flag != 8:
                readData = {"readName":rec.query, "chrom":rec.rnam,
                "alignmentScore":rec.align_score, "numMismatches":rec.mismatches,
                "insertSize": None, "pos_1": rec.pos, "pos_2": -1,
                "numMatches":numMatches, "numDeletions":numDeletions}

                fastqData.append(readData)
    return fastqData

def getReadDataFromUnmergedSam(unmergedSam):
    fastqData = []
    for rec in ParseSam(samFile):
        if rec.cigar != "*":
            numMatches = extractFromCigar("M", rec.cigar)
            numDeletions = extractFromCigar("D", rec.cigar)
            readData = {"readName":rec.query, "samFlag":rec.flag, "chrom":rec.rnam,
            "alignmentScore":rec.align_score, "numMismatches":rec.mismatches,
            "mateStart":rec.pnext, "insertSize":rec.tlen, "pos":rec.pos,
            "numMatches":numMatches, "numDeletions":numDeletions}
            fastqData.append(readDdata)
    return fastqData


def readRescueMerged(mergedSam, unmergedSam, diffFromMaxPercent=0.025, 
    minMatch=74, mismatchRatio=0.025, numDeletions=3):

    fastqData = getReadDataFromMergedSam(mergedSam)

    if fastqData:
        #Transform the fastqData into a pandas dataframe
        df = pd.DataFrame(fastqData)
        df_unfilt = df.copy(deep=True)

        #Filter reads that arent within a certian quality percentage
        # from the maximum scoring alignment for each read
        df["maxScoreReadName"] = df.groupby("readName")["alignmentScore"].transform(max)
        df["diffPercentMaxRound"] = np.round(df["maxScoreReadName"] * diffFromMaxPercent)
        df["scoreMaxDiff"] = df["maxScoreReadName"] - df["alignmentScore"]
        df = df[df["scoreMaxDiff"] <= df["diffPercentMaxRound"]]

        # Filter out reads with less than 75 bp matching.
        df = df[df["numMatches"] > minMatch]
        # Filter out reads with too many mismatches (more than 1 mismatch every 40 bp).
        df["mismatchRatio"] = df["numMismatches"].astype(float) / df["numMatches"].astype(float)
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


def readRescueUnmerged(mergedSam, diffFromMaxPercent=0.025, minMatch=100, numDeletions=3):

    read_pair = []
    fastq_data = []
    fastqData = getReadDataFromUnmergedSam(mergedSam)

    for readData in fastqData:
        read_pair.append(readData)
        if not read_pair:
            if readData["samFlag"] & 2 == 2 and readData["samFlag"] & 64 == 64:
                read_pair.append(readData)

        else:
            if (abs(readData["insertSize"]) == abs(read_pair[0]["insertSize"]) and readData["chrom"] == read_pair[0]["chrom"] and readData["readName"] == read_pair[0]["readName"]):
                read_pair.append(readData)
                total_num_matches = read_pair[0]["numMatches"] + read_pair[1]["numMatches"]
                total_num_mismatches = read_pair[0]["numMismatches"] + read_pair[1]["numMismatches"]
                total_num_deletions = read_pair[0]["numDeletions"] + read_pair[1]["numDeletions"]
                total_alignment_score = read_pair[0]["alignmentScore"] + read_pair[1]["alignmentScore"]
                fastq_data.append({"readName": readData["readName"], "numMismatches": total_num_mismatches, "chrom": readData["chrom"],
                                   "numMatches": total_num_matches, "numDeletions": total_num_deletions,
                                   "insertSize": abs(readData["insertSize"]), "alignmentScore": total_alignment_score,
                                   "pos_1": read_pair[0]["pos"], "pos_2": read_pair[1]["pos"]})
                read_pair = []

    df_pair = pd.DataFrame(fastq_data)
    df_pair_unfilt = df_pair.copy(deep=True)

    # Filter reads with big insert size (gt 999)
    df_pair = df_pair[df_pair["insertSize"] < 1000]

    # Filter out reads less than 2.5% max score.
    df_pair["maxScoreReadName"] = df_pair.groupby("readName")["alignmentScore"].transform(max)
    df_pair["diffPercentMaxRound"] = np.round(df_pair["maxScoreReadName"] * diffFromMaxPercent)
    df_pair["scoreMaxDiff"] = df_pair["maxScoreReadName"] - df_pair["alignmentScore"]
    df_pair = df_pair[df_pair["scoreMaxDiff"] <= df_pair["diffPercentMaxRound"]]

    # Filter out reads with less than 100 bp matching.
    df_pair = df_pair[df_pair["numMatches"] >= minMatch]

    # Filter out reads with too many mismatches (more than 1 mismatch every 40 bp).
    df_pair["mismatchRatio"] = df_pair["numMismatches"].astype(float) / df_pair["numMatches"].astype(float)
    df_pair = df_pair[df_pair["mismatchRatio"] < 0.025]

    # Filter out reads with too many deletions.
    df_pair = df_pair[df_pair["numDeletions"] <= numDeletions]

    # Finally keep reads that have only one annotated chromosome.
    mapper = df_pair.groupby("readName")["chrom"].nunique().to_dict()
    df_pair["uniqueChroms"] = df_pair["readName"].map(mapper)
    df_pair = df_pair[df_pair["uniqueChroms"] == 1]

    


