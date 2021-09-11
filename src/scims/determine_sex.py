# Author: Kobie Kirven, Kyle McGovern
# Davenport Lab - Penn State University
# Date: 9-2-2021

import os
from Bio import SeqIO
import gzip
from .build_index import TestFile
from .sam_files import *
import subprocess
import pandas as pd
import numpy as np
import tempfile
import os
global filesList


filesList = []

def tempFilesList(fileName):
    filesList.append(fileName)

def deleteTempFileList():
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
    return tempdir

def alignMergedWithBowtie2(index, merged, threads=1, noOfAlignPerRead=50):
"""Align merged reads with bowtie2"""
with tempfile.NamedTemporaryFile(delete=False) as f:
    subprocess.run(["bowtie2","-p",threads,"-k", noOfAlignPerRead,"--local",
        "-x", index, "-U", merged], stdout=f)
    tempFilesList(f.name)
    return f.name

def pairedReadsWithBowtie2(index, forwardReads, reverseReads, threads=1, noOfAlignPerRead=50):
    """Align paired reads with bowtie2"""
    with tempfile.NamedTemporaryFile(delete=False) as f:
        subprocess.run(["bowtie2","-p",threads,"-k", noOfAlignPerRead,"--local",
            "-x", index, "-1", forwardReads, "-2", reverseReads], stdout=f)
        return f.name

def bowtie2ReadRescue(mergedSam, unmergedSam, diffFromMaxPercent=0.025):

    fastqData = {}
    for rec in ParseSam(mergedSam):
        numMatches = extractFromCigar("M", rec.cigar)
        numDeletions = extractFromCigar("D", rec.cigar)

        readData = {"readName":rec.query, "chrom":rec.rnam,
        "alignmentScore":query.align_score, "mismatches":rec.mismatches,
        "insert_size": None, "pos_1": rec.pos, "pos_2": -1,
        "numMatches":numMatches, "numDeletions":numDeletions}

        fastqData.append(read_data)

    #Transform the fastqData into a pandas dataframe
    df = pd.DataFrame(fastq_data)
    df_unfilt = df.copy(deep=True)

    #Filter reads that arent within a certian quality percentage
    # from the maximum scoring alignment for each read
    df["maxScoreReadName"] = df.groupby("readName")["alignmentScore"].transform(max)
    df["2point5_percent_max_round"] = np.round(df["read_name_score_max"] * diffFromMaxPercent)
    df["score_max_diff"] = df["read_name_score_max"] - df["alignment_score"]
    df = df[df["score_max_diff"] <= df["2point5_percent_max_round"]]



