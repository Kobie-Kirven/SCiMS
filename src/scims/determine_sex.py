# Author: Kobie Kirven
# Davenport Lab - Penn State University
# Date: 9-2-2021

import os
from Bio import SeqIO
import gzip
from .build_index import TestFile
import subprocess


class GetProportion:
    #

    def __init__(self, index, forward, reverse, threads):
        self.index = index
        self.forward = forward
        self.reverse = reverse
        self.threads = threads

    def getHumanSequences(self):
        # Get the sequences that map to the human genome
        print("Now Aligning Reads")
        subprocess.call(
            "bwa mem -t "
            + self.threads
            + " "
            + self.index
            + " "
            + self.forward
            + " "
            + self.reverse
            + " | samtools view -F 4 -F 8 -b -o temp.bam",
            shell=True,
        )

    def getSexSequences(self, homogamete, heterogamete, outFile):
        # Find
        subprocess.call(
            "samtools view -F 2048 -F 256 temp.bam"
            + " | awk '$5 < 50 && ($3 =="
            + '"'
            + homogamete
            + '"'
            + "|| $3 =="
            + '"'
            + heterogamete
            + '"'
            + ") && ($9 > 75 || $9 < -75) {print}' > "
            + outFile,
            shell=True,
        )

    def getFastqReadsInSam(self):
        # Get the unique sequence names from the SAM file
        with open("reads.sam") as fn:
            lines = fn.readlines()
            readIds = set([line.split("\t")[0] for line in lines])
            # print(readIds)

        # Open output files
        outFastq1, outFastq2 = open("align1.fastq", "w"), open("align2.fastq", "w")

        # Get whether the file is fasta or fastq
        fileType = TestFile.fastaOrFastq(self.forward)

        # Open the input read files and get those FASTA sequences from the SAM file
        if TestFile.isFileZip(self.forward) == True:
            forward = gzip.open(self.forward, "rb")
            reverse = gzip.open(self.reverse, "rb")

        else:
            for rec1 in SeqIO.parse(self.forward, fileType):
                if str(rec1.id) in readIds:
                    SeqIO.write(rec1, outFastq1, fileType)

            for rec2 in SeqIO.parse(self.reverse, fileType):
                if str(rec2.id) in readIds:
                    SeqIO.write(rec2, outFastq2, fileType)

        outFastq1.close()
        outFastq2.close()

    def mergeWithFlash(self):
        # Merge paired-end reads with flash
        subprocess.call(
            "flash -x .1 -M 150 -m 40 " + "align1.fastq align2.fastq -o flash",
            shell=True,
        )

    def alignWithBowtie2(self, merged, input1, input2=None):
        # Align reads with bowtie2
        if merged == True:
            subprocess.call(
                "bowtie2 -p "
                + self.threads
                + " -k 50 --local -x "
                + self.index
                + " -U "
                + input1
                + " > bowtie2_merged.sam",
                shell=True,
            )

        elif merged == False:
            subprocess.call(
                "bowtie2 -p "
                + self.threads
                + " -k 50 --local -x "
                + self.index
                + " -1 "
                + input1
                + " -2 "
                + input2
                + " > bowtie2_unmerged.sam",
                shell=True,
            )

    # def bowtie2ReadRescue():
