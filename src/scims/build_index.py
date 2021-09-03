# Author: Kobie Kirven
# Davenport Lab - Penn State University
# Date: 9-2-2021

from Bio import SeqIO
import gzip
import os


class BuildIndex:
    """Class for building bwa and bowtie2 indexes"""

    def __init__(self, reference, homogamete, heterogamete, output):
        self.reference = reference
        self.homo = homogamete
        self.hetero = heterogamete
        self.output = output

    def isFileZip(self):
        # Check if the file is gzipped
        if self.reference[-3:] == ".gz":
            return True
        else:
            return False

    def testFileType(self):
        # Check if the input file is fasta or fastq
        if BuildIndex.isFileZip(self) == True:
            with gzip.open(self.reference, "rb") as fn:
                line = fn.readline()
                if str(line)[2:].startswith("@"):
                    return "fastq"
                elif str(line)[2:].startswith(">"):
                    return "fasta"
                else:
                    raise IOError
        else:
            with open(self.reference) as fn:
                line = fn.readline()
                if line.startswith("@"):
                    return "fastq"
                elif line.startswith(">"):
                    return "fasta"
                else:
                    raise IOError

    def getSexChromosomes(self):
        # Extract the sex chromosomes from the reference genome
        # and output them into a fasta file
        outputFile = open(self.output + ".fasta", "w")

        if BuildIndex.isFileZip(self) == True:
            with gzip.open(self.reference, "rt") as handle:
                for rec in SeqIO.parse(handle, BuildIndex.testFileType(self)):
                    if rec.id == self.homo or rec.id == self.hetero:
                        SeqIO.write(rec, outputFile, "fasta")
        else:
            for rec in SeqIO.parse(self.reference, BuildIndex.testFileType(self)):
                if rec.id == self.homo or rec.id == self.hetero:
                    SeqIO.write(rec, outputFile, "fasta")
        return self.output + ".fasta"

    def buildBwaIndex(self):
        # Build the index for bwa
        os.system("bwa index " + BuildIndex.getSexChromosomes(self))

    def buildBowtie2Index(self):
        # Build the index for Bowtie2
        os.system(
            "bowtie2-build "
            + BuildIndex.getSexChromosomes(self)
            + " "
            + BuildIndex.getSexChromosomes(self)
        )
