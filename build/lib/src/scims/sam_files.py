# Author: Kobie Kirven, Kyle McGovern
# Davenport Lab - Penn State University
# Date: 9-2-2021
import gzip
import subprocess
import tempfile

class ParseSam:
    """Parse a SAM file and return a SAM record
    for each sequence in the SAM file

    usage example:
    for sam_record in ParseSam("sam_file.sam"):
            print(rec.query)

    """

    def __init__(self, samFile):
        self.samFile = samFile

    def __iter__(self):
        self.fn = open(self.samFile)
        return self

    def __next__(self):
        line = self.fn.__next__()
        while line.startswith("@"):
            line = self.fn.__next__()
        samData = ReadSamLine(line)
        samData.getFeatures()
        return samData


class ReadSamLine:
    """
    Read a line from a SAM file and get each of the fields.
    Naming conventions are based on the SAM format

    """

    def __init__(self, samLine):
        self.samLine = samLine

    def getFeatures(self):
        """
        Takes a line from a SAM file as
        a string collects information
        from the fields
        """

        if not self.samLine.startswith("@"):
            samFields = self.samLine.strip("\n").split("\t")
            self.query = samFields[0]
            self.flag = int(samFields[1])
            self.rnam = samFields[2]
            self.pos = int(samFields[3])
            self.mapq = int(samFields[4])
            self.cigar = samFields[5]
            self.rnext = samFields[6]
            self.pnext = samFields[7]
            self.tlen = int(samFields[8])
            self.seq = samFields[9]
            self.qual = samFields[10]
            self.align_score = None
            self.align_score_best = None
            self.align_score_mate = None
            self.num_ambig = None
            self.mismatches = None
            self.gap_opens = None
            self.gap_ext = None
            self.edit_distance = None
            self.why_filtered = None
            self.yt = None
            self.mismatch_ref_bases = None


def parseCigarStr(cigarStr):
    """Parse thought a cigar string from the SAM file
    and return a list for each cigar feature"""
    lengthOfCigarFeature = ""
    cigarList = []
    for char in cigarStr:
        if char in {"S", "I", "M", "D", "X", "N", "H", "P", "="}:
            cigarList.append((char, int(lengthOfCigarFeature)))
            lengthOfCigarFeature = ""
        else:
            try:
                lengthOfCigarFeature = lengthOfCigarFeature + char
            except:
                raise Exception("Unexpected char in cigar string: {}".format(char))
    return cigarList


def extractFromCigar(feature, cigar):
    # Exract a speficied feature from a cigar string
    return sum([match[1] for match in parseCigarStr(cigar) if match[0] == feature])


def isFileGzip(fileName):
    """
    Inputs a file name and returns a boolean 
    of if the input file is gzipped
    """
    if fileName[-3:] == ".gz":
        return True
    else:
        return False


def fastaOrFastq(fileName):
    """
    Checks if the input file is in fasta or fastq format

    args:
    file_name(str): Name of the file to be checked

    returns: A string "fasta" if the file is a fasta
            file and "fastq" if if file is fastq
    """
    if isFileGzip(fileName):
        with gzip.open(fileName, "rb") as fn:
            line = fn.readline()
            if str(line)[2:].startswith("@"):
                return "fastq"
            elif str(line)[2:].startswith(">"):
                return "fasta"
            else:
                raise IOError
    else:
        with open(fileName) as fn:
            line = fn.readline()
            if line.startswith("@"):
                return "fastq"
            elif line.startswith(">"):
                return "fasta"
            else:
                raise IOError


def decompose_sam_flag(flag):
    """
    Decompose SAM flag into its component parts

    Parameters:
        flag(int): Sam flag

    Returns:
        (list): Elements of sam flag
    """
    out_list = []
    flag_list = ["PAIRED", "PROPER_PAIR", "UNMAP", "MUNMAP",
                 "REVERSE", "MREVERSE", "READ1", "READ2",
                 "SECONDARY", "QCFAIL", "DUP", "SUPPLEMENTARY"]

    binary = str(f"{flag:b}"[::-1])
    for i in range(len(binary)):
        if binary[i] == "1":
            out_list.append(flag_list[i])
    return out_list