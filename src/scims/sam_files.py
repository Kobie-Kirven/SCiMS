# Author: Kobie Kirven, Kyle McGovern
# Davenport Lab - Penn State University
# Date: 9-2-2021


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
        while line.startswith("@") == True:
            line = self.fn.__next__()
        samData = ReadSamLine(line)
        samData.getFeatures()
        return samData


class ReadSamLine:
    """
    Read a line from a SAM file and get each of the fields
    """

    def __init__(self, samLine):
        self.samLine = samLine

    def getFeatures(self):
        if self.samLine.startswith("@") == False:
            samFields = self.samLine.strip("\n").split("\t")
            self.query = samFields[0]
            self.flag = samFields[1]
            self.rnam = samFields[2]
            self.pos = samFields[3]
            self.mapq = samFields[4]
            self.cigar = samFields[5]
            self.rnext = samFields[6]
            self.pnext = samFields[7]
            self.tlen = samFields[8]
            self.seq = samFields[9]
            self.qual = samFields[10]


def parseCigarStr(self, cigarStr):
    """Parse thought a cigar string from the SAM file
    and return a list for each cigar feature"""
    lengthOfCigarFeature = 0
    cigarList = []
    for char in cigarStr:
        if char in {"S", "I", "M", "D", "X", "N", "H", "P", "="}:
            cigarList.append((char, int(lengthOfCigarFeature)))
            lengthOfCigarFeature = 0
        else:
            try:
                lengthOfCigarFeature += int(char)
            except:
                raise Exception("Unexpected char in cigar string: {}".format(char))
    return cigarList
