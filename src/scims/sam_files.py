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
    Read a line from a SAM file and get each of the fields.
    Naming conventions are based on the SAM format

    """

    def __init__(self, samLine):
        self.samLine = samLine

    def getFeatures(self,bowtie2=False):
        '''
        Takes a line from a SAM file as
        a string collects information
        from the fields 
        '''

        if self.samLine.startswith("@") == False:
            samFields = self.samLine.strip("\n").split("\t")
            self.query = samFields[0]
            self.flag = samFields[1]
            self.rnam = samFields[2]
            self.pos = samFields[3]
            self.mapq = int(samFields[4])
            self.cigar = samFields[5]
            self.rnext = samFields[6]
            self.pnext = samFields[7]
            self.tlen = samFields[8]
            self.seq = samFields[9]
            self.qual = int(samFields[10])
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

        if bowtie2 == True:
            #If the file is a bowtie file, get the bowtie-specific
            # fields in the SAM file
            for field in samFields[10:]:
                if "AS:i" in field:
                    self.align_score = int(field)
                elif "XS:i" in field:
                    self.align_score_best = int(field)
                elif "YS:i" in field:
                    self.align_score_mate = field
                elif "XN:i" in field:
                    self.num_ambig = field
                elif "XM:i" in field:
                    self.mismatches = field
                elif "XO:i" in field:
                    self.gap_opens = field
                elif "XG:i" in field:
                    self.gap_ext = field
                elif "NM:i" in field:
                    self.edit_distance = field
                elif "YF:Z" in field:
                    self.why_filtered = field
                elif "YT:Z" in field:
                    self.yt = field
                elif "MD:Z" in field:
                    self.mismatch_ref_bases = field



def parseCigarStr(cigarStr):
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

def extractFromCigar(feature, cigar):
    #Exract a speficied feature from a cigar string
    return sum([match[1] for  in parseCigarStr(cigar) if match[0] == feature])
