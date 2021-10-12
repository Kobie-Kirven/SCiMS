from Bio import SeqIO
import subprocess


def changeids(input, out):
    count = 0
    fn = open(out, "w")
    for rec in SeqIO.parse(input, "fastq"):
        fn.write("@seq" + str(count) + "\n")
        fn.write(str(rec.seq) + "\n")
        fn.write("+\n")
        quality = ""
        # print(rec.letter_annotations["phred_quality"])
        for qual in rec.letter_annotations["phred_quality"]:
            quality = quality + str(chr(qual + 33))
        fn.write(str(quality + "\n"))
        count += 1
    fn.close()


changeids("maleTest1.fastq", "male1.fastq")
changeids("maleTest2.fastq", "male2.fastq")
