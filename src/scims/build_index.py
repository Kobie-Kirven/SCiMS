# Author: Kobie Kirven
# Davenport Lab - Penn State University
# Date: 9-2-2021

# Imports
import subprocess


def buildBwaIndex(inputFile, outputName):
    """
    Build BWA index from reference genome

    Parameters:
        inputFile (str): path to reference genome FASTA file
        outputName (str): prefix for the BWA index
    """
    subprocess.run(["bwa index " +  inputFile + " " + outputName], shell=True)

def buildBowtie2Index(inputFile, outputName):
    """
    Build index for Bowtie2 from the reference genome FASTA file

    Parameters:
         inputFile (str): path to reference genome FASTA file
        outputName (str): prefix for the Bowtie2 index

    """
    subprocess.run(["bowtie2-build " + inputFile + " " + outputName], shell=True)
