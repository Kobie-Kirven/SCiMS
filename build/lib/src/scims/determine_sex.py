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
import gzip
import matplotlib.pyplot as plt
import scipy.stats as stats

files_list = []


def temp_files_list(file_name):
    """
    Adds the temporary file to the list of temporary files
    
    Args: file_name(str): path to the temporary file
    """
    if type(file_name) in [int, list]:
        raise TypeError
    elif file_name in files_list:
        pass
    else:
        files_list.append(file_name)


def delete_temp_file_list():
    """
    Delete the list of temporary files
    """
    try:
        for file in files_list:
            os.unlink(file)
    finally:
        pass


def check_index_files(endings, index):
    """
    Check whether index files are present

    Parameters:
        endings (list): list of endings for
        index (str): Prefix of the index
    """
    with tempfile.NamedTemporaryFile(delete=False) as f:
        if "/" in index:
            path = "/".join(index.split("/")[:-1])
            index = index.split("/")[-1]
            subprocess.run(["ls", path], stdout=f)
        else:
            subprocess.run(["ls"], stdout=f)
    bool_count = 0
    with open(f.name) as fn:
        lines = fn.readlines()
        for pattern in endings:
            for i in range(len(lines)):
                if lines[i].strip("\n").endswith(pattern) and index in lines[i]:
                    bool_count += 1
    os.unlink(f.name)
    if bool_count == len(endings):
        return True
    else:
        return False

def verify_sam_file(sam_file):
    """
    Verify that the SAM file is valid

    Parameters:
        sam_file (str): Path to SAM file

    Returns:
        (bool): Whether the SAM file is valid or not
    """
    with open(sam_file) as fn:
        flag = False
        for line in fn:
            if line.startswith("@"):
                continue
            else:
                flag = True
                if len(line.split("\t")) >= 11:
                    return True
                else:
                    return False
        return flag

def get_human_sequences(input_sam):
    """
    Get the sequences from BWA alignment output
    that map to the human genome

    Parameters:
        input_sam (str): Path to input SAM file

    Returns:
        f.name (str): Path to output BAM file
    """
    with tempfile.NamedTemporaryFile(delete=False) as f:
        if verify_sam_file(input_sam):
            subprocess.run(
                ["samtools", "view", "-F", "4", "-F", "8", "-b", input_sam],
                stdout=f,
                stderr=subprocess.DEVNULL,
            )
            temp_files_list(f.name)
            return f.name
        else:
            raise Exception("Not a valid SAM file")


def get_unique_ids_in_sam(input_sam):
    """
    Get the unique query names from the SAM file

    Parameters:
        input_sam (str): Path to input SAM file

    returns:
        (list): Set of unique FASTA IDs
    """
    read_ids = []
    with open(input_sam) as fn:
        for rec in fn:
            line = ReadSamLine(rec)
            line.getFeatures()
            read_ids.append(line.query)
    return set(read_ids)


def get_fastq_reads_in_sam(read_ids, forward_reads, reverse_reads):
    """
    Get FASTQ reads for sequences in SAM file

    Parameters:
        read_ids (list): List of FASTA IDs
        forward_reads (str): Path to forward FASTQ reads
        reverse_reads (str): Path to reverse FASTQ reads

    Returns:
        out_fastq1.name (str): Path to output forward FASTQ file
        out_fastq2.name (str): Path to output reverse FASTQ file
    """

    # Open output files
    out_fastq1, out_fastq2 = (
        tempfile.NamedTemporaryFile(delete=False),
        tempfile.NamedTemporaryFile(delete=False),
    )

    # Get whether the file is FASTA or FASTQ
    file_type = fastaOrFastq(forward_reads)

    # Open the input read files and get those FASTQ sequences from the SAM file
    if isFileGzip(forward_reads):
        forward_reads = gzip.open(forward_reads, "rt")
        reverse_reads = gzip.open(reverse_reads, "rt")

    for rec1 in SeqIO.parse(forward_reads, file_type):
        if str(rec1.id)[:-2] in read_ids:
            out_fastq1.write(("@" + str(rec1.description)[:-2] + "\n").encode())
            out_fastq1.write((str(rec1.seq) + "\n").encode())
            out_fastq1.write("+\n".encode())
            quality = ""
            for qual in rec1.letter_annotations["phred_quality"]:
                quality = quality + str(chr(qual + 33))
            out_fastq1.write((str(quality) + "\n").encode())

    for rec2 in SeqIO.parse(reverse_reads, file_type):
        if str(rec2.id)[:-2] in read_ids:
            out_fastq2.write(("@" + str(rec2.description)[:-2] + "\n").encode())
            out_fastq2.write((str(rec2.seq) + "\n").encode())
            out_fastq2.write("+\n".encode())
            quality = ""
            for qual in rec2.letter_annotations["phred_quality"]:
                quality = quality + str(chr(qual + 33))
            out_fastq2.write((str(quality) + "\n").encode())

    out_fastq1.close()
    out_fastq2.close()
    temp_files_list(out_fastq1.name)
    temp_files_list(out_fastq2.name)
    return out_fastq1.name, out_fastq2.name


def paired_reads_with_bowtie2(
        index, forward_reads, reverse_reads, method, threads=1
):
    """
    Align paired reads with bowtie2

    Parameters:
        index (str): Prefix for Bowite2 index
        forward_reads (str): Path to forward FASTQ reads
        reverse_reads (str): Path to reverse FASTq reads
        threads (int): Number of threads to use
        no_of_align_per_read (int): Maximum number of alignments per-read

    Returns:
        f.name (str): Path to output SAM file
    """
    with tempfile.NamedTemporaryFile(delete=False) as f:
        subprocess.run(
            [
                "bowtie2",
                "-p",
                str(threads),
                "--" + method,
                "-x",
                index,
                "-1",
                forward_reads,
                "-2",
                reverse_reads,
            ],
            stdout=f,
            stderr=subprocess.DEVNULL,
        )
        temp_files_list(f.name)
        return f.name

def not_in_list(not_list, input_list):
    flag = True
    for ele in not_list:
        if ele in input_list:
            return False
    return True

def determine_chrom_lengths(reference_genome):
    chrom_lengths = {}
    for record in SeqIO.parse(reference_genome, "fasta"):
        if "NW" not in record.id and "NT" not in record.id:
            chrom_lengths[record.id] = len(str(record.seq))
    return chrom_lengths


def count_chrom_alignments(sam_file):
    hit_counts = {}
    for rec in ParseSam(sam_file):
        if not_in_list(["SECONDARY", "UNMAP", "MUNMAP", "SUPPLEMENTARY"], decompose_sam_flag(rec.flag)) == True:
            if "NW" not in rec.rnam and "NT" not in rec.rnam:
                if rec.rnam not in hit_counts:
                    hit_counts[rec.rnam] = 1
                else:
                    hit_counts[rec.rnam] += 1
    return hit_counts

def normalize_by_chrom_lengths(counts_dict, chrom_lengths_dict):
    for ele in counts_dict:
        counts_dict[ele] = counts_dict[ele] / float(chrom_lengths_dict[ele])
    return counts_dict

def compare_to_homogametic(counts_dict_normal, homogametic_element):
    counts_dict = {}
    for ele in counts_dict_normal:
        counts_dict[ele] =  counts_dict_normal[homogametic_element] / counts_dict_normal[ele]
    return counts_dict

def generate_plot(counts_dict_plot, output, homogametic_element, heterogametic_element):
    ids, counts = [], []
    for ele in counts_dict_plot:
        if ele != heterogametic_element:
            ids.append(ele)
            counts.append(counts_dict_plot[ele])
    plt.scatter(ids, counts)
    plt.ylim([0, 2])
    plt.title(homogametic_element + ":Autosomal Coverage Ratio")
    plt.ylabel('Normalized Coverage')
    plt.xlabel('Chromosome ID')
    plt.xticks(rotation=-90)
    plt.savefig(output + '.png', bbox_inches='tight')

def stats_test(counts_dict, homogametic_element, heterogametic_element):
    counts_list = []
    for ele in counts_dict:
        if ele != homogametic_element and ele != heterogametic_element:
            counts_list.append(counts_dict[ele])
    male = stats.ttest_1samp(counts_list, 0.5, alternative='two-sided')
    female = stats.ttest_1samp(counts_list, 1, alternative='two-sided')
    return male, female