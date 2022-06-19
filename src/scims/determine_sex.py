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

    Returns:
        f.name (str): Path to output SAM file
    """
    with tempfile.NamedTemporaryFile(delete=False) as f:
        subprocess.run(
            [
                "bowtie2",
                "-p" +
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

def single_reads_with_bowtie2(
        reads, index, method, threads=1
):
    """
    Align paired reads with bowtie2

    Parameters:
        reads (str): Path to read file for alignment
        index (str): Prefix for Bowite2 index
        threads (int): Number of threads to use
    Returns:
        f.name (str): Path to output SAM file
    """
    with tempfile.NamedTemporaryFile(delete=False) as f:
        subprocess.run(
            [
                "bowtie2",
                "-p" +
                str(threads),
                "--" + method,
                "-x",
                index,
                "-U",
                reads
            ],
            stdout=f,
            stderr=subprocess.DEVNULL,
        )
        temp_files_list(f.name)
        return f.name

def not_in_list(not_list, input_list):
    """
    Check to make sure that no elements in one list are in another list

    Parameters:
        not_list (list): List of elements that should not be in input_list
        input_list (list): List of elements that are not wanted

    Returns:
        (bool): Returns True if no elements in not_list are in input_list
    """
    flag = True
    for ele in not_list:
        if ele in input_list:
            return False
    return True

def read_scaffold_names(scaffolds_file):
    scaffold_list = []
    with open(scaffolds_file) as fn:
        lines = fn.readlines()
        for line in lines:
            if line.startswith(">"):
                scaffold_list.append(line[1:].strip("\n"))
            else:
                scaffold_list.append(line.strip("\n"))
    return scaffold_list

def read_scaffold_lengths(scaffolds_file_lengths):
    """
    Read the lengths of the scaffolds from a file

    Parameters:
        scaffolds_file_lengths (str): Path to file containing scaffold lengths

    Returns:
        chrom_lengths_dict (dict): Dictionary of scaffold lengths
    """
    chrom_lengths_dict = {}
    with open(scaffolds_file_lengths) as fn:
        lines = fn.readlines()
        for line in lines:
            line = line.split("\t")
            chrom_lengths_dict[line[0]] = int(line[1])
    return chrom_lengths_dict
    
def determine_chrom_lengths(reference_genome, scaffold_list=None):
    """
    Compute the lengths of the chromosome scaffolds

    Parameters:
        reference_genome (str): path to the reference genome

    Returns:
        chrom_lengths(dict): Scaffold names are keys and values are scaffold lengths
    """
    chrom_lengths = {}

    # If the reference genome is gzipped, open it as a gzip file
    if reference_genome[-3:]==".gz":
        with gzip.open(reference_genome, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if scaffold_list is None:
                    count = 0
                    for let in str(record.seq):
                        if let != "N":
                            count += 1
                    chrom_lengths[record.id] = count

                elif record.id in scaffold_list:
                    count = 0
                    for let in str(record.seq):
                        if let != "N":
                            count += 1
                    chrom_lengths[record.id] = count
    else:
        for record in SeqIO.parse(reference_genome, "fasta"):
            if scaffold_list is None:
                count = 0
                for let in str(record.seq):
                    if let != "N":
                        count += 1
                chrom_lengths[record.id] = count
            elif record.id in scaffold_list:
                count = 0
                for let in str(record.seq):
                    if let != "N":
                        count += 1
                chrom_lengths[record.id] = count
    return chrom_lengths

def count_chrom_alignments(sam_file, scaffold_list=[]):
    """
    Count all alignments that match mapping criteria

    Parameters:
        sam_file (str): path to SAM file to be used for analysis
        scaffold_list(list): List of scaffold names to be used

    Returns:
        hit_counts (dict): Keys are scaffold names and number of counts are values
    """
    hit_counts = {}
    for rec in ParseSam(sam_file):
        if not_in_list(["SECONDARY", "UNMAP", "MUNMAP", "SUPPLEMENTARY"], decompose_sam_flag(rec.flag)) == True:
            if scaffold_list:
                if rec.rnam in scaffold_list:
                    if rec.rnam not in hit_counts:
                        hit_counts[rec.rnam] = 1
                    else:
                        hit_counts[rec.rnam] += 1
            else:
                if rec.rnam not in hit_counts:
                    hit_counts[rec.rnam] = 1
                else:
                    hit_counts[rec.rnam] += 1
    return hit_counts

def normalize_by_chrom_lengths(counts_dict, chrom_lengths_dict):
    """
    Normalize the number of counts by the length of the chromosome

    Parameters:
        counts_dict(dict): count_chrom_alignments
        chrom_lengths_dict(dict): output from determine_chrom_lengths()

    Returns:
        counts_dict (dict):
    """
    for ele in counts_dict:
        counts_dict[ele] = counts_dict[ele] / float(chrom_lengths_dict[ele])
    return counts_dict

def compare_to_homogametic(counts_dict_normal, homogametic_element):
    """
    Compare the coverages to the coverage of the homogametic element

    Parameters:
        counts_dict_normal(dict): output from normalize_by_chrom_lengths()
        chrom_lengths_dict(dict): output from determine_chrom_lengths()

    Returns:
        counts_dict(dict): keys are chromosome names and values are normalized coverages
    """
    counts_dict = {}
    for ele in counts_dict_normal:
        counts_dict[ele] =  counts_dict_normal[homogametic_element] / counts_dict_normal[ele]
    return counts_dict

def count_seqs(count_dict):
    """
    Counts the total number of valid alignments

    Parameters:
        count_dict(dict): output from count_chrom_alignments()

    Returns:
        sum (int): Total number of valid alignments
    """
    sum = 0
    for id in count_dict:
        sum += int(count_dict[id])
    return sum


def generate_plot(counts_dict_plot, output, homogametic_element, heterogametic_element):
    """
    Create a plot from chromosomal coverage data
    """
    ids = []
    for ele in counts_dict_plot:
        if ele != heterogametic_element:
            ids.append(ele)
    ids.sort()
    counts = []
    for ele in ids:
        if ele != heterogametic_element:
            counts.append(counts_dict_plot[ele])
    plt.scatter(ids, counts, color="black")
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
    return sum(counts_list) / len(counts_list)