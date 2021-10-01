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


def check_bwa_index(index):
    """
    Verify index files are present for BWA

    Parameters:
        index (str): Prefix for BWA index
    """
    if check_index_files(["amb", "ann", "bwt", "pac", "sa"], index):
        return True
    else:
        raise Exception("BWA index files can not be found")


def align_with_bwa(index, forward_reads, reverse_reads, threads):
    """
        Align paired-end reads with BWA

        Parameters:
            index (str): prefix for BWA index
            forward_reads (str): path to file with forward FASTQ reads
            reverse_reads (str): path to file with reverse FASTQ reads
            threads (int): number of threads to use

        Returns:
            f.name (str): Path to output SAM file
    """
    check_bwa_index(index)
    with tempfile.NamedTemporaryFile(delete=False) as f:
        subprocess.run(["bwa", "mem", "-t", threads, index, forward_reads, reverse_reads],
                       stdout=f, stderr=subprocess.DEVNULL)
        temp_files_list(f.name)
        return f.name

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


def get_sex_sequences(
        input_bam, homogamete, heterogamete, is_filter=True, mapq_min=50, tlen=75
):
    """
    Get the sequences that map to the sex chromosomes

    Parameters:
        input_bam (str): Path to input BAM file
        homogamete (str): FASTA ID for homogametic element
        heterogamete (str): FASTA ID for heterogametic element
        is_filter (bool): Filter reads based on quality thresholds
        mapq_min (int): Minimum mapq score
        tlen (int): Template length

    Returns:
        g.name (str): Path to SAM file with reads that mapped to sex chromosomes

    """
    with tempfile.NamedTemporaryFile(delete=False) as f:
        subprocess.run(["samtools", "view", "-F", "2048", "-F", "256", input_bam], stdout=f)
        temp_files_list(f.name)

    if is_filter:
        with tempfile.NamedTemporaryFile(delete=False) as g:
            for line in open(f.name):
                rec = ReadSamLine(line)
                rec.getFeatures()
                if rec.mapq < mapq_min:
                    if rec.rnam == homogamete or rec.rnam == heterogamete:
                        if rec.tlen > tlen or rec.tlen < (tlen * -1):
                            g.write(line.encode())
            temp_files_list(g.name)
            return g.name
    else:
        with tempfile.NamedTemporaryFile(delete=False) as g:
            for line in open(f.name):
                rec = ReadSamLine(line)
                rec.getFeatures()
                if rec.rnam == homogamete or rec.rnam == heterogamete:
                    if rec.tlen > tlen or rec.tlen < (tlen * -1):
                        g.write(line.encode())
            temp_files_list(g.name)
            return g.name


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


def merge_with_flash(pair1, pair2, mismatch_ratio=0.1, max_overlap=150, min_overlap=40):
    """
    Merge paired-end reads with FLASH

    Parameters:
        pair1 (str): Path to forward FASTQ file
        pair2 (str): Path to reverse FASTQ file
        mismatch_ratio (float): Maximum mismatch ratio
        max_overlap (int): Maximum overlap for read merger
        min_overlap (int): Minimum overlap for read merger

    Returns:
        merged (str): Path to merged FASTQ file
        unmerged1 (str): Path to forward, unmerged FASTQ file
        unmerged2 (str): Path to reverse, unmerged FASTQ file
    """
    tempdir = tempfile.TemporaryDirectory()
    subprocess.run(
        [
            "flash",
            "-x",
            str(mismatch_ratio),
            "-M",
            str(max_overlap),
            "-m",
            str(min_overlap),
            pair1,
            pair2,
            "-d",
            tempdir.name,
        ],
        stdout=subprocess.DEVNULL,
    )
    merged = convert_to_temp(tempdir.name + "/out.extendedFrags.fastq")
    unmerged1 = convert_to_temp(tempdir.name + "/out.notCombined_1.fastq")
    unmerged2 = convert_to_temp(tempdir.name + "/out.notCombined_2.fastq")
    return merged, unmerged1, unmerged2


def convert_to_temp(old_file):
    """
    Convert a file into a temporary file

    Parameters:
        old_file (str): Path to input file

    Returns:
        f.name (str): Path to output temporary file
    """
    with open(old_file) as fn:
        with tempfile.NamedTemporaryFile(delete=False) as f:
            for line in fn:
                f.write(line.encode())
            return f.name


def align_merged_with_bowtie2(index, merged, threads=1, no_of_align_per_read=50):
    """
    Align merged reads with bowtie2

    Parameters:
        index (str): Prefix for Bowite2 index
        merged (str): Path to FLASH merged, FASTQ reads
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
                "-k",
                str(no_of_align_per_read),
                "--local",
                "-x",
                index,
                "-U",
                merged,
            ],
            stdout=f,
            stderr=subprocess.DEVNULL,
        )
        temp_files_list(f.name)
        return f.name


def paired_reads_with_bowtie2(
        index, forward_reads, reverse_reads, threads=1, no_of_align_per_read=50
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
                "-k",
                str(no_of_align_per_read),
                "--local",
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
        return f.name


def get_read_data_from_merged_sam(sam_file):
    """
    Get data from the SAM file of merged reads

    Parameters:
        sam_file (str): Path to input SAM file

    Returns:
        fastq_data (list): A list of dictionaries where each element
                            corresponds to each line in the SAM file
    """
    fastq_data = []
    for rec in ParseSam(sam_file):
        if rec.cigar != "*":
            num_matches = extractFromCigar("M", rec.cigar)
            num_deletions = extractFromCigar("D", rec.cigar)

            if "UNMAP" not in decompose_sam_flag(rec.flag) and "MUNMAP" not in decompose_sam_flag(rec.flag):
                read_data = {
                    "read_name": rec.query,
                    "chrom": rec.rnam,
                    "alignment_score": rec.align_score,
                    "num_mismatches": rec.mismatches,
                    "insert_size": None,
                    "pos_1": rec.pos,
                    "pos_2": -1,
                    "num_matches": num_matches,
                    "num_deletions": num_deletions,
                }

                fastq_data.append(read_data)
    return fastq_data


def get_read_data_from_unmerged_sam(unmerged_sam):
    """
    Get data from SAM file of unmerged reads

    Parameters:
        unmerged_sam (str): Path to input SAM file

    Returns:
        fastq_data (list): List of dictionaries where each element
                            in the list corresponds to a line in
                            the SAM file
    """
    fastq_data = []
    for rec in ParseSam(unmerged_sam):
        if rec.cigar != "*":
            num_matches = extractFromCigar("M", rec.cigar)
            num_deletions = extractFromCigar("D", rec.cigar)
            read_data = {
                "read_name": rec.query,
                "sam_flag": rec.flag,
                "chrom": rec.rnam,
                "alignment_score": rec.align_score,
                "num_mismatches": rec.mismatches,
                "mate_start": rec.pnext,
                "insert_size": rec.tlen,
                "pos": rec.pos,
                "num_matches": num_matches,
                "num_deletions": num_deletions,
            }
            fastq_data.append(read_data)
    return fastq_data


###########################################
# Filter reads based on mapping info
###########################################

def not_close_to_max(pandas_df, diff_from_max_percent):
    """
    Filter sequences that are not within a certain
    percentage to the quality of the highest scoring
    alignment for that read

    Parameters:
        pandas_df (pandas data frame):
        diff_from_max_percent (float): Maximum difference in alignment quality between
                                    a given read and the best alignment for that read
    Return:
        Filtered data frame
    """

    pandas_df["max_score_read_name"] = pandas_df.groupby("read_name")["alignment_score"].transform(max)
    pandas_df["diff_percent_max_round"] = np.round(
        pandas_df["max_score_read_name"] * diff_from_max_percent
    )
    pandas_df["score_max_diff"] = pandas_df["max_score_read_name"] - pandas_df["alignment_score"]
    return pandas_df[pandas_df["score_max_diff"] <= pandas_df["diff_percent_max_round"]]


def min_num_of_matches(pandas_df, min_match):
    """
    Filter out reads with less than a certain number of reads matching

    Parameters:
        pandas_df (pandas data frame):
        min_match (int)

    Returns:
        (pandas data frame): Filtered data frame
    """
    return pandas_df[pandas_df["num_matches"] > min_match]


def too_many_mismatch(pandas_df, mismatch_ratio):
    """
    Filter out alignments that have too many mismatches

    Parameters:
        pandas_df (pandas data frame): Pandas data frame with read information
        mismatch_ratio (float): Ratio of mismatches to matches

    Returns:
        (pandas data frame): Filtered pandas data frame
    """
    pandas_df["mismatch_ratio"] = pandas_df["num_mismatches"].astype(float) / pandas_df[
        "num_matches"
    ].astype(float)
    return pandas_df[pandas_df["mismatch_ratio"] < mismatch_ratio]


def too_many_deletions(pandas_df, num_deletions):
    """
    Filter out alignments that have too many deletions
    """
    return pandas_df[pandas_df["num_deletions"] < num_deletions]


def map_to_one_chrom(pandas_df):
    """
    Filter out alignments that only map to one chromosome
    """
    mapper = pandas_df.groupby("read_name")["chrom"].nunique().to_dict()
    pandas_df["unique_chroms"] = pandas_df["read_name"].map(mapper)
    return pandas_df[pandas_df["unique_chroms"] == 1]

# ----------------------------------------------


def filter_merged_sam(
        merged_sam,
        diff_from_max_percent=0.025,
        min_match=75,
        mismatch_ratio=0.025,
        num_deletions=3,
):
    """
    Filter reads from the SAM file of merged reads based on certain filtering stats

    Parameters:
         merged_sam (str): Path to SAM file of merged reads
         diff_from_max_percent (float): Maximum difference in alignment quality between
                                    a given read and the best alignment for that read
        min_match (int): Minimum number of bases matching between the query and reference
        mismatch_ratio (float): Maximum ratio of matches to mismatches between query
                                and reference
        num_deletions (float): Maximum number of deletions

    Returns:
        df (pandas data frame): Data frame with info for sequences that survived read rescue
    """
    fastq_data = get_read_data_from_merged_sam(merged_sam)

    if fastq_data:
        # Transform the fastq_data into a pandas dataframe
        df = pd.DataFrame(fastq_data)

        df = not_close_to_max(df, diff_from_max_percent)

        # Filter out reads with less than a certain number of bases matching.
        df = min_num_of_matches(df, min_match)

        df = too_many_mismatch(df, mismatch_ratio)

        # Filter out reads with too many deletions.
        df = too_many_deletions(df, num_deletions)

        # Finally keep reads that have only one annotated chromosome.
        df = map_to_one_chrom(df)
        return df
    else:
        return False


def read_rescue_unmerged(
        unmerged_sam, diff_from_max_percent=0.025, min_match=100, num_deletions=3,
        mismatch_ratio=0.025
):
    """"
    Preform 'read rescue' for unmerged (paired-end) sequences

    Parameters:
        unmerged_sam (str): Path to SAM file of unmerged reads
        diff_from_max_percent (float): Maximum difference in alignment quality between
                                    a given read and the best alignment for that read
        min_match (int): Minimum number of bases matching between the query and reference
        num_deletions (float): Maximum number of deletions
        mismatch_ratio (float): Ratio of mismatches to matches

    Returns:
        df_pair (pandas data frame): Data frame with data for reads that survived
                                    read rescue

    """
    read_pair = []
    fastq_data = []
    unmerged_read_data = get_read_data_from_unmerged_sam(unmerged_sam)

    for read_data in unmerged_read_data:
        read_pair.append(read_data)
        if read_pair:
            if (
                    abs(read_data["insert_size"]) != abs(read_pair[0]["insert_size"])
                    or read_pair[0]["chrom"] != read_data["chrom"]
                    or read_data["read_name"] != read_pair[0]["read_name"]
            ):
                continue
            read_pair.append(read_data)
            total_num_matches = read_pair[0]["num_matches"] + read_pair[1]["num_matches"]
            total_num_mismatches = (
                    read_pair[0]["num_mismatches"] + read_pair[1]["num_mismatches"]
            )
            total_num_deletions = (
                    read_pair[0]["num_deletions"] + read_pair[1]["num_deletions"]
            )
            total_alignment_score = (
                    read_pair[0]["alignment_score"] + read_pair[1]["alignment_score"]
            )
            fastq_data.append(
                {
                    "read_name": read_data["read_name"],
                    "num_mismatches": total_num_mismatches,
                    "chrom": read_data["chrom"],
                    "num_matches": total_num_matches,
                    "num_deletions": total_num_deletions,
                    "insert_size": abs(read_data["insert_size"]),
                    "alignment_score": total_alignment_score,
                    "pos_1": read_pair[0]["pos"],
                    "pos_2": read_pair[1]["pos"],
                }
            )
            read_pair = []

        else:
            if read_data["sam_flag"] & 2 == 2 and read_data["sam_flag"] & 64 == 64:
                read_pair.append(read_data)

    df_pair = pd.DataFrame(fastq_data)

    # Filter reads with big insert size (gt 999)
    df_pair = df_pair[df_pair["insert_size"] < 1000]

    # Filter out reads less than 2.5% max score.
    df_pair = not_close_to_max(df_pair, diff_from_max_percent)

    # Filter out reads with less than 100 bp matching.
    df_pair = min_num_of_matches(df_pair, min_match)

    # Filter out reads with too many mismatches (more than 1 mismatch every 40 bp).
    df_pair = too_many_mismatch(df_pair, mismatch_ratio)

    # Filter out reads with too many deletions.
    df_pair = df_pair[df_pair["num_deletions"] <= num_deletions]

    # Finally keep reads that have only one annotated chromosome.
    df_pair = map_to_one_chrom(df_pair)


def combine_df(df1, df2):
    """
    Combine 2 pandas data frames

    Parameters:
        df1 (pandas data frame): Data frame 1
        df2 (pandas data frame): Data frame 2

    Returns:
        (pandas df): Concatenated data frame

    """
    return pd.concat([df1, df2])


def bam2sam(bam_file):
    """
    Convert BAM files to SAM files

    Parameters:
        bam_file (str): path to BAM file

    Returns:
        f.name (str): path to SAM file
    """
    with tempfile.NamedTemporaryFile(delete=False) as f:
        subprocess.run(["samtools", "view", bam_file], stdout=f)
        temp_files_list(f.name)
        return f.name


def get_start_and_stop(in_sam, min_mapq):
    """
    Get the start and stop positions for reads in a SAM file

    Parameters:
        in_sam (str): Path to input SAM file
        min_mapq (int): Minimum mapq score

    Returns:
        seen_starts (list):
    """
    seen_starts, seen_ends = set(), set()
    for rec in ParseSam(in_sam):
        if "SECONDARY" in decompose_sam_flag(rec.flag) or "SUPPLEMENTARY" in decompose_sam_flag(rec.flag):
            continue
        if "READ1" in decompose_sam_flag(rec.flag) and rec.tlen >= 0:
            if rec.mapq >= min_mapq:
                seen_starts.add("{}:{}".format(rec.rnam, rec.pos))
                seen_ends.add("{}:{}".format(rec.rnam, rec.pos + rec.tlen - 1))
    return seen_starts, seen_ends


def read_rescue_update(df, in_sam, min_mapq):
    """
    Preform 'read rescue update'

    Parameters:
        df (pandas data frame): data frame from combine_df()
        in_sam (str): path to input BAM file
        min_mapq (int): Minimum mapq score

    Returns:
        out_sam.name (str): Path to output SAM
    """
    ok_reads = set()
    with tempfile.NamedTemporaryFile(delete=False) as out_sam:
        starts_and_stops = get_start_and_stop(in_sam, min_mapq)
        seen_starts, seen_ends = starts_and_stops[0], starts_and_stops[1]
        group_df_dict = {g: gdf for g, gdf in df.groupby("read_name")}

        for line in open(in_sam, "r"):
            if line.startswith("@"):
                out_sam.write(line.encode())
                continue
            sam_fields = line.strip().split("\t")
            header = sam_fields[0]
            chrom = sam_fields[2]
            mapq_score = int(sam_fields[4])
            # Check if read has been added already. if it has write and
            if header in ok_reads:
                out_sam.write(line.encode())
                continue
            elif mapq_score >= min_mapq:
                out_sam.write(line.encode())
                ok_reads.add(header)
                continue
            else:
                try:
                    sub_df = group_df_dict[header]
                except:
                    continue
                do_not_add = False
                if not sub_df.empty:
                    if not sub_df["chrom"].iloc[0] == chrom:
                        continue
                    to_add_start = []
                    to_add_end = []
                    for i, row in sub_df.iterrows():
                        coord_start = "{}:{}".format(row["chrom"], row["pos_1"])
                        if int(row["pos_2"]) != -1:
                            coord_end = "{}:{}".format(row["chrom"], row["pos_2"])
                        else:
                            coord_end = False
                        if coord_start in seen_starts:
                            do_not_add = True
                        else:
                            if coord_end:
                                if coord_end in seen_ends:
                                    do_not_add = True
                        if not do_not_add:
                            to_add_start.append(coord_start)
                            if coord_end:
                                to_add_end.append(coord_end)
                else:
                    do_not_add = True
                if not do_not_add:
                    for coord in to_add_start:
                        seen_starts.add(coord)
                    for coord in to_add_end:
                        seen_ends.add(coord)
                    ok_reads.add(header)
                    out_sam.write(line.encode())
        return out_sam.name


def count_chrom(sam, hom, het):
    """
    Count the number of reads that map to the homogametic and
    heterogametic elements

    Parameters:
        sam (str): path to SAM File
        hom (str): FASTA identifier for homogametic element
        het (str): FASTA identifier for heterogametic element

    Returns:
        hom_counts (int): Count of reads that mapped to homogametic element
        het_counts (int): Count of reads that mapped to heterogametic element
    """
    het_counts, hom_counts = 0, 0
    for rec in ParseSam(sam):
        if rec.rnam == hom:
            hom_counts += 1
        elif rec.rnam == het:
            het_counts += 1
    return hom_counts, het_counts


def calculate_stats(homogametic_counts, heterogametic_counts):
    """
    Calculate statistics for sequence counts

    Parameters:
        homogametic_counts (int): Number of reads that mapped to homogametic element
        heterogametic_counts (int): Number of reads that mapped to heterogametic element

    Returns:
        round(prop, 3) (float): Proportion of heterogametic to homogametic reads
                                rounded to 3 decimal places
        round(ci, 3) (float): Confidence interval for read proportions
    """
    prop = heterogametic_counts / (homogametic_counts + heterogametic_counts)
    ci = 1.96 * math.sqrt(
        (prop * (1 - prop)) / (heterogametic_counts + homogametic_counts)
    )
    return round(prop, 3), round(ci, 3)
