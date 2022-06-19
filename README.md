
<img src="https://github.com/Kobie-Kirven/SCiMS/blob/main/docs/_static/logo.png" width="300">
<h1>Sex Calling for Metagenomic Sequences</h1>

- [What is SCiMS?](#What-is-SCiMS?)
- [Installation](#Installation)
- [Usage](#Usage)
- [Install Requirements](#Install-Requirements)
- [Quick Tutorial](#Quick-Tutorial)
---
## What is SCiMS?
SCiMS (Sex Calling for Metagenomic Sequences) is a tool for determining 
host sex using shotgun metagenomic data. Currently, this tool supports sex
determination for host organisms that contain two sex chromosomes.
(Ex. X and Y in humans). Briefly, SCiMS works by aligning shotgun metagenomic
data to a reference genome and comparing the coverage of the homogametic chromosome
to the coverage of the autosomes. 
## Installation

### Requirements
- bowtie2==2.4.4
- Samtools==1.13

All of the required tools can be installed with conda (see [Install Requirements](#Install-Requirements))

SCiMS can be easily installed with:
```bash
pip3 install git+https://github.com/Kobie-Kirven/SCiMS
```
## Usage

- SCiMS is designed to work with both single and paired-end metagenomic data. The command-line usage is:
    ```
    usage: scims [-h] [-v] [--bowtie-index INDEX] [--ref-genome REF] [--scaffold-names SCAFFOLD] [-1 FORWARD] [-2 REVERSE] [-s SINGLE] [--x HOMOGAMETIC]
             [--y HETEROGAMETIC] --o OUTPUT [--t THREADS] [--from-sam FROM_SAM] [--get-scaffold-lengths LENGTHS] [--scaffold-lengths SCAF_LENGTHS]
   Sex Calling for Metagenomic Sequences
   options:
   -h, --help            show this help message and exit
   -v, --version         show program's version number and exit
   --bowtie-index INDEX  Path to Bowtie2 index
   --ref-genome REF      Reference genome in FASTA or FASTA.gz format
   --scaffold-names SCAFFOLD
                           Path to file with scaffold names
   -1 FORWARD            Forward reads in FASTA or FASTQ format (PE mode only)
   -2 REVERSE            Reverse reads in FASTA or FASTQ format (PE mode only)
   -s SINGLE             Single-end reads in FASTA or FASTQ format (SE mode only)
   --x HOMOGAMETIC       ID of homogametic sex chromosome (ex. X)
   --y HETEROGAMETIC     ID of heterogametic sex chromesome (ex. Y)
   --o OUTPUT            Output plot prefix
   --t THREADS           Number of threads to use
   --from-sam FROM_SAM   Use sam file instead of preforming the alignment
   --get-scaffold-lengths LENGTHS
                           Determine scaffold lengths
   --scaffold-lengths SCAF_LENGTHS
                           Path to file with scaffold lengths
    ```

## Install Requirements
To install the required tools with conda, run the following code. 
```bash
conda install -c bioconda bowtie2
conda install -c bioconda samtools
```
Note: For some reason, the current version of Bowtie2 has an issue with tbb. To fix this, 
run 
```bash
conda install tbb=2020.2
```
## Quick Tutorial
This tutorial is intended to ensure that your SCiMS installation is working correctly:

### Directly from SAM file 
1. Download test files:
   ```bash
   wget https://github.com/Kobie-Kirven/SCiMS/raw/main/test_data/male_test.sam
   wget https://github.com/Kobie-Kirven/SCiMS/raw/main/test_data/scaffold_lengths.txt
   wget https://github.com/Kobie-Kirven/SCiMS/raw/main/test_data/scaffolds.txt
   ```
   ```
5. Run SCiMS in alignment-free mode:
   ```text
   scims_test]$ scims --scaffold-names scaffolds.txt --x NC_000023.11 --y NC_000024.10 --o test -
-t 1 --from-sam male_test.sam --scaffold-lengths scaffold_lengths.txt "
   ```

6. You should see output similar to the following:
   ```text
   Using SAM file
   Extracting propper alignments:
         A total of 18449 alignments met the criteria
   The average NC_000024.10:NC_000023.11 ratio for the file was 0.5285776962795797
   ```
   <img src="https://github.com/Kobie-Kirven/SCiMS/blob/main/docs/_static/test.png" width="300">