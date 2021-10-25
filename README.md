
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
data to a reference genome and calculating the proportion of reads that map 
to the heterogametic chromosome versus both sex chromosomes together. 
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
    usage: scims [-h] [-v] -i INDEX -r REF -sca SCAFFOLD [-1 FORWARD] [-2 REVERSE] [-s SINGLE] [-t THREADS] -hom HOMOGAMETIC
             [-het HETEROGAMETIC] -o OUTPUT

  Sex Calling for Metagenomic Sequences
  
  optional arguments:
    -h, --help          show this help message and exit
    -v, --version       show program's version number and exit
    -i INDEX            Path to Bowtie2 index
    -r REF              Reference genome in FASTA or FASTA.gz format
    -sca SCAFFOLD       Path to file with scaffold names
    -1 FORWARD          Forward reads in FASTA or FASTQ format (PE mode only)
    -2 REVERSE          Reverse reads in FASTA or FASTQ format (PE mode only)
    -s SINGLE           Single-end reads in FASTA or FASTQ format (SE mode only)
    -t THREADS          Number of threads to use
    -hom HOMOGAMETIC    ID of homogametic sex chromosome (ex. X)
    -het HETEROGAMETIC  ID of heterogametic sex chromesome (ex. Y)
    -o OUTPUT           Output plot prefix
    ```

- What do you need to run SCiMS?
  - INDEX: Database created from ```bowtie2-build``` using the appropriate reference genome
  - REF: Reference genome that was used to generate bowtie2 database
  - SCA: Text file containing the FASTA IDs for the scaffolds to be included in the analysis (see example data)
  - FORWARD: For paired-end mode, the file containing the forward shotgun metagenomic sequences
  - REVERSE: For paired-end mode, the file containing the reverse shotgun metagenomic sequences
  - SINGLE: For single-end mode, the file containing the shotgun metagenomic sequences
  - THREADS: The number of threads to use
  - HOMOGAMETIC: The FASTA ID of the sex chromosome in which two coppies denotes a particular genetic sex
  - HETEROGAMETIC: The FASTA ID of the sex chromosome in which one copy denotes a particular genetic sex
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

1. Download test FASTQ files:
   ```bash
   wget https://github.com/Kobie-Kirven/SCiMS/blob/main/test_data/female_10000_1.fa
   wget https://github.com/Kobie-Kirven/SCiMS/blob/main/test_data/female_10000_2.fa
   ```
2. Download human reference genome from NCBI:
   ```bash
   wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
   ```
3. Build reference indices for Bowite2:
   ```text
   bowtie2-build GRCh38_latest_genomic.fna.gz GRCh38_latest_genomic.fna.gz
   ```
4. Get the scaffold names that we want to include in our analysis:
    ```
   zcat GRCh38_latest_genomic.fna.gz | grep NC | cut -d " " -f 1 > chroms.txt
   ```
5. Run SCiMS:
   ```text
   scims -i GRCh38_latest_genomic.fna.gz -r GRCh38_latest_genomic.fna.gz \
   -1 female_10000_1.fa -2 female_10000_2.fa -t 2 -sca chroms.txt \
   -hom "NC_000023.11" -het "NC_000024.10"
   ```
6. You should see output similar to the following:
   ```text
    ***COMING SOON***
   ```
   <img src="https://github.com/Kobie-Kirven/SCiMS/blob/main/docs/_static/female.png" width="300">