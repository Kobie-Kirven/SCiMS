
<img src="https://github.com/Kobie-Kirven/SCiMS/blob/main/docs/_static/logo.png" width="300">
<h1>Sex Calling for Metagenomic Sequences</h1>

- [What is SCiMS?](#What-is-SCiMS?)
- [Installation](#Installation)
- [Usage](#Usage)
  - [Building Indicies](#Building-Indices)
  - [Determine Sex](#Determine-Sex)
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

<h4>Requirements:</h4>
* bwa 
* bowtie2 
* Flash 
* Samtools

All of the required tools can be installed with conda (see [Install Requirements](#Install Requirements:))

SCiMS can be easily installed with:
```bash
pip3 install git+https://github.com/Kobie-Kirven/SCiMS
```
## Usage
- ### Building Indices
  The first step in using SCiMS is to build indices of the reference genome for
  BWA and Bowtie2. The ```build-index``` module makes it easy to build the indices
  for both tools with one command.
    ```bash
    usage: scims build-index [-h] [-r REFERENCE] [-o OUTPUT]

    optional arguments:
      -h, --help            show this help message and exit
      -r REFERENCE, --reference REFERENCE
                            Host reference genome in FASTA or FASTQ format
                            (default: None)
      -o OUTPUT, --output OUTPUT

    ```
- ### Determine Sex
  The ```determine-sex``` module will give an output of the proportion
  of reads that map to the homogametic chromosome versus both sex chromosomes. 
  This module requires that you already have the BWA/Bowtie2 indices built. 
    ```bash
    usage: scims determine-sex [-h] [-i INDEX] [-1 FORWARD] [-2 REVERSE] [-t THREADS] [-hom HOMOGAMETIC]
                               [-het HETEROGAMETIC]
    
    optional arguments:
      -h, --help            show this help message and exit
      -i INDEX, --index-name INDEX
                            Name of bowtie2 and bwa index (default: None)
      -1 FORWARD, --forward-reads FORWARD
                            Forward reads in fasta or fastq format (default: None)
      -2 REVERSE, --reverse-reads REVERSE
                            Reverse reads in fasta or fastq format (default: None)
      -t THREADS, --number-of-threads THREADS
                            Number of threads to use (default: None)
      -hom HOMOGAMETIC, --homogametic HOMOGAMETIC
                            ID of homogametic sex chromesome (ex. X) (default: None)
      -het HETEROGAMETIC, --heterogametic HETEROGAMETIC
                            ID of heterogametic sex chromesome (ex. Y) (default: None)
    ```
## Install Requirements
To install the required tools with conda, run the following code. 
```bash
conda install -c bioconda bwa
conda install -c bioconda flash
conda install -c bioconda bowtie2
conda install -c bioconda samtools
```
## Quick Tutorial