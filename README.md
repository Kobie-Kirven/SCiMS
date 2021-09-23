
<img src="https://github.com/Kobie-Kirven/SCiMS/blob/main/docs/_static/logo.png" width="300">
<h1>Sex Calling for Metagenomic Sequences</h1>

+ [What is SCiMS?](#What-is-SCiMS?:)
+ [Installation](#Installation:)
+ [Usage](#Usage:)
  + [Building Indicies](#Building-Indices:)
  + [Determine Sex](#Determine-Sex:)
+ [Install Requirements](#Install-Requirements:)
+ [Quick Tutorial](#Quick-Tutorial:)
---
# What is SCiMS?:
SCiMS (Sex Calling for Metagenomic Sequences) is a tool for determining 
host sex using shotgun metagenomic data. Currently, this tool supports sex
determination for host organisms that contain two sex chromosomes.
(Ex. X and Y in humans). Briefly, SCiMS works by aligning shotgun metagenomic
data to a reference genome and calculating the proportion of reads that map 
to the heterogametic chromosome versus both sex chromosomes together. 
## Installation:

<h4>Requirements:</h4>
- bwa
- bowtie2
- Flash
- Samtools

All of the required tools can be installed with conda (see [Install Requirements](#Install Requirements:))

SCiMS can be easily installed with:
```bash
pip3 install git+https://github.com/Kobie-Kirven/SCiMS
```
## Usage:
- ### Building Indices:
  The first step in using SCiMS is to build indices of the reference genome for
  BWA and Bowtie2. The ```build-index``` module makes it easy to build the indices
  for both tools with one command. 
- #### Requirements:
  * Reference genome in FASTA format (can be gzip)
  * Name of output index
    ```bash
    scims build-index -r <reference_genome> -o <index_name>
    ```
- ### Determine Sex:
## Install Requirements:
To install the required tools with conda, run the following code. 
```bash
conda install -c bioconda bwa
conda install -c bioconda flash
conda install -c bioconda bowtie2
conda install -c bioconda samtools
```
## Quick Tutorial: