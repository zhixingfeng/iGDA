# iGDA

## Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Usage](#usage)
- [Demo](#demo)
- [License](./LICENSE)

# Overview 
Cellular genetic heterogeneity is common in many biological conditions including cancer, microbiome, co-infection of multiple pathogens. Detecting and phasing minor variants, which is to determine whether multiple variants are from the same haplotype, play an instrumental role in deciphering cellular genetic heterogeneity, but are still difficult because of technological limitations. Recently, long-read sequencing technologies, including those by Pacific Biosciences and Oxford Nanopore, have provided an unprecedented opportunity to tackle these challenges. However, high error rates make it difficult to take full advantage of these technologies. To fill this gap, we introduce iGDA, an open-source tool that can accurately detect and phase minor single-nucleotide variants (SNVs), whose frequencies are as low as 0.2%, from raw long-read sequencing data. We also demonstrated that iGDA can accurately reconstruct haplotypes in closely-related strains of the same species (divergence >= 0.011%) from long-read metagenomic data. Our approach, therefore, presents a significant advance towards the complete deciphering of cellular genetic heterogeneity. 

# System Requirements

## Hardware Requirements

For optimal performance, we recommend a computer with the following specs:

RAM: 64+ GB

CPU: 16+ cores

## Software Requirements

### OS Requirements

The iGDA package is tested on *Linux* operating systems.

Linux: CentOS Linux release 7.6.1810 (Core)

# Installation Guide

## Install via *Conda* (recommended, need *Conda*)
conda install -c zhixingfeng igda

## Compile from source code (not recommended)
### Install dependencies
1. xgboost 0.90
2. bioawk 1.0
3. samtools 1.10

### Compile source code (GCC version >= 5, or other C++ compilers supporting c++14 standard)
Download the source code of iGDA from *Release*, unzip it, enter the directory and type "Make". Add the "bin" directory to your PATH or create a soft link to ./bin/igda in a directory that the system can find.

### Download script 
Download https://github.com/zhixingfeng/shell/archive/0.9.3.tar.gz, unzip it and add the folder to your PATH.

# Usage 
**To detect low-frequency SNVs, use:**

(PacBio data) igda_pipe_detect_pb infile(bam or sam file) reffile contextmodel outdir

(Nanopore data) igda_pipe_detect_ont infile(bam or sam file) reffile contextmodel outdir

**To detect combination of  low-frequency SNVs (phasing), use:**

(PacBio data) igda_pipe_phase_pb indir(output of igda_pipe_detect_pb) reffile outdir

(Nanopore data) igda_pipe_phase_ont indir(output of igda_pipe_detect_ont) reffile outdir

**Please note that reffile is the reference fasta file. Current version assumes there is only one contig in samfile and reffile.**

contextmodel is the context effect model trained on independent data. They can be download in https://github.com/zhixingfeng/igda_contextmodel

**Output format:**

**For detecting low-frequency SNVs, realign.var in outdir is the final result.**

Column 1 is the locus of detected SNVs.

Column 2 is the alternative base of detected SNVs.

The other columns are reserved for internal use


**For detecting combinations of low-frequency SNVs, realign.ann.tested.ft.count.ft.assembled.count.nc.ft in outdir is the final result.**

Each row is a contig

Column 1 is the SNVs of the contigs. It is encoded, for each integer x, floor(x/4) = 0-based locus, and x modulo 4 = base

Column 2 is start locus (0-based)

Column 3 is end locus (0-based)

Column 4 is number of reads aligned to the contig

Column 5 is coverage of the contig

The other columns are reserved for internal use

The file realign.ann.tested.ft.count.ft.assembled.count.nc.ft.sam in outdir can be loaded directly into IGV.

**For any questions, contact zxfeng.thu@gmail.com or zhixing.feng@mssm.edu**

# Demo

Example data and code can be download from the following links:

1. [A mixture of 186 *Bordetella* spp. samples (PacBio sequencing)](https://www.dropbox.com/sh/uusx8modggni96m/AAAxjKEa7YdG-HYpKnzousKBa?dl=0)
2. [A mixture of 155 *Escherichia coli* samples (PacBio sequencing)](https://www.dropbox.com/sh/uusx8modggni96m/AAAxjKEa7YdG-HYpKnzousKBa?dl=0)
3. [A mixture of 65 *Klebsiella pneumoniae* samples (Oxford Nanopore sequencing)](https://www.dropbox.com/sh/uusx8modggni96m/AAAxjKEa7YdG-HYpKnzousKBa?dl=0)

Alternative links:

link: https://pan.baidu.com/s/1nGk7ptqQyJq56gk-cmtPuw  

password: rbw4


The run time of demo 1~3 is expected to be < 1 hour with 16 cores, but may vary depending on the specific machines and their working conditions. 
