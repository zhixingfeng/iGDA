# iGDA

## Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Usage](#usage)
- [Demo](#demo)
- [Maintainer](#maintainer)
- [License](./LICENSE)

# Overview 
Cellular genetic heterogeneity is common in many biological conditions including cancer, microbiome, co-infection of multiple pathogens. Detecting and phasing minor variants, which is to determine whether multiple variants are from the same haplotype, play an instrumental role in deciphering cellular genetic heterogeneity, but are still difficult because of technological limitations. Recently, long-read sequencing technologies, including those by Pacific Biosciences and Oxford Nanopore, have provided an unprecedented opportunity to tackle these challenges. However, high error rates make it difficult to take full advantage of these technologies. To fill this gap, we introduce iGDA, an open-source tool that can accurately detect and phase minor single-nucleotide variants (SNVs), whose frequencies are as low as 0.2%, from raw long-read sequencing data. We also demonstrated that iGDA can accurately reconstruct haplotypes in closely-related strains of the same species (divergence >= 0.011%) from long-read metagenomic data. Our approach, therefore, presents a significant advance towards the complete deciphering of cellular genetic heterogeneity. 

# System Requirements

## Hardware Requirements

For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB

CPU: 8+ cores

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
**To detect minor SNVs, use:**

(PacBio data) igda_pipe_detect -m pb bamfile reffile contextmodel outdir

(Nanopore data) igda_pipe_detect -m ont bamfile reffile contextmodel outdir

**To phase minor SNVs, use:**

(PacBio data) igda_pipe_phase -m pb indir(outdir of igda_pipe_detect) reffile outdir

(Nanopore data) igda_pipe_phase -m ont indir(outdir of igda_pipe_detect) reffile outdir

**Please note**

bamfile is the aligned bam file (indexed).

reffile is the reference fasta file.

contextmodel is the context effect model trained on independent data. They can be download in https://github.com/zhixingfeng/igda_contextmodel

**Output format:**

**For detecting minor SNVs, detected_snv.vcf in outdir is the final result.**

**For phasing minor SNVs, contigs.sam, contigs.fa, and contigs.ann are the final results.**

In the contigs.ann file, each row is a contig.

Column 1 is chromosome name.

Column 2 is the SNVs of the contig. It is encoded, for each integer x, floor(x/4) = 0-based locus, and x modulo 4 = base (0=A, 1=C, 2=G, 3=T)

Column 3 is start locus (0-based)

Column 4 is end locus (0-based)

The other columns are reserved for internal use

# Demo

Example data and code can be download from the following link, which includes:

1. A mixture of 186 *Bordetella* spp. samples (PacBio sequencing)
2. A mixture of 155 *Escherichia coli* samples (PacBio sequencing)
3. A mixture of 65 *Klebsiella pneumoniae* samples (Oxford Nanopore sequencing)
4. A mixture of 11 Borrelia burgdorferi strains and 744 other species to mimic a metagenome (PacBio sequencing)

https://www.dropbox.com/sh/umi03t3eendcktf/AABhXy2cNfQAr0wBike87Kc0a?dl=0

The run time of demo 1~3 is expected to be several minutes with 16 cores, but may vary depending on the specific machines and their working conditions. It takes about 5 hours to finish demo 4 with 32 cores due to its size.

# Maintainer

Zhixing Feng (冯智星)

**For any questions, contact zxfeng.thu@gmail.com or zhixing.feng@mssm.edu**

