# iGDA

## Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Usage](#usage)
- [Demo](#demo)
- [Citation](#citation)
- [Maintainer](#maintainer)
- [License](./LICENSE)

# Overview 
Cellular genetic heterogeneity is common in many biological conditions including cancer, microbiome, and co-infection of multiple pathogens. Detecting and phasing minor variants, which is to determine whether multiple variants are from the same haplotype, play an instrumental role in deciphering cellular genetic heterogeneity, but are still difficult because of technological limitations. Recently, long-read sequencing technologies, including those by Pacific Biosciences and Oxford Nanopore, have provided an unprecedented opportunity to tackle these challenges. However, high error rates make it difficult to take full advantage of these technologies. To fill this gap, we introduce iGDA, an open-source tool that can accurately detect and phase minor single-nucleotide variants (SNVs), whose frequencies are as low as 0.2%, from raw long-read sequencing data. We also demonstrated that iGDA can accurately reconstruct haplotypes in closely-related strains of the same species (divergence >= 0.011%) from long-read metagenomic data. 

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

## Install via Conda (recommended, need Conda)
```conda install -c zhixingfeng igda```

## Compile from source code (not recommended)
### Install dependencies
1. xgboost 0.90
2. bioawk 1.0
3. samtools 1.10

### Compile source code (GCC version >= 5, or other C++ compilers supporting c++14 standard)
Download the source code of iGDA from *Release*, unzip it, enter the directory and type "Make". Add the "bin" directory to your PATH or create a soft link to ./bin/igda in a directory that the system can find.

### Download script 
Download https://github.com/zhixingfeng/shell/archive/refs/tags/1.0.1.tar.gz, unzip it and add the folder to your PATH.

# Usage 
## Preprocessing:

Convert lower case letters to upper case in reference fasta file:

```fasta2upper infasta outfasta```

Convert wildcard letters to N in reference fasta file:

```fastaclean fastafile outfafile```

Realign reads aligned to the negative strand:

(PacBio data) ```igda_align_pb infile(bam or sam file) reffile outfile nthread```

(PacBio data with no QV, like Sequel or Sequel II) ```igda_align_pb_fa infile(bam or sam file) reffile outfile nthread```

(Nanopore data) ```igda_align_ont infile(bam or sam file) reffile outfile nthread```

Sort and convert realign samfile to bamfile:

```sam2bam samfile nthread```

## Detect minor SNVs:

(PacBio data) ```igda_pipe_detect -m pb bamfile reffile contextmodel outdir```

(Nanopore data) ```igda_pipe_detect -m ont bamfile reffile contextmodel outdir```

**Please note:**

```bamfile``` is the aligned bam file (sorted and indexed).

```reffile``` is the reference fasta file.

```contextmodel``` is the context effect model trained on independent data. They can be download in https://github.com/zhixingfeng/igda_contextmodel

## Phase minor SNVs:

(PacBio data) ```igda_pipe_phase -m pb indir(outdir of igda_pipe_detect) reffile outdir```

(Nanopore data) ```igda_pipe_phase -m ont indir(outdir of igda_pipe_detect) reffile outdir```

## Output format:

For detecting minor SNVs, ```detected_snv.vcf``` in outdir is the final result.

For phasing minor SNVs, ```contigs.sam```, ```contigs.fa```, and ```contigs.ann are``` the final results.

In the ```contigs.ann file```, each row is a contig.

Column 1 is chromosome name.

Column 2 is the SNVs of the contig. It is encoded, for each integer x, floor(x/4) = 0-based locus, and x modulo 4 = base (0=A, 1=C, 2=G, 3=T)

Column 3 is start locus (0-based)

Column 4 is end locus (0-based)

The other columns are reserved for internal use

## Parameter tuning:

It is difficult to find an universally optimal parameter setting. Here are some tips for parameter tuning if the default one does not have the expected performance:

* By default, ```igda_pipe_detect``` discards reads with aligned length < 1000 because it can reduce the impact read-maping ambiguity. Use ```igda_pipe_detect -l 0 ``` instead if the data have lots of aligned reads shorter than 1000.

* In ```igda_pipe_detect```, ```-r``` and ```-c``` are the two major parameters affecting the accuracy. ```-r``` is "Minimal depth for each SNV" and ```-c``` is "Minimal maximal conditional substitution rate". By default, ```igda_pipe_detect``` uses ```-r 25 -c 0.65```, which is a conservative setting aiming to achieve a low false discover rate (FDR). These parameters might have a low sensitivity if the sequencing depth corresponding to a minor SNV is lower or close to 25. It is possible to increase sensitivity by decreasing ```-r```, but FDR might increase. 

# Demo

Example data and code can be download from the following link, which includes:

1. A mixture of 186 *Bordetella* spp. samples (PacBio sequencing)
2. A mixture of 155 *Escherichia coli* samples (PacBio sequencing)
3. A mixture of 65 *Klebsiella pneumoniae* samples (Oxford Nanopore sequencing)
4. A mixture of 11 Borrelia burgdorferi strains and 744 other species to mimic a metagenome (PacBio sequencing)

https://www.dropbox.com/sh/umi03t3eendcktf/AABhXy2cNfQAr0wBike87Kc0a?dl=0

The run time of demo 1~3 is expected to be several minutes with 16 cores, but may vary depending on the specific machines and their working conditions. It might take about several hours to finish demo 4 with 32 cores due to its size.

# Citation
Feng, Z., Clemente, J.C., Wong, B. et al. Detecting and phasing minor single-nucleotide variants from long-read sequencing data. Nat Commun 12, 3032 (2021). 

DOI: https://doi.org/10.1038/s41467-021-23289-4

# Maintainer

Zhixing Feng (冯智星)

**For any questions, contact zxfeng.thu@gmail.com or zhixing.feng@mssm.edu**

