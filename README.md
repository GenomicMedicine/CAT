# CAT (CricRNA All Types) Tutorial
![fig1](https://github.com/user-attachments/assets/bc17357d-4907-4cbf-b431-5ccf863488b3)

## Overview
CAT is a comprehensive tool for detecting and analyzing all types of CRIC-RNA detection. This tutorial will guide you through the installation, setup, and usage of CAT.

## Dependencies

### Python Packages
- tqdm
- pandas
- arrow
- pyensembl
- pysam
- biopython (for SeqIO)

### External Tools
- bedtools
- STAR
- bowtie2

## Installation

Clone the repository
```bash
git clone https://github.com/username/CAT.git
cd CAT
```
Create a conda environment (recommended)
```bash
conda create -n cat_env python=3.8
conda activate cat_env
```
Install Python dependencies
```bash
pip install tqdm pandas pyarrow pyensembl pysam biopython
```
Install external tools
```bash
conda install -c bioconda bedtools star bowtie
```

## Data Preparation

1. Prepare your sequencing data (FASTQ files)
2. Prepare reference files(Recommended ensembel version of the genome and annotations):
   - GTF annotation file
   - Genome FASTA file
   - Build STAR and Bowtie2 indices

Build STAR index
```bash
STAR --runMode genomeGenerate --genomeDir star_index \
--genomeFastaFiles genome.fa --sjdbGTFfile annotation.gtf \
--runThreadN 8
```
Build Bowtie2 index
```bash
bowtie2-build genome.fa bowtie2_index
```

## Usage

See wiki
