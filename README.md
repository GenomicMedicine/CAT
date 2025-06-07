![logo_l](https://github.com/user-attachments/assets/71ca4f1d-d4d8-4713-aa25-a30215245d6d)

## Overview
CAT is a comprehensive tool for detecting and analyzing all types of circular RNA detection. CAT does not rely on annotation of splicing sites and uses filters driven by biological foundation models to reduce false positives.

## Benchmark
We evaluated the CAT method in identifying c-circRNAs – its performance on c-circRNAs can be used as a proxy for its performance on all types of circRNAs. For this analysis, we leveraged data and results from a recent benchmarking study* that compared 16 existing circRNA detection tools
CAT detected 92.7% of the results that had been experimentally verified in previous studies, and an additional 7.2% of novo circRNAs that were not detected by any of the 16 previous methods.
![fig2](https://github.com/user-attachments/assets/51131639-592b-49e3-9400-1da35b71b457)

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

## Reference
* Vromman, M., et al., Large-scale benchmarking of circRNA detection tools reveals large differences in sensitivity but not in precision. Nat Methods, 2023. 20(8): p. 1159-1169.
