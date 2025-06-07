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

### Basic Command
```bash
python cat.py --project_dir /path/to/output \
--seq_dir /path/to/fastq_files \
--length 150 \
--gtf /path/to/annotation.gtf \
--bt2 /path/to/bowtie2_index \
--star /path/to/star_index \
--strand F2R1
```

### Parameters

- `--project_dir`: Output directory for results
- `--seq_dir`: Directory containing input FASTQ files
- `--length`: Average read length (default: 150)
- `--gtf`: Path to GTF annotation file
- `--bt2`: Path to Bowtie2 reference index
- `--star`: Path to STAR reference index
- `--strand`: Strand specificity (F2R1 or F1R2)
- `--threads`: Number of threads running cat(default: 4)

## Workflow

1. **Read Alignment**: CAT first aligns reads to the reference genome using STAR
2. **Filtering**: Removes low-quality alignments and PCR duplicates
3. **CRISPR Detection**: Identifies potential CRISPR-RNA interactions
4. **Classification**: Categorizes different types of CRISPR events
5. **Quantification**: Calculates statistics for each event type
6. **Visualization**: Generates plots and summary reports

## Output Files
- `project_dir/`
   - `MappingResults/`: Contains first two alignment files (BAM)
   - `ProcessedResults/`: Breakpoint information and quantitative results
   - `circ_raw_result.tsv`: Raw results not filtered by breakpoint perimeter overlays
   - `circleRNA_predict_result.tsv`: The final output of CAT
  

## Troubleshooting

- **Memory Issues**: For large datasets, increase available memory or reduce threads
- **Alignment Errors**: Ensure reference indices are correctly built 
- **Missing Dependencies**: Verify all tools are properly installed and in PATH


