import pandas as pd
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

def create_circ_reference():
    # Read circRNA prediction result file
    df = pd.read_csv('/data3/H23/ribo-seq-a549/A549_circleRNA_predict_result.tsv', sep='\t')
    
    # Create list of fasta records
    records = []
    for _, row in df.iterrows():
        ref_seq = row['Reference'][100:200]  # Get sequence from position 100-200
        if len(ref_seq) > 0:  # Ensure sequence is not empty
            record = SeqRecord(
                Seq(ref_seq),
                id=row['CircID'],
                description=""
            )
            records.append(record)
    
    # Write to fasta file
    output_fasta = "a549_noncanonial_circle_refer.fasta"
    SeqIO.write(records, output_fasta, "fasta")
    return output_fasta

def quality_control(input_fastq):
    """
    Perform quality control on input fastq file
    """
    output_dir = os.path.dirname(input_fastq)
    cmd = f"trim_galore -j 5 -q 20 --no_report_file --length 20 --max_length 35 -o {output_dir} {input_fastq}"
    subprocess.run(cmd, shell=True)
    return f"{os.path.splitext(input_fastq)[0]}_trimmed.fq"

def filter_reads(input_fastq):
    """
    Filter rRNA, tRNA, DNA and RNA sequentially
    """
    # Filter rRNA
    cmd_rRNA = (f"bowtie2 -p 100 -x /home/lab/lk/code/bowtie2refer/human_rRNA/rRNA "
                f"--un-gz rmrRNA.fq.gz -U {input_fastq} "
                f"-S rRNA.sam --score-min C,0,0")
    subprocess.run(cmd_rRNA, shell=True)
    subprocess.run("rm rRNA.sam", shell=True)

    # Filter tRNA
    cmd_tRNA = (f"bowtie2 -p 100 -x /home/lab/lk/code/bowtie2refer/human_tRNA/tRNA "
                f"--un rmtRNA.fq -U rmrRNA.fq.gz "
                f"-S tRNA.sam --score-min C,0,0")
    subprocess.run(cmd_tRNA, shell=True)
    subprocess.run("rm tRNA.sam", shell=True)

    # Filter DNA
    cmd_DNA = (f"bowtie2 -p 100 -x /home/lab/lk/code/bowtie2refer/dnarefer "
               f"--un unmapped1.fq -U rmtRNA.fq "
               f"-S DNA.sam --score-min C,0,0")
    subprocess.run(cmd_DNA, shell=True)
    subprocess.run("rm DNA.sam", shell=True)

    # Filter RNA
    cmd_RNA = (f"bowtie2 -p 100 -x /home/lab/lk/code/bowtie2refer/cdnarefer "
               f"--un unmapped2.fq -U unmapped1.fq "
               f"-S RNA.sam --score-min C,0,0")
    subprocess.run(cmd_RNA, shell=True)
    subprocess.run("rm RNA.sam", shell=True)

    return "unmapped2.fq"

def merge_fastq_files():
    """Merge three SRR files"""
    output_file = "merged_riboseq.fastq"
    srr_files = ["SRR9332878.fastq", "SRR9332879.fastq", "SRR9332880.fastq"]
    
    with open(output_file, 'w') as outfile:
        for srr_file in srr_files:
            file_path = os.path.join("/data3/H23/ribo-seq-a549/converted_fasta", srr_file)
            # First perform quality control
            trimmed_file = quality_control(file_path)
            # Then write to merged file
            with open(trimmed_file) as infile:
                outfile.write(infile.read())
    
    return output_file

def process_sam_file(sam_file, allow_mismatch=True):
    """Process SAM file and return matched CircIDs and their supporting read counts"""
    circ_reads = {}
    
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            
            fields = line.strip().split('\t')
            if fields[2] == '*':  # Unmatched
                continue
                
            # Check NM:i tag (mismatch count)
            nm_count = 0
            for tag in fields[11:]:
                if tag.startswith('NM:i:'):
                    nm_count = int(tag.split(':')[2])
                    break
            
            if allow_mismatch:
                if nm_count <= 1:  # Allow 1 mismatch
                    circ_reads[fields[2]] = circ_reads.get(fields[2], 0) + 1
            else:
                if nm_count == 0:  # Only allow perfect match
                    circ_reads[fields[2]] = circ_reads.get(fields[2], 0) + 1
    
    return circ_reads

def add_circ_type():
    """
    Look up CircType from original CSV file and add to perfect_matches.txt
    """
    # Read original circRNA prediction result file, keep only needed columns
    circ_df = pd.read_csv('/data3/H23/ribo-seq-a549/A549_circleRNA_predict_result.tsv', 
                         sep='\t',
                         usecols=['CircID', 'CircType'])
    
    # Convert CircID and CircType to dictionary for lookup
    circ_type_dict = dict(zip(circ_df['CircID'], circ_df['CircType']))
    
    with open("perfect_matches.txt", 'r') as infile, open("perfect_matches_with_type.txt", 'w') as outfile:
        # Write new header
        header = infile.readline().strip()
        outfile.write(f"{header}\tCircType\n")
        
        # Process each line
        for line in infile:
            fields = line.strip().split('\t')
            circ_id = fields[0]
            # Look up CircType from dictionary, mark as "Unknown" if not found
            circ_type = circ_type_dict.get(circ_id, "Unknown")
            # Write all fields plus CircType
            outfile.write(f"{line.strip()}\t{circ_type}\n")

    return "perfect_matches_with_type.txt"

def main():
    # Merge fastq files and perform quality control
    merged_fastq = merge_fastq_files()
    
    # Filter reads
    filtered_fastq = filter_reads(merged_fastq)
    
    # Create circRNA reference sequence
    ref_fasta = create_circ_reference()
    
    # Build index
    subprocess.run(f"bowtie2-build --threads 100 {ref_fasta} noncanonicalcirc", shell=True)
    
    # Align to circRNA reference sequence (allow 1 mismatch)
    subprocess.run(
        f"bowtie2 -p 100 -x noncanonicalcirc --un unmappedcirc.fq -U {filtered_fastq} "
        f"-S ribocircRNA.sam --score-min L,0,-1", 
        shell=True
    )
    
    # Process results
    perfect_matches = process_sam_file('ribocircRNA.sam', allow_mismatch=False)
    one_mismatch = process_sam_file('ribocircRNA.sam', allow_mismatch=True)
    
    # Output results
    with open('perfect_matches.txt', 'w') as f:
        f.write("CircID\tSupport_Reads\n")
        for circ_id, count in perfect_matches.items():
            f.write(f"{circ_id}\t{count}\n")
            
    with open('one_mismatch.txt', 'w') as f:
        f.write("CircID\tSupport_Reads\n")
        for circ_id, count in one_mismatch.items():
            f.write(f"{circ_id}\t{count}\n")

    # Add CircType information
    final_result = add_circ_type()
    print(f"Processing complete, results saved in: {final_result}")

if __name__ == "__main__":
    main()