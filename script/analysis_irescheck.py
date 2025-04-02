import pandas as pd
import subprocess
from collections import defaultdict

def run_blastn(query_fasta, ires_fasta, output_file):
    """Run BLASTN alignment"""
    cmd = f"blastn -query {query_fasta} -subject {ires_fasta} " \
          f"-out {output_file} -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue' " \
          f"-perc_identity 80 -word_size 7"
    
    subprocess.run(cmd, shell=True)
    return output_file

def parse_blast_results(blast_output):
    """Parse BLAST results, return circRNA IDs that meet criteria"""
    ires_circs = set()
    
    with open(blast_output, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            circ_id = fields[0]
            identity = float(fields[2])
            length = int(fields[3])
            
            if identity >= 80 and length >= 30:
                ires_circs.add(circ_id)
    
    return ires_circs

def analyze_circ_types(circ_ids, tsv_file):
    """Analyze types of circRNAs with IRES"""
    # Read TSV file
    df = pd.read_csv(tsv_file, sep='\t')
    
    # Initialize counters
    type_counts = {
        'with_ires': defaultdict(int),
        'without_ires': defaultdict(int)
    }
    
    # Count each type
    for _, row in df.iterrows():
        circ_type = row['CircType']
        if row['CircID'] in circ_ids:
            type_counts['with_ires'][circ_type] += 1
        else:
            type_counts['without_ires'][circ_type] += 1
    
    return type_counts

def main():
    # File paths
    query_fasta = "/data3/a549_circle_sequences.fasta"
    ires_fasta = "/data3/H23/ires_A549/Human_IRES.fa"
    tsv_file = "/data3/H23/ribo-seq-a549/A549_circleRNA_predict_result.tsv"
    blast_output = "ires_blast_results.txt"
    
    # Run BLASTN
    print("Running BLASTN alignment...")
    run_blastn(query_fasta, ires_fasta, blast_output)
    
    # Parse BLAST results
    print("Parsing BLAST results...")
    ires_circs = parse_blast_results(blast_output)
    
    # Analyze circRNA types
    print("Analyzing circRNA types...")
    type_counts = analyze_circ_types(ires_circs, tsv_file)
    
    # Output results
    print("\n=== Analysis Results ===")
    print("\nCircRNA types with IRES:")
    for circ_type, count in type_counts['with_ires'].items():
        print(f"{circ_type}: {count}")
    
    print("\nCircRNA types without IRES:")
    for circ_type, count in type_counts['without_ires'].items():
        print(f"{circ_type}: {count}")
    
    # Save information of circRNAs with IRES to file
    df = pd.read_csv(tsv_file, sep='\t')
    ires_df = df[df['CircID'].isin(ires_circs)]
    ires_df.to_csv('circrnas_with_ires.tsv', sep='\t', index=False)
    print(f"\nInformation of circRNAs with IRES saved to: circrnas_with_ires.tsv")

if __name__ == "__main__":
    main()
