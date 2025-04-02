import pandas as pd
from Bio import SeqIO
import pysam
from Bio.Seq import Seq

def get_sequence(chr, start, end, strand, ref_file="/home/lab/lk/code/bowtie2refer/Homo_sapiens.GRCh38.dna.primary_assembly.fa"):
    """Get sequence for specified region, regardless of strand"""
    try:
        with pysam.FastaFile(ref_file) as fasta:
            seq = fasta.fetch(str(chr), start, end)
            return seq
    except Exception as e:
        print(f"Error fetching sequence for {chr}: {e}")
        return None

def get_motif_pair(area1_start, area1_end, area2_start, area2_end, chr, strand):
    """Get all possible motif pairs"""
    if strand == '+':
        area2_seq = get_sequence(chr, area1_start, area1_end+2, strand)
        area1_seq = get_sequence(chr, area2_start-2, area2_end, strand)
    else:
        area1_seq = get_sequence(chr, area1_start-2, area1_end, strand)
        area2_seq = get_sequence(chr, area2_start, area2_end+2, strand)
    
    if not area1_seq or not area2_seq:
        return None, []
    
    motif_pairs = []
    max_length = min(len(area1_seq)-1, len(area2_seq)-1)
    for i in range(max_length):
        pair = (area1_seq[i:i+2].upper(), area2_seq[i:i+2].upper())
        if strand == '-':
            pair = (str(Seq(pair[1]).reverse_complement()), 
                   str(Seq(pair[0]).reverse_complement()))
        motif_pairs.append(f"{pair[0]}-{pair[1]}")
    
    # Return AG-GT if found, along with all motif pairs
    if "AG-GT" in motif_pairs:
        return "AG-GT", motif_pairs
    # Otherwise return first motif pair found, along with all motif pairs
    elif motif_pairs:
        return motif_pairs[0], motif_pairs
    
    return None, []

def analyze_sheet_motifs(sheet_df):
    """Analyze motif distribution in a sheet"""
    motif_counts = {"AG-GT": 0}  # Initialize AG-GT count
    
    # Print results for first 10 entries
    print("\nAnalyzing first 10 circRNAs:")
    for i in range(10):
        row = sheet_df.iloc[i]
        print(f"\nProcessing circRNA #{i+1}:")
        print(f"Position: chr{row['chr']}:{row['left boundary of shs area1']}-{row['right boundary of shs area1']} "
              f"and {row['left boundary of shs area2']}-{row['right boundary of shs area2']}")
        print(f"Strand: {row['Strand(circ/gene)'].split('/')[0]}")
        
        try:
            strand = row['Strand(circ/gene)'].split('/')[0]
            motif, all_motifs = get_motif_pair(
                int(row['left boundary of shs area1']),
                int(row['right boundary of shs area1']),
                int(row['left boundary of shs area2']),
                int(row['right boundary of shs area2']),
                str(row['chr']),
                strand
            )
            if motif:
                print(f"Selected Motif: {motif}")
                print(f"All possible Motifs: {all_motifs}")
                motif_counts[motif] = motif_counts.get(motif, 0) + 1
            
        except Exception as e:
            print(f"Processing error: {e}")
    
    # Process remaining data
    total = len(sheet_df)
    processed = 0
    
    for idx, row in sheet_df.iterrows():
        if idx < 10:  # Skip first 10
            continue
        try:
            strand = row['Strand(circ/gene)'].split('/')[0]
            motif, _ = get_motif_pair(  # Only take first return value (motif)
                int(row['left boundary of shs area1']),
                int(row['right boundary of shs area1']),
                int(row['left boundary of shs area2']),
                int(row['right boundary of shs area2']),
                str(row['chr']),
                strand
            )
            if motif:  # motif is string, can be used as dictionary key
                motif_counts[motif] = motif_counts.get(motif, 0) + 1
            processed += 1
            
            if processed % 1000 == 0:
                print(f"Processed {processed}/{total} records")
                
        except Exception as e:
            print(f"Error processing row {idx}: {e}")
            continue
    
    # Calculate percentages
    total_count = sum(motif_counts.values())
    if total_count > 0:  # Avoid division by zero
        motif_percentages = {k: (v/total_count)*100 for k, v in motif_counts.items()}
    else:
        motif_percentages = {k: 0 for k in motif_counts}
    return motif_percentages

def main():
    excel_path = "/data6/H23_public466/466_newcat_test/circleRNA_predict_result_paper_supplement_cutoff5_1_15.xlsx"
    sheets = ['circle_canonical', 'circle_half_interior', 'circle_lariat', 
              'circle_intergenic', 'circle_complete_interior']
    
    print("\nMotif Distribution Analysis Results:")
    print("-" * 70)
    
    # Create dictionary to store motif data for all sheets
    all_motif_data = {}
    
    for sheet in sheets:
        try:
            df = pd.read_excel(excel_path, sheet_name=sheet)
            df.name = sheet
            print(f"\nProcessing {sheet}...")
            motif_percentages = analyze_sheet_motifs(df)
            all_motif_data[sheet] = motif_percentages
            
            # Output percentages in fixed order
            print(f"\n{sheet} motif distribution:")
            temlist = []
            for motif in ["AG-GT", "AG-GC", "AC-CT", "GC-CT", "AC-AT", "AT-GT"]:
                print(f"{motif}: {motif_percentages.get(motif, 0):.3f}")
                temlist.append(motif_percentages.get(motif, 0))
            temlist.append(100 - sum(temlist))
            print(f"Other: {temlist[-1]:.3f}")
            
        except Exception as e:
            print(f"Error processing sheet {sheet}: {e}")
    
    # Convert all motif data to DataFrame and save to CSV
    all_data = []
    for sheet, motifs in all_motif_data.items():
        for motif, percentage in motifs.items():
            all_data.append({
                'Sheet': sheet,
                'Motif': motif,
                'Percentage': percentage
            })
    
    motif_df = pd.DataFrame(all_data)
    # Sort by Sheet and percentage descending
    motif_df = motif_df.sort_values(['Sheet', 'Percentage'], ascending=[True, False])
    # Save to CSV file
    output_path = "motif_distribution.csv"
    motif_df.to_csv(output_path, index=False)
    print(f"\nAll motif distribution data saved to: {output_path}")
    
    print("-" * 70)

if __name__ == "__main__":
    main()