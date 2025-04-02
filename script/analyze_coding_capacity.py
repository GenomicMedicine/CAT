import pandas as pd
import re
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

def process_circid(circ_id):
    """Process CircID to get host gene, chr, fusion junctions"""
    parts = circ_id.split('/')
    host_gene = parts[0]
    loc_parts = parts[1].split(':')
    chr_name = loc_parts[0]
    positions = loc_parts[1].split('-')
    junction1 = int(positions[0])
    junction2 = int(positions[1])
    
    # Ensure junction1 is smaller than junction2
    if junction1 > junction2:
        junction1, junction2 = junction2, junction1
        
    return host_gene, chr_name, junction1, junction2

def process_circtype(circ_type):
    """Process CircType, standardize naming"""
    if 'antisense' in circ_type.lower() or 'interior' in circ_type.lower():
        return 'complete interior'
    elif 'partial' in circ_type.lower():
        return 'half interior'
    elif 'canonical' in circ_type.lower():
        return 'canonical'
    return circ_type

def create_output_df(df):
    """Create output DataFrame"""
    results = []
    for _, row in df.iterrows():
        circ_type = process_circtype(row['CircType'])
        host_gene, chr_name, junction1, junction2 = process_circid(row['CircID'])
        
        # If host gene is intergenic, set to empty value
        if host_gene.lower() == 'intergenic':
            host_gene = ''
            
        results.append({
            'type': circ_type,
            'host gene': host_gene,
            'chr': chr_name,
            'fusion junction1': junction1,
            'fusion junction2': junction2,
            'strand': row['Strand(circ/gene)'][0]  # Take first character as strand
        })
    return pd.DataFrame(results)

def count_canonical_noncanonical(df):
    """Count canonical and noncanonical numbers"""
    canonical = 0
    noncanonical = 0
    for _, row in df.iterrows():
        if 'canonical' in row['CircType'].lower():
            canonical += 1
        elif any(x in row['CircType'].lower() for x in ['intergenic', 'interior', 'partial', 'antisense']):
            noncanonical += 1
    return canonical, noncanonical

def plot_venn_diagram(ribo_ids, ires_ids, orf_ids, output_file='venn_diagram.png'):
    """Draw Venn diagram for three analysis methods"""
    plt.figure(figsize=(10, 10))
    venn3([set(ribo_ids), set(ires_ids), set(orf_ids)], 
          set_labels=('Ribo-seq', 'IRES', 'ORF'))
    plt.title('Intersection of Coding Capacity Analysis Methods')
    plt.savefig(output_file)
    plt.close()

def analyze_intersections(ribo_ids, ires_ids, orf_ids):
    """Analyze intersections of three methods"""
    # Calculate all possible intersections
    ribo_ires = ribo_ids & ires_ids
    ribo_orf = ribo_ids & orf_ids
    ires_orf = ires_ids & orf_ids
    all_three = ribo_ids & ires_ids & orf_ids
    
    # Print intersection information
    print("\n=== Intersection Analysis ===")
    print(f"Ribo-seq & IRES: {len(ribo_ires)} circRNAs")
    print(f"Ribo-seq & ORF: {len(ribo_orf)} circRNAs")
    print(f"IRES & ORF: {len(ires_orf)} circRNAs")
    print(f"All three methods: {len(all_three)} circRNAs")
    
    # Save intersection results
    with open('intersection_results.txt', 'w') as f:
        f.write("=== Ribo-seq & IRES ===\n")
        f.write("\n".join(ribo_ires) + "\n\n")
        
        f.write("=== Ribo-seq & ORF ===\n")
        f.write("\n".join(ribo_orf) + "\n\n")
        
        f.write("=== IRES & ORF ===\n")
        f.write("\n".join(ires_orf) + "\n\n")
        
        f.write("=== All three methods ===\n")
        f.write("\n".join(all_three))
    
    return all_three

def validate_orf_position(circ_id, start_pos, end_pos):
    """
    Validate if ORF position matches circRNA length
    
    Args:
        circ_id: circRNA ID
        start_pos: ORF start position
        end_pos: ORF end position
    
    Returns:
        bool: True if ORF position is valid, False otherwise
    """
    try:
        # Extract junction positions from circRNA ID
        positions = circ_id.split(':')[1].split('-')
        junction1 = int(positions[0])
        junction2 = int(positions[2])  # Use third number as second junction
        
        # Calculate circRNA length
        circ_length = abs(junction2 - junction1)
        
        # Calculate ORF length
        orf_length = end_pos - start_pos + 1
        
        # Check if 1x, 2x or 3x length is within ORF range
        for multiplier in [1, 2, 3]:
            length = circ_length * multiplier
            if start_pos <= length <= end_pos:
                return True
            
        return False
    except:
        return False

def process_orf_file(orf_file):
    """Process ORF file, keep only valid ORFs"""
    valid_orf_circs = set()
    
    with open(orf_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 5:
                circ_id = fields[0]
                start_pos = int(fields[3])
                end_pos = int(fields[4])
                
                if validate_orf_position(circ_id, start_pos, end_pos):
                    valid_orf_circs.add(circ_id)
    
    return valid_orf_circs

def main():
    # Read original circRNA data
    all_circs = pd.read_csv('/data3/H23/ribo-seq-a549/A549_circleRNA_predict_result.tsv', sep='\t')
    
    # 1. Count all circRNA types
    print("\n=== All circRNA Type Statistics ===")
    canonical_count, noncanonical_count = count_canonical_noncanonical(all_circs)
    print(f"Canonical circRNAs: {canonical_count}")
    print(f"Non-canonical circRNAs: {noncanonical_count}")
    
    # 2. Process Ribo-seq results
    print("\n=== Ribo-seq Result Analysis ===")
    # Read Ribo-seq result file
    ribo_df = pd.read_csv('/data3/A549_riboseq/perfect_matches_with_type.txt', sep='\t')
    ribo_circs = all_circs[all_circs['CircID'].isin(set(ribo_df['CircID']))]
    canonical_ribo, noncanonical_ribo = count_canonical_noncanonical(ribo_circs)
    print(f"Ribo-seq supported canonical circRNAs: {canonical_ribo}")
    print(f"Ribo-seq supported non-canonical circRNAs: {noncanonical_ribo}")
    ribo_output = create_output_df(ribo_circs)
    ribo_output.to_csv('ribo_seq_circs.csv', index=False)
    ribo_ids = set(ribo_df['CircID'])
    
    # 3. Process IRES results
    print("\n=== IRES Result Analysis ===")
    ires_df = pd.read_csv('/data3/A549_coding_capacity/circrnas_with_ires.tsv', sep='\t')
    canonical_ires, noncanonical_ires = count_canonical_noncanonical(ires_df)
    print(f"IRES-containing canonical circRNAs: {canonical_ires}")
    print(f"IRES-containing non-canonical circRNAs: {noncanonical_ires}")
    ires_output = create_output_df(ires_df)
    ires_output.to_csv('ires_circs.csv', index=False)
    ires_ids = set(ires_df['CircID'])
    
    # 4. Process ORF results
    print("\n=== ORF Result Analysis ===")
    # First get valid ORF circRNAs
    valid_orf_circ_ids = process_orf_file('/data3/A549orf.gtf')
    orf_circs = all_circs[all_circs['CircID'].isin(valid_orf_circ_ids)]
    canonical_orf, noncanonical_orf = count_canonical_noncanonical(orf_circs)
    print(f"Valid ORF-containing canonical circRNAs: {canonical_orf}")
    print(f"Valid ORF-containing non-canonical circRNAs: {noncanonical_orf}")
    orf_output = create_output_df(orf_circs)
    orf_output.to_csv('valid_orf_circs.csv', index=False)
    
    # 5. Analyze intersections and draw Venn diagram
    common_circs = analyze_intersections(ribo_ids, ires_ids, valid_orf_circ_ids)
    plot_venn_diagram(ribo_ids, ires_ids, valid_orf_circ_ids)
    
    # 6. Analyze common circRNA types
    if common_circs:
        common_df = all_circs[all_circs['CircID'].isin(common_circs)]
        canonical_common, noncanonical_common = count_canonical_noncanonical(common_df)
        print("\n=== Statistics of circRNAs Supported by All Three Methods ===")
        print(f"Canonical circRNAs: {canonical_common}")
        print(f"Non-canonical circRNAs: {noncanonical_common}")
        
        # Save detailed information of common circRNAs
        common_output = create_output_df(common_df)
        common_output.to_csv('common_circs.csv', index=False)

if __name__ == "__main__":
    main() 