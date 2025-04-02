import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import pairwise2
from Bio.Seq import Seq
import numpy as np
from collections import defaultdict
from multiprocessing import Pool
from itertools import product
import time
from functools import partial

def process_file(args):
    """Process a single circRNA file"""
    root, file = args
    if file == "circleRNA_predict_result.tsv":
        file_path = os.path.join(root, file)
        
        # Get full path
        abs_path = os.path.abspath(file_path)
        path_parts = abs_path.split('/')
        
        # Find directory containing tissue info (third from last)
        tissue_dir = path_parts[-3]
        
        # Split directory name and extract tissue part
        parts = tissue_dir.split('_')
        if len(parts) > 2:
            # Remove first two elements (ID and Human/Mouse), merge remaining
            tissue = ''.join(parts[2:])
            # Remove possible spaces
            tissue = tissue.replace(' ', '')
            if tissue.lower() != 'lung':  # Only process lung tissue
                return None
        else:
            return None  # Skip if tissue cannot be determined
            
        df = pd.read_csv(file_path, sep='\t')
        df['tissue'] = tissue
        df = df[df['LenShs'] <= 20]
        return df[['CircID', 'Reference', 'ReadCount', 'tissue', 'CircType']]
    return None

def load_circ_data(base_path, species):
    """Load all circRNA data for specified species"""
    print(f"Loading {species} data...")
    
    res_path = os.path.join(base_path, f"{species}_res")
    
    # Collect all file paths
    file_args = [
        (root, file)
        for root, _, files in os.walk(res_path)
        for file in files
        if file == "circleRNA_predict_result.tsv"
    ]
    
    # Print found file paths and extracted tissues (for debugging)
    print(f"Found {len(file_args)} files:")
    for root, file in file_args:
        full_path = os.path.join(root, file)
        tissue_dir = full_path.split('/')[-3]
        tissue = ''.join(tissue_dir.split('_')[2:]).replace(' ', '')
        print(f"Path: {full_path}")
        print(f"Extracted tissue: {tissue}")
        print("---")
    
    # Process files using multiprocessing
    with Pool() as pool:
        all_data = pool.map(process_file, file_args)
    
    # Filter None and merge data
    all_data = [df for df in all_data if df is not None]
    result = pd.concat(all_data, ignore_index=True)
    
    # All unique tissues found (for validation)
    unique_tissues = result['tissue'].unique()
    print(f"\nFound tissues: {sorted(unique_tissues)}")
    
    print(f"{species} data loading complete, {len(result)} records total")
    return result

def analyze_conservation_within_species(circ_data, output_dir, species):
    """Analyze circRNA conservation between tissues within species and save detailed data"""
    print(f"Analyzing {species} tissue conservation...")
    
    # Get tissue-specific circRNA sets
    tissue_circs = {
        tissue: set(group['CircID']) 
        for tissue, group in circ_data.groupby('tissue')
    }
    
    # Save tissue-specific circRNA lists
    for tissue, circs in tissue_circs.items():
        pd.DataFrame(list(circs), columns=['CircID']).to_csv(
            os.path.join(output_dir, f'{species}_{tissue}_specific_circrnas.csv'),
            index=False
        )
    
    tissues = list(tissue_circs.keys())
    overlap_matrix = np.zeros((len(tissues), len(tissues)))
    overlap_details = []
    
    # Calculate conservation matrix and save overlap details
    for i, t1 in enumerate(tissues):
        for j in range(i, len(tissues)):
            t2 = tissues[j]
            overlap_circs = tissue_circs[t1] & tissue_circs[t2]
            overlap = len(overlap_circs)
            overlap_matrix[i,j] = overlap
            overlap_matrix[j,i] = overlap
            
            # Save overlap details
            overlap_details.append({
                'tissue1': t1,
                'tissue2': t2,
                'overlap_count': overlap,
                'tissue1_total': len(tissue_circs[t1]),
                'tissue2_total': len(tissue_circs[t2])
            })
    
    # Save overlap details
    pd.DataFrame(overlap_details).to_csv(
        os.path.join(output_dir, f'{species}_conservation_details.csv'),
        index=False
    )
    
    return pd.DataFrame(overlap_matrix, index=tissues, columns=tissues)

def analyze_expression_patterns(circ_data, output_dir, species):
    """Analyze circRNA expression patterns and save detailed data"""
    # Basic statistics
    expr_summary = circ_data.groupby('tissue')['ReadCount'].agg(['mean', 'median', 'std'])
    
    # Detailed expression analysis
    detailed_stats = []
    for tissue, group in circ_data.groupby('tissue'):
        stats = {
            'tissue': tissue,
            'total_circs': len(group),
            'mean_expression': group['ReadCount'].mean(),
            'median_expression': group['ReadCount'].median(),
            'std_expression': group['ReadCount'].std(),
            'min_expression': group['ReadCount'].min(),
            'max_expression': group['ReadCount'].max(),
            'q1_expression': group['ReadCount'].quantile(0.25),
            'q3_expression': group['ReadCount'].quantile(0.75)
        }
        detailed_stats.append(stats)
    
    # Save detailed statistics
    pd.DataFrame(detailed_stats).to_csv(
        os.path.join(output_dir, f'{species}_expression_detailed_stats.csv'),
        index=False
    )
    
    # Save raw expression data by tissue
    for tissue, group in circ_data.groupby('tissue'):
        group[['CircID', 'ReadCount']].to_csv(
            os.path.join(output_dir, f'{species}_{tissue}_expression.csv'),
            index=False
        )
    
    return expr_summary

def calculate_sequence_similarity(seq1, seq2):
    """Calculate sequence similarity between two sequences"""
    try:
        alignment = pairwise2.align.globalxx(seq1, seq2, score_only=True)
        return alignment / max(len(seq1), len(seq2))
    except Exception:
        return 0

def process_pair(pair, species1_seqs, species2_seqs, similarity_threshold):
    """Process a pair of circRNAs from two species, calculate similarity"""
    species1_id, species2_id = pair
    species1_seq = species1_seqs[species1_id]
    species2_seq = species2_seqs[species2_id]
    
    similarity = calculate_sequence_similarity(species1_seq, species2_seq)
    if similarity >= similarity_threshold:
        return {
            'species1_circ': species1_id,
            'species2_circ': species2_id,
            'similarity': similarity
        }
    return None

def find_homologous_circs(species1_data, species2_data, species1_name, species2_name, similarity_threshold=0.8):
    """Find homologous circRNAs between two species"""
    print(f"Starting homology search between {species1_name} and {species2_name}...")
    start_time = time.time()
    
    # Pre-calculate sequence objects
    species1_seqs = {row['CircID']: Seq(row['Reference'][125:175]) 
                    for _, row in species1_data.iterrows()}
    species2_seqs = {row['CircID']: Seq(row['Reference'][125:175]) 
                    for _, row in species2_data.iterrows()}
    
    # Create all possible pairs
    pairs = list(product(species1_seqs.keys(), species2_seqs.keys()))
    
    # Use multiprocessing
    with Pool() as pool:
        process_func = partial(process_pair, 
                             species1_seqs=species1_seqs, 
                             species2_seqs=species2_seqs, 
                             similarity_threshold=similarity_threshold)
        chunk_size = max(1, len(pairs) // (os.cpu_count() * 4))
        results = pool.map(process_func, pairs, chunksize=chunk_size)
    
    homologs = [r for r in results if r is not None]
    
    print(f"Homology analysis complete, took {time.time() - start_time:.2f} seconds")
    print(f"Found {len(homologs)} homologous circRNA pairs")
    
    return pd.DataFrame(homologs)

def categorize_circ_type(circ_type):
    """Categorize circRNA type into five main categories"""
    if 'canonical' in circ_type.lower():
        return 'Canonical'
    elif 'partial' in circ_type.lower():
        return 'Half Interior'
    elif 'lariat' in circ_type.lower():
        return 'Lariat'
    elif 'interior' in circ_type.lower() or 'antisense' in circ_type.lower():
        return 'Complete Interior'
    elif 'intergenic' in circ_type.lower():
        return 'Intergenic'

def analyze_type_distribution(data, species, output_dir):
    """Analyze circRNA type distribution"""
    # Add category column
    data['Category'] = data['CircType'].apply(categorize_circ_type)
    
    # Calculate type counts and percentages
    type_counts = data['Category'].value_counts()
    type_percentages = (type_counts / len(data) * 100).round(2)
    
    # Save distribution data
    distribution_df = pd.DataFrame({
        'Count': type_counts,
        'Percentage': type_percentages
    })
    distribution_df.to_csv(os.path.join(output_dir, f'{species}_type_distribution.csv'))
    
    # Plot pie chart
    plt.figure(figsize=(10, 8))
    plt.pie(type_percentages, labels=type_percentages.index, 
            autopct='%1.1f%%', startangle=90)
    plt.title(f'{species.capitalize()} Lung circRNA Type Distribution')
    plt.axis('equal')
    plt.savefig(os.path.join(output_dir, f'{species}_type_distribution_pie.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    return distribution_df

def find_homologous_circs_by_type(species1_data, species2_data, species1_name, species2_name):
    """Find homologous circRNAs based on type"""
    # Add category column
    species1_data['Category'] = species1_data['CircType'].apply(categorize_circ_type)
    species2_data['Category'] = species2_data['CircType'].apply(categorize_circ_type)
    
    # Separate canonical and noncanonical
    species1_canonical = species1_data[species1_data['Category'] == 'Canonical']
    species2_canonical = species2_data[species2_data['Category'] == 'Canonical']
    
    species1_noncanonical = species1_data[species1_data['Category'].isin(['Half Interior', 'Complete Interior', 'Intergenic'])]
    species2_noncanonical = species2_data[species2_data['Category'].isin(['Half Interior', 'Complete Interior', 'Intergenic'])]
    
    # Find homology
    canonical_homologs = find_homologous_circs(species1_canonical, species2_canonical, 
                                              species1_name, species2_name)
    noncanonical_homologs = find_homologous_circs(species1_noncanonical, species2_noncanonical, 
                                                 species1_name, species2_name)
    
    return canonical_homologs, noncanonical_homologs

def plot_homology_counts(homologs_df, title, output_path):
    """Plot homology counts for different similarity thresholds"""
    similarity_thresholds = np.arange(0.8, 1.00, 0.05)
    counts = []
    
    for threshold in similarity_thresholds:
        count = len(homologs_df[homologs_df['similarity'] >= threshold])
        counts.append(count)
    
    plt.figure(figsize=(10, 6))
    plt.plot(similarity_thresholds, counts, marker='o')
    plt.xlabel('Similarity Threshold')
    plt.ylabel('Number of Homologous Pairs')
    plt.title(title)
    plt.grid(True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def save_combined_homology_stats(human_macaque_can, human_macaque_noncan,
                               macaque_mouse_can, macaque_mouse_noncan,
                               human_mouse_can, human_mouse_noncan,
                               output_dir):
    """Save all homology analysis results to a CSV file"""
    # Define similarity thresholds
    similarity_thresholds = np.arange(0.8, 1.01, 0.05)
    
    # Create results dictionary
    results = {
        'Similarity_Threshold': similarity_thresholds,
        'Human_Macaque_Canonical': [],
        'Human_Macaque_Noncanonical': [],
        'Macaque_Mouse_Canonical': [],
        'Macaque_Mouse_Noncanonical': [],
        'Human_Mouse_Canonical': [],
        'Human_Mouse_Noncanonical': []
    }
    
    # Calculate homology counts for each threshold
    for threshold in similarity_thresholds:
        # Human-Macaque
        hm_can_count = len(human_macaque_can[human_macaque_can['similarity'] >= threshold]) if human_macaque_can is not None else 0
        hm_noncan_count = len(human_macaque_noncan[human_macaque_noncan['similarity'] >= threshold]) if human_macaque_noncan is not None else 0
        
        # Macaque-Mouse
        mm_can_count = len(macaque_mouse_can[macaque_mouse_can['similarity'] >= threshold]) if macaque_mouse_can is not None else 0
        mm_noncan_count = len(macaque_mouse_noncan[macaque_mouse_noncan['similarity'] >= threshold]) if macaque_mouse_noncan is not None else 0
        
        # Human-Mouse
        h_m_can_count = len(human_mouse_can[human_mouse_can['similarity'] >= threshold]) if human_mouse_can is not None else 0
        h_m_noncan_count = len(human_mouse_noncan[human_mouse_noncan['similarity'] >= threshold]) if human_mouse_noncan is not None else 0
        
        results['Human_Macaque_Canonical'].append(hm_can_count)
        results['Human_Macaque_Noncanonical'].append(hm_noncan_count)
        results['Macaque_Mouse_Canonical'].append(mm_can_count)
        results['Macaque_Mouse_Noncanonical'].append(mm_noncan_count)
        results['Human_Mouse_Canonical'].append(h_m_can_count)
        results['Human_Mouse_Noncanonical'].append(h_m_noncan_count)
    
    # Create DataFrame and save
    df = pd.DataFrame(results)
    df.to_csv(os.path.join(output_dir, 'combined_homology_statistics.csv'), index=False)
    print(f"Combined homology statistics saved to {output_dir}/combined_homology_statistics.csv")

def main():
    start_time = time.time()
    base_path = "/data5/lk/result/"
    output_dir = "lung_results"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    human_data = load_circ_data(base_path, "human")
    mouse_data = load_circ_data(base_path, "mouse")
    macaque_data = load_circ_data(base_path, "macaque")
    
    # Analyze type distribution
    human_dist = analyze_type_distribution(human_data, 'human', output_dir)
    mouse_dist = analyze_type_distribution(mouse_data, 'mouse', output_dir)
    macaque_dist = analyze_type_distribution(macaque_data, 'macaque', output_dir)
    
    # Analyze homology between species
    human_macaque_can, human_macaque_noncan = find_homologous_circs_by_type(
        human_data, macaque_data, "human", "macaque")
    macaque_mouse_can, macaque_mouse_noncan = find_homologous_circs_by_type(
        macaque_data, mouse_data, "macaque", "mouse")
    human_mouse_can, human_mouse_noncan = find_homologous_circs_by_type(
        human_data, mouse_data, "human", "mouse")
    
    # Save homology analysis results
    human_macaque_can.to_csv(os.path.join(output_dir, 'human_macaque_canonical_homologs.csv'), index=False)
    human_macaque_noncan.to_csv(os.path.join(output_dir, 'human_macaque_noncanonical_homologs.csv'), index=False)
    macaque_mouse_can.to_csv(os.path.join(output_dir, 'macaque_mouse_canonical_homologs.csv'), index=False)
    macaque_mouse_noncan.to_csv(os.path.join(output_dir, 'macaque_mouse_noncanonical_homologs.csv'), index=False)
    human_mouse_can.to_csv(os.path.join(output_dir, 'human_mouse_canonical_homologs.csv'), index=False)
    human_mouse_noncan.to_csv(os.path.join(output_dir, 'human_mouse_noncanonical_homologs.csv'), index=False)
    
    # Save combined statistics
    save_combined_homology_stats(
        human_macaque_can, human_macaque_noncan,
        macaque_mouse_can, macaque_mouse_noncan,
        human_mouse_can, human_mouse_noncan,
        output_dir
    )
    
    # Plot homology counts
    plot_homology_counts(human_macaque_can, 
                        'Human-Macaque Canonical Homologs',
                        os.path.join(output_dir, 'human_macaque_canonical_counts.png'))
    plot_homology_counts(human_macaque_noncan,
                        'Human-Macaque Non-canonical Homologs',
                        os.path.join(output_dir, 'human_macaque_noncanonical_counts.png'))
    plot_homology_counts(macaque_mouse_can,
                        'Macaque-Mouse Canonical Homologs',
                        os.path.join(output_dir, 'macaque_mouse_canonical_counts.png'))
    plot_homology_counts(macaque_mouse_noncan,
                        'Macaque-Mouse Non-canonical Homologs',
                        os.path.join(output_dir, 'macaque_mouse_noncanonical_counts.png'))
    plot_homology_counts(human_mouse_can,
                        'Human-Mouse Canonical Homologs',
                        os.path.join(output_dir, 'human_mouse_canonical_counts.png'))
    plot_homology_counts(human_mouse_noncan,
                        'Human-Mouse Non-canonical Homologs',
                        os.path.join(output_dir, 'human_mouse_noncanonical_counts.png'))
    
    print(f"Analysis completed! Total time: {(time.time() - start_time) / 60:.2f} minutes")
    print(f"Results saved in {output_dir}/")

if __name__ == "__main__":
    main()