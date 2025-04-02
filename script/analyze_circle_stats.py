import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

def load_data(excel_path):
    """Load data from all sheets"""
    sheets = ['circle_canonical', 'circle_half_interior', 'circle_lariat', 
              'circle_intergenic', 'circle_complete_interior']
    data = {}
    for sheet in sheets:
        data[sheet] = pd.read_excel(excel_path, sheet_name=sheet)
    return data

def analyze_chromosomes(data, output_dir):
    """Analyze chromosome distribution"""
    # Sheets to analyze
    analyze_sheets = ['circle_canonical', 'circle_half_interior', 
                     'circle_intergenic', 'circle_complete_interior']
    
    # Create storage for results
    chr_stats = {}
    
    # Define chromosome order
    chr_order = [str(i) for i in range(1, 23)] + ['X', 'Y']
    
    for sheet_name in analyze_sheets:
        df = data[sheet_name]
        # Extract chromosome numbers
        chr_counts = df['chr'].value_counts()
        # Normalize chromosome names
        chr_counts.index = chr_counts.index.astype(str).str.replace('chr', '')
        
        # Keep only main chromosomes
        main_chr_counts = chr_counts[chr_counts.index.isin(chr_order)]
        chr_stats[sheet_name] = main_chr_counts
    
    # Create plot
    plt.figure(figsize=(15, 8))
    bar_width = 0.2
    x = np.arange(len(chr_order))
    
    for i, (sheet_name, counts) in enumerate(chr_stats.items()):
        values = [counts.get(chr_name, 0) for chr_name in chr_order]
        plt.bar(x + i*bar_width, values, bar_width, label=sheet_name)
    
    plt.xlabel('Chromosome')
    plt.ylabel('Count')
    plt.title('Chromosome Distribution')
    plt.xticks(x + bar_width*1.5, chr_order)
    plt.legend()
    plt.tight_layout()
    
    # Save plot and data
    plt.savefig(f'{output_dir}/chromosome_distribution.png')
    plt.close()
    
    # Save data to Excel
    pd.DataFrame(chr_stats).to_excel(f'{output_dir}/chromosome_stats.xlsx')

def analyze_circle_length(data, output_dir):
    """Analyze CircLen distribution"""
    length_stats = {}
    
    plt.figure(figsize=(15, 8))
    for sheet_name, df in data.items():
        # Calculate quantiles
        stats = df['CircLen'].describe()
        length_stats[sheet_name] = stats
        
        # Create density plot
        sns.kdeplot(data=df['CircLen'], label=sheet_name)
    
    plt.xlabel('Circle Length')
    plt.ylabel('Density')
    plt.title('CircRNA Length Distribution')
    plt.legend()
    plt.tight_layout()
    
    # Save plot and data
    plt.savefig(f'{output_dir}/length_distribution.png')
    plt.close()
    
    # Save data to Excel
    pd.DataFrame(length_stats).to_excel(f'{output_dir}/length_stats.xlsx')

def analyze_read_counts(data, output_dir):
    """Analyze ReadCount distribution"""
    thresholds = [5, 10, 20, 50, 100, 200, 500, 1000]
    count_stats = {}
    
    for sheet_name, df in data.items():
        counts = []
        for threshold in thresholds:
            count = len(df[df['ReadCount'] >= threshold])
            counts.append(count)
        count_stats[sheet_name] = counts
    
    # Create plot
    plt.figure(figsize=(15, 8))
    x = np.arange(len(thresholds))
    for sheet_name, counts in count_stats.items():
        plt.plot(x, counts, marker='o', label=sheet_name)
    
    plt.xlabel('Read Count Threshold')
    plt.ylabel('Number of circRNAs')
    plt.title('circRNA Read Count Distribution')
    plt.xticks(x, thresholds)
    plt.legend()
    plt.tight_layout()
    
    # Save plot and data
    plt.savefig(f'{output_dir}/read_count_distribution.png')
    plt.close()
    
    # Save data to Excel
    pd.DataFrame(count_stats, index=thresholds).to_excel(f'{output_dir}/read_count_stats.xlsx')

def analyze_gene_types(data, output_dir):
    """Analyze GeneType distribution"""
    # Sheets to analyze
    analyze_sheets = ['circle_canonical', 'circle_half_interior', 'circle_lariat', 
                     'circle_complete_interior']
    
    gene_type_stats = {}
    
    for sheet_name in analyze_sheets:
        df = data[sheet_name]
        gene_counts = df['GeneType'].value_counts()
        gene_type_stats[sheet_name] = gene_counts
    
    # Create percentage plot
    plt.figure(figsize=(15, 8))
    df_percentages = pd.DataFrame(gene_type_stats).fillna(0)
    df_percentages = df_percentages.div(df_percentages.sum()) * 100
    df_percentages.plot(kind='bar', stacked=True)
    
    plt.xlabel('Gene Type')
    plt.ylabel('Percentage')
    plt.title('Gene Type Distribution (Percentage)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/gene_type_percentage.png')
    plt.close()
    
    # Create value line plot
    plt.figure(figsize=(15, 8))
    df_absolute = pd.DataFrame(gene_type_stats).fillna(0)
    for column in df_absolute.columns:
        plt.plot(df_absolute.index, df_absolute[column], marker='o', label=column)
    
    plt.xlabel('Gene Type')
    plt.ylabel('Count')
    plt.title('Gene Type Distribution')
    plt.xticks(rotation=45)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{output_dir}/gene_type.png')
    plt.close()
    
    # Save data to Excel
    writer = pd.ExcelWriter(f'{output_dir}/gene_type_stats.xlsx')
    df_percentages.to_excel(writer, sheet_name='Percentages')
    df_absolute.to_excel(writer, sheet_name='Absolute')
    writer.close()

def analyze_canonical_vs_noncanonical(data, output_dir):
    """Compare canonical and non-canonical circRNAs"""
    # Merge data
    canonical = data['circle_canonical'].copy()
    canonical['Type'] = 'Canonical'
    
    non_canonical = pd.concat([
        data['circle_half_interior'],
        data['circle_complete_interior'],
        data['circle_intergenic']
    ])
    non_canonical['Type'] = 'Non-canonical'
    
    combined_data = pd.concat([canonical, non_canonical])
    
    # 1. Analyze chromosome distribution
    chr_order = [str(i) for i in range(1, 23)] + ['X', 'Y']
    chr_stats = {}
    
    for type_name, group in combined_data.groupby('Type'):
        chr_counts = group['chr'].value_counts()
        chr_counts.index = chr_counts.index.astype(str).str.replace('chr', '')
        main_chr_counts = chr_counts[chr_counts.index.isin(chr_order)]
        chr_stats[type_name] = main_chr_counts
    
    plt.figure(figsize=(15, 8))
    bar_width = 0.35
    x = np.arange(len(chr_order))
    
    for i, (type_name, counts) in enumerate(chr_stats.items()):
        values = [counts.get(chr_name, 0) for chr_name in chr_order]
        plt.bar(x + i*bar_width, values, bar_width, label=type_name)
    
    plt.xlabel('Chromosome')
    plt.ylabel('Count')
    plt.title('Chromosome Distribution (Canonical vs Non-canonical)')
    plt.xticks(x + bar_width/2, chr_order)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{output_dir}/chromosome_distribution_canonical_vs_noncanonical.png')
    plt.close()
    
    pd.DataFrame(chr_stats).to_excel(f'{output_dir}/chromosome_stats_canonical_vs_noncanonical.xlsx')
    
    # 2. Analyze CircLen distribution
    plt.figure(figsize=(15, 8))
    for type_name, group in combined_data.groupby('Type'):
        sns.kdeplot(data=group['CircLen'], label=type_name)
    
    plt.xlabel('Circle Length')
    plt.ylabel('Density')
    plt.title('CircRNA Length Distribution (Canonical vs Non-canonical)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{output_dir}/length_distribution_canonical_vs_noncanonical.png')
    plt.close()
    
    length_stats = combined_data.groupby('Type')['CircLen'].describe()
    length_stats.to_excel(f'{output_dir}/length_stats_canonical_vs_noncanonical.xlsx')
    
    # 3. Analyze ReadCount distribution
    thresholds = [5, 10, 20, 50, 100, 200, 500, 1000]
    count_stats = {}
    
    for type_name, group in combined_data.groupby('Type'):
        counts = []
        for threshold in thresholds:
            count = len(group[group['ReadCount'] >= threshold])
            counts.append(count)
        count_stats[type_name] = counts
    
    plt.figure(figsize=(15, 8))
    x = np.arange(len(thresholds))
    for type_name, counts in count_stats.items():
        plt.plot(x, counts, marker='o', label=type_name)
    
    plt.xlabel('Read Count Threshold')
    plt.ylabel('Number of circRNAs')                                                                                                                                             
    plt.title('circRNA Read Count Distribution (Canonical vs Non-canonical)')
    plt.xticks(x, thresholds)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{output_dir}/read_count_distribution_canonical_vs_noncanonical.png')
    plt.close()
    
    pd.DataFrame(count_stats, index=thresholds).to_excel(
        f'{output_dir}/read_count_stats_canonical_vs_noncanonical.xlsx'
    )

def export_length_comparison(data, output_dir):
    """Export length comparison table for canonical and non-canonical"""
    # Prepare data
    canonical_lengths = data['circle_canonical']['CircLen'].tolist()
    
    # Merge non-canonical data
    non_canonical_lengths = pd.concat([
        data['circle_half_interior']['CircLen'],
        data['circle_complete_interior']['CircLen'],
        data['circle_intergenic']['CircLen']
    ]).tolist()
    
    # Create DataFrame
    # Since the two datasets may have different lengths, create list with max length
    max_length = max(len(canonical_lengths), len(non_canonical_lengths))
    length_df = pd.DataFrame({
        'canonical': canonical_lengths + [None] * (max_length - len(canonical_lengths)),
        'non_canonical': non_canonical_lengths + [None] * (max_length - len(non_canonical_lengths))
    })
    
    # Save as CSV
    length_df.to_csv(f'{output_dir}/circle_lengths.csv', index=False)

def main():
    # Set input and output paths
    excel_path = "/data6/H23_public466/466_newcat_test/Pub_h23_classify.xlsx"
    excel_path = "/data6/H23_public466/466_newcat_test/circleRNA_predict_result_paper_supplement_cutoff5_1_15.xlsx" 
    output_dir = "analysis_results"
    Path(output_dir).mkdir(exist_ok=True)
    
    # Load data
    data = load_data(excel_path)
    
    # Perform analyses
    analyze_chromosomes(data, output_dir)
    analyze_circle_length(data, output_dir)
    analyze_read_counts(data, output_dir)
    analyze_gene_types(data, output_dir)
    analyze_canonical_vs_noncanonical(data, output_dir)
    export_length_comparison(data, output_dir)

if __name__ == "__main__":
    main() 