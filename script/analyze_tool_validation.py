import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def get_validated_circles():
    """Read validated circRNAs"""
    # Read validated circRNA data
    df1 = pd.read_csv("/data6/circleRNA/circRNA_benchmarking/data/Supplementary_Table_3_selected_circRNAs.txt", sep="\t")
    
    # Filter circRNAs that passed all validations
    df1 = df1[(df1["amp_seq_val"]=="pass") & 
              (df1["amp_seq_val_detail"]=="pass") & 
              (df1["compound_val"]=="pass") & 
              (df1["RR_val"]=="pass") & 
              (df1["RR_val_detail"]=="pass") & 
              (df1["qPCR_val"]=="pass") & 
              (df1["qPCR_val_detail"]=="pass")]
    
    # Extract circRNAs from H23 cell line with count ≥ 5
    df1_h23 = df1[(df1["cell_line"]=="NCI-H23") & 
                  (df1["count_group_median"]=="count ≥ 5")]
    
    # Convert to set format
    validated_circles = set(df1_h23[["chr","start","end"]].apply(
        lambda x: (x[0][3:], int(x[1]), int(x[2])), axis=1).tolist())
    
    return validated_circles

def check_circle_validated(chr_name, start, end, validated_circles):
    """Check if a circRNA has been validated"""
    # Generate list of possible positions (allowing 2bp error)
    possible_positions = []
    for i in range(start-2, start+3):
        for j in range(end-2, end+3):
            if i < j:
                possible_positions.append((chr_name, i, j))
    
    # Check if any possible position exists in validation set
    return bool(set(possible_positions).intersection(validated_circles))

def analyze_tool_validation():
    """Analyze validation status of circRNAs detected by different tools"""
    # Read data
    df = pd.read_csv("/data6/circleRNA/circRNA_benchmarking/data/Supplementary_Table_4_all_circRNAs_treated.txt", 
                     sep="\t")
    
    # Get validated circRNAs
    validated_circles = get_validated_circles()
    total_validated = 234  # Total number of validated circRNAs
    
    # Get list of all tools
    all_tools = set()
    for tools in df['tool'].str.split('/'):
        all_tools.update(tools)
    
    # Calculate validation rate for each tool
    tool_stats = {}
    for tool in all_tools:
        # Get all circRNAs detected by this tool
        tool_circles = df[df['tool'].str.contains(tool)]
        
        total = len(tool_circles)
        validated = 0
        
        for _, circle in tool_circles.iterrows():
            if check_circle_validated(str(circle['chr']), 
                                   int(circle['start']), 
                                   int(circle['end']), 
                                   validated_circles):
                validated += 1
        
        recall = validated / total_validated if total_validated > 0 else 0
        tool_stats[tool] = {
            'total_predictions': total,
            'validated_found': validated,
            'recall': recall
        }
    
    # Convert to DataFrame
    results_df = pd.DataFrame.from_dict(tool_stats, orient='index')
    
    # Create output directory
    output_dir = "analysis_results"
    Path(output_dir).mkdir(exist_ok=True)
    
    # Save results to Excel
    results_df.to_excel(f"{output_dir}/tool_validation_stats.xlsx")
    
    # Create plot
    plt.figure(figsize=(12, 6))
    
    # Sort and create bar plot
    sorted_tools = sorted(tool_stats.items(), key=lambda x: x[1]['recall'], reverse=True)
    tools = [x[0] for x in sorted_tools]
    recalls = [x[1]['recall'] * 100 for x in sorted_tools]
    
    plt.bar(range(len(tools)), recalls)
    plt.xticks(range(len(tools)), tools, rotation=45, ha='right')
    plt.ylabel('Recall (%)')
    plt.title('Tool Recall (% of validated circRNAs found)')
    
    # Add specific values on top of each bar
    for i, v in enumerate(recalls):
        plt.text(i, v + 1, f'{v:.1f}%', ha='center')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/tool_validation_recall.png")
    plt.close()
    
    return results_df

def main():
    print("Starting tool validation analysis...")
    results = analyze_tool_validation()
    print("\nAnalysis complete! Results saved to analysis_results directory")
    print("\nTool validation statistics:")
    print(results)

if __name__ == "__main__":
    main() 