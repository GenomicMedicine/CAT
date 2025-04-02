import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def get_database_circles():
    """Read circRNAs from known database"""
    database_circles = set()
    with open("/data6/circleRNA/circRNA_benchmarking/data/details/circ_db_hg38.txt") as f:
        next(f)  # Skip header
        for line in f:
            line = line.strip().split("\t")
            database_circles.add((line[0], int(line[1]), int(line[2])))
    return database_circles

def check_circle_in_database(chr_name, start, end, database_circles):
    """Check if a circRNA exists in the database"""
    # Generate list of possible positions (allowing 2bp error)
    possible_positions = []
    for i in range(start-2, start+3):
        for j in range(end-2, end+3):
            if i < j:
                possible_positions.append((chr_name, i, j))
    
    # Check if any possible position exists in database
    return bool(set(possible_positions).intersection(database_circles))

def analyze_tool_precision():
    """Analyze precision of different tools"""
    # Read data
    df = pd.read_csv("/data6/circleRNA/circRNA_benchmarking/data/Supplementary_Table_4_all_circRNAs_treated.txt", 
                     sep="\t")
    
    df = df[df["cell_line"]=="NCI-H23"]
    # Filter rows where BSJ_count is greater than or equal to 5
    df = df[df["BSJ_count"]>=5]

    # Get circRNAs from database
    database_circles = get_database_circles()
    
    # Get list of all tools
    all_tools = set()
    for tools in df['tool'].str.split('/'):
        all_tools.update(tools)
    
    # Calculate precision for each tool
    tool_stats = {}
    for tool in all_tools:
        # Get all circRNAs detected by this tool
        tool_circles = df[df['tool'].str.contains(tool)]
        
        total = len(tool_circles)
        in_database = 0
        
        for _, circle in tool_circles.iterrows():
            if check_circle_in_database(str(circle['chr']), 
                                     int(circle['start']), 
                                     int(circle['end']), 
                                     database_circles):
                in_database += 1
        
        precision = in_database / total if total > 0 else 0
        tool_stats[tool] = {
            'total': total,
            'in_database': in_database,
            'precision': precision
        }
    
    # Convert to DataFrame
    results_df = pd.DataFrame.from_dict(tool_stats, orient='index')
    
    # Create output directory
    output_dir = "analysis_results"
    Path(output_dir).mkdir(exist_ok=True)
    
    # Save results to Excel
    results_df.to_excel(f"{output_dir}/tool_precision_stats.xlsx")
    
    # Create plot
    plt.figure(figsize=(12, 6))
    
    # Sort and create bar plot
    sorted_tools = sorted(tool_stats.items(), key=lambda x: x[1]['precision'], reverse=True)
    tools = [x[0] for x in sorted_tools]
    precisions = [x[1]['precision'] * 100 for x in sorted_tools]
    
    plt.bar(range(len(tools)), precisions)
    plt.xticks(range(len(tools)), tools, rotation=45, ha='right')
    plt.ylabel('Precision (%)')
    plt.title('Tool Precision (% of predictions in database)')
    
    # Add specific values on top of each bar
    for i, v in enumerate(precisions):
        plt.text(i, v + 1, f'{v:.1f}%', ha='center')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/tool_precision.png")
    plt.close()
    
    return results_df

def main():
    print("Starting tool precision analysis...")
    results = analyze_tool_precision()
    print("\nAnalysis complete! Results saved to analysis_results directory")
    print("\nTool precision statistics:")
    print(results)

if __name__ == "__main__":
    main() 