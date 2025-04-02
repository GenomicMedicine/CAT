import pandas as pd
import re
import pickle
import pyliftover
import bisect
from collections import Counter
import logomaker
import warnings
pd.options.mode.chained_assignment = None

def get_dict(txt_path):
    """Load gene dictionary from GO annotation file"""
    dic_gene = {}
    with open(txt_path) as f:
        for line in f:
            line = line.strip().split(")")
            genelist = line[1].strip().split("\t")
            dic_gene[line[0]] = tuple(genelist)
    dic_gene2 = {}
    for i in dic_gene.keys():
        for j in dic_gene[i]:
            if j not in dic_gene2.keys():
                dic_gene2[j] = [i]
            else:
                dic_gene2[j].append(i)
                
    return dic_gene2

def add_column(dic_gene, df, name):
    """Add GO annotation column to dataframe"""
    df[name] = df.apply(lambda row: dic_gene[row[0]] if row[0] in dic_gene.keys() else "Non", axis=1)
    return df

def find_database(result_df, compare_list, new_row):
    """Find circRNAs in database with position matching"""
    # Convert position columns to integer type
    result_df['left boundary of shs area1'] = result_df['left boundary of shs area1'].astype(int)
    result_df['right boundary of shs area1'] = result_df['right boundary of shs area1'].astype(int)
    result_df['left boundary of shs area2'] = result_df['left boundary of shs area2'].astype(int)
    result_df['right boundary of shs area2'] = result_df['right boundary of shs area2'].astype(int)
    
    # Get position lists
    result_df_left1 = result_df['left boundary of shs area1'].tolist()
    result_df_right1 = result_df['right boundary of shs area1'].tolist()
    result_df_left2 = result_df['left boundary of shs area2'].tolist()
    result_df_right2 = result_df['right boundary of shs area2'].tolist()  
    result_df_chr = result_df["chr"].astype(str).tolist()
    
    # Initialize new column
    result_df[new_row] = "no"
    found_circRNA = []

    # Check each circRNA position
    for i in range(len(result_df)):
        list3 = []
        list1 = list(range(result_df_left1[i]-5, result_df_right1[i]+5))
        list2 = list(range(result_df_left2[i]-5, result_df_right2[i]+5))
        chr_name = result_df_chr[i]
        for j in list1:
            for k in list2:
                if j < k:
                    list3.append((chr_name, j, k))
                else:
                    list3.append((chr_name, k, j))

        # If any position matches database, mark as found
        if set(list3).intersection(compare_list):
            result_df.iloc[i,-1] = "yes"
            found_circRNA.append(next(iter(set(list3).intersection(compare_list))))
    found_circRNA = set(found_circRNA)

    return result_df, found_circRNA

def main():
    # Load configuration and paths
    cat_res = "/data6/circleRNA/circRNA_benchmarking/data/circleRNA_predict_result1_14.tsv"
    project_dir = "/data6/H23_public466"
    output_file2 = project_dir+"/output_classify_result_12_19_end2end.xlsx"
    writer = pd.ExcelWriter("/data6/H23_public466/466_newcat_test/circleRNA_predict_result_paper_supplement_cutoff5.xlsx")

    # Load and filter validation data
    df1 = pd.read_csv("/data6/circleRNA/circRNA_benchmarking/data/Supplementary_Table_3_selected_circRNAs.txt", sep="\t")
    df1 = df1[(df1["amp_seq_val"]=="pass") & 
              (df1["amp_seq_val_detail"]=="pass") & 
              (df1["compound_val"]=="pass") & 
              (df1["RR_val"]=="pass") & 
              (df1["RR_val_detail"]=="pass") & 
              (df1["qPCR_val"]=="pass") & 
              (df1["qPCR_val_detail"]=="pass")]

    # Filter by count and cell line
    df1_1 = df1[df1["count_group_median"]=="count ≥ 5"]
    df1_h23 = df1_1[df1_1["cell_line"]=="NCI-H23"]
    
    # Get ground truth data
    ground_truth = df1_h23[["chr","start","end"]].apply(
        lambda x: (x[0][3:], int(x[1]), int(x[2])), axis=1).tolist()

    # Load and process circlebank data
    circlebank_df = pd.read_csv("/data6/circleRNA/circlebank/circBank_circrna_annotation.txt", sep="\t")
    
    # Load database circles
    circle_result = []
    with open("/data6/circleRNA/circRNA_benchmarking/data/details/circ_db_hg38.txt") as f:
        next(f)  # Skip header
        for line in f:
            line = line.strip().split("\t")
            circle_result.append((line[0][3:], int(line[1]), int(line[2])))
    circle_result = set(circle_result)

    # Load and filter prediction data
    circle_predict_df = pd.read_csv(cat_res, sep="\t")
    circle_predict_df = circle_predict_df[circle_predict_df["LenShs"]<=20]
    circle_predict_df = circle_predict_df[circle_predict_df["ReadCount"]>=5]

    # Split by circRNA types
    type_dfs = {
        'canonical': circle_predict_df[circle_predict_df["CircType"].str.contains("canonical")],
        'interior': circle_predict_df[circle_predict_df["CircType"].str.contains("interior")],
        'lariat': circle_predict_df[circle_predict_df["CircType"].str.contains("lariat")],
        'intergenic': circle_predict_df[circle_predict_df["CircType"].str.contains("intergenic")],
        'antisense': circle_predict_df[circle_predict_df["CircType"].str.contains("antisense")],
        'partial': circle_predict_df[circle_predict_df["CircType"].str.contains("partial")]
    }

    # Process each type
    processed_dfs = {}
    for type_name, df in type_dfs.items():
        # Check database presence
        df_known, found_circs = find_database(df, circle_result, "in known database")
        # Check validation status
        df_validated, _ = find_database(df_known, ground_truth, "is validated")
        processed_dfs[type_name] = df_validated
        
        # Print statistics
        print(f"\n{type_name} circRNA statistics:")
        print("In known database:", df_validated["in known database"].value_counts(), 
              df_validated["in known database"].value_counts()/len(df_validated))
        print("Is validated:", df_validated["is validated"].value_counts(),
              df_validated["is validated"].value_counts()/len(df_validated))

    # Save results to Excel
    for type_name, df in processed_dfs.items():
        df.to_excel(writer, sheet_name=f"circle_{type_name}")
    writer.save()

    # Process tool comparison
    df_tool = pd.read_csv("/data6/circleRNA/circRNA_benchmarking/data/Supplementary_Table_4_all_circRNAs_treated.txt", sep="\t")
    df_tool = df_tool[df_tool["cell_line"]=="NCI-H23"]
    df_tool = df_tool[df_tool["count_group"]=="count ≥ 5"]

    # Analyze each tool
    all_tools = set()
    for tools in df_tool['tool'].str.split('/'):
        all_tools.update(tools)

    for tool in all_tools:
        tool_circles = df_tool[df_tool['tool'].str.contains(tool)]
        total = len(tool_circles)

        # Process positions
        tool_circles["left boundary of shs area1"] = tool_circles["start"].astype(int)
        tool_circles["right boundary of shs area1"] = tool_circles["start"].astype(int)
        tool_circles["left boundary of shs area2"] = tool_circles["end"].astype(int)
        tool_circles["right boundary of shs area2"] = tool_circles["end"].astype(int)
        tool_circles["chr"] = tool_circles["chr"].str.slice(3)

        # Check validation
        tool_circles_known, tool_found_circRNA = find_database(tool_circles, ground_truth, "in known database")
        print(f"Tool {tool} found {len(tool_found_circRNA)} validated circRNAs",
              f"(Recall: {len(tool_found_circRNA)/len(ground_truth):.2%})")

if __name__ == "__main__":
    main()


