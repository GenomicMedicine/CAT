import pandas as pd
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
import tempfile
import re
def clean_sequence(sequence):
    """清理和验证序列"""
    if pd.isna(sequence):
        return None
    
    # 转换为字符串
    sequence = str(sequence).strip()
    
    # 移除非法字符，只保留 ATCG
    sequence = re.sub(r'[^ATCG]', '', sequence.upper())
    
    return sequence if sequence else None

def create_reference_fasta(excel_file, output_fasta):
    """从Excel文件创建参考序列FASTA文件"""
    print("Creating reference FASTA file...")
    df = pd.read_excel(excel_file)
    
    # 验证reference列存在
    if 'reference' not in df.columns:
        raise ValueError("Column 'reference' not found in Excel file")
    
    valid_sequences = 0
    invalid_sequences = 0
    
    with open(output_fasta, 'w') as f:
        for idx, row in df.iterrows():
            sequence = clean_sequence(row['reference'])
            if sequence:
                seq_id = f"seq_{idx}"
                f.write(f">{seq_id}\n{sequence}\n")
                valid_sequences += 1
            else:
                invalid_sequences += 1
                print(f"Warning: Invalid sequence at row {idx + 1}")
    
    print(f"Processed {valid_sequences + invalid_sequences} sequences:")
    print(f"  Valid sequences: {valid_sequences}")
    print(f"  Invalid sequences: {invalid_sequences}")
    
    if valid_sequences == 0:
        raise ValueError("No valid sequences found in the Excel file")
    
    return df

def run_mapping(ref_fasta, read1, read2, output_bam, temp_dir):
    """运行bowtie2进行比对"""
    print("Running bowtie2 mapping...")
    
    # 创建bowtie2索引
    index_base = os.path.join(temp_dir, "reference")
    print("Building bowtie2 index...")
    subprocess.run([
        'bowtie2-build',
        '--quiet',
        ref_fasta,
        index_base
    ], check=True)
    
    # 运行bowtie2比对,
    print("Running bowtie2 alignment...")
    # bowtie2_cmd = [
    #     'bowtie2',
    #     '-p', '100',
    #     '-x', index_base,
    #     '-1', read1,
    #     '-2', read2,
    #     '--very-sensitive-local',
    #     '--no-unal'
    # ]
    #     # 运行bowtie2比对,不允许错配
    # print("Running bowtie2 alignment...")
    # bowtie2_cmd = [
    #     'bowtie2',
    #     '-p', '100',              # 使用100个线程
    #     '-x', index_base,         # 索引前缀
    #     '-1', read1,              # 第一个fastq文件
    #     '-2', read2,              # 第二个fastq文件
    #     '--no-unal',             # 不输出未比对reads
    #     '--no-1mm-upfront',      # 禁用1错配的快速比对
    #     '--no-mixed',            # 禁用不完全配对
    #     '--no-discordant',       # 禁用不一致配对
    #     '-N', '0',               # seed中不允许错配
    #     '-L', '22',              # seed长度
    #     '--score-min', 'L,0,0',  # 最小分数阈值，不允许任何惩罚
    #     '--np', '0',             # N的惩罚为0
    #     '--rdg', '1000,1000',    # read gap开启和延伸的高惩罚
    #     '--rfg', '1000,1000',    # reference gap开启和延伸的高惩罚
    #     '--mp', '1000,1000'      # 错配的高惩罚
    # ]
    
    
    
        # 运行bowtie2比对，允许1个错配
    print("Running bowtie2 alignment...")
    # bowtie2_cmd = [
    #     'bowtie2',
    #     '-p', '100',              # 使用100个线程
    #     '-x', index_base,         # 索引前缀
    #     '-1', read1,              # 第一个fastq文件
    #     '-2', read2,              # 第二个fastq文件
    #     '--no-unal',             # 不输出未比对reads
    #     '--no-mixed',            # 禁用不完全配对
    #     '--no-discordant',       # 禁用不一致配对
    #     '-N', '1',               # seed中允许1个错配
    #     '-L', '22',              # seed长度
    #     '--score-min', 'L,-0.6,0',  # 允许一个错配的评分阈值
    #     '--mp', '6,6',           # 错配惩罚（降低惩罚以允许错配）
    #     '--rdg', '1000,1000',    # 保持gap开启和延伸的高惩罚
    #     '--rfg', '1000,1000',    # 保持gap开启和延伸的高惩罚
    # ]
    # 运行bowtie2比对，允许1个错配
    print("Running bowtie2 alignment...")
    bowtie2_cmd = [
        'bowtie2',
        '-p', '100',              # 使用100个线程
        '-x', index_base,         # 索引前缀
        '-1', read1,              # 第一个fastq文件
        '-2', read2,              # 第二个fastq文件
        '--score-min', 'C,-6,0',  # 更宽松的评分阈值

        '--rfg', '6,6',      # 降低gap惩罚但仍保持较高

    ]
    
    #--score-min C,{min_score},0 --rfg 6,6
    # 首先转换为SAM文件
    sam_file = os.path.join(temp_dir, "mapped.sam")
    with open(sam_file, 'w') as sam_out, open(os.path.join(temp_dir, "bowtie2.log"), 'w') as log_out:
        subprocess.run(bowtie2_cmd, stdout=sam_out, stderr=log_out, check=True)
    
    # 然后转换为BAM文件
    print("Converting SAM to BAM...")
    subprocess.run([
        'samtools', 'view',
        '-bS',
        '-o', output_bam,
        sam_file
    ], check=True)
    
    # 对BAM文件进行排序
    sorted_bam = os.path.join(temp_dir, "mapped.sorted.bam")
    print("Sorting BAM file...")
    subprocess.run([
        'samtools', 'sort', "-T",temp_dir+"/sort_temp/temp",
        '-@', '100',
        '-o', sorted_bam,
        output_bam
    ], check=True)
    
    # 替换原始BAM文件
    os.replace(sorted_bam, output_bam)
    
    # 创建BAM索引
    print("Creating BAM index...")
    subprocess.run(['samtools', 'index', output_bam], check=True)
    
    # 删除中间SAM文件
    os.remove(sam_file)
    
def count_reads(bam_file):
    """计算每个参考序列的比对reads数"""
    print("Counting mapped reads...")
    
    # 使用samtools idxstats获取比对统计
    result = subprocess.run(
        ['samtools', 'idxstats', bam_file],
        capture_output=True, text=True, check=True
    )
    
    # 解析结果
    counts = {}
    for line in result.stdout.strip().split('\n'):
        ref, length, mapped, unmapped = line.split()
        if ref != '*':  # 忽略未比对的reads
            counts[ref] = int(mapped)
    
    return counts
def main():
    # 设置输入输出路径
    excel_file = "/data4/H23_amplicon/candidate.xlsx"
    read1 = "/data4/H23_amplicon/clean_data/clean_1.fq.gz"
    read2 = "/data4/H23_amplicon/clean_data/clean_2.fq.gz"
    
    # 创建固定的临时目录
    temp_dir = "/data4/H23_amplicon/temp_mapping2"
    os.makedirs(temp_dir, exist_ok=True)
    print(f"Working directory: {temp_dir}")
    
    ref_fasta = os.path.join(temp_dir, "reference.fa")
    output_bam = os.path.join(temp_dir, "mapped.bam")
    
    try:
        # 读取Excel并创建参考序列
        print(f"Reading Excel file: {excel_file}")
        df = create_reference_fasta(excel_file, ref_fasta)
        print(f"Reference FASTA created: {ref_fasta}")
        
        # 运行比对
        run_mapping(ref_fasta, read1, read2, output_bam, temp_dir)
        
        # 统计比对结果
        read_counts = count_reads(output_bam)
        
        # 将结果添加到DataFrame
        df['mapped_reads'] = [read_counts.get(f"seq_{i}", 0) for i in range(len(df))]
        
        # 保存结果
        output_excel = "/data4/H23_amplicon/candidate_with_expression_mis1.xlsx"
        df.to_excel(output_excel, index=False)
        print(f"Results saved to {output_excel}")
        
        # 打印统计信息
        print("\nMapping statistics:")
        print(f"Total sequences: {len(df)}")
        print(f"Sequences with mapped reads: {sum(1 for count in df['mapped_reads'] if count > 0)}")
        print(f"Total mapped reads: {sum(df['mapped_reads'])}")
        
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        raise
    finally:
        # 询问是否删除临时文件
        response = input(f"\nDelete temporary files in {temp_dir}? (y/n): ")
        if response.lower() == 'y':
            import shutil
            shutil.rmtree(temp_dir)
            print("Temporary files cleaned up")
        else:
            print(f"Temporary files kept in: {temp_dir}")
            print("Files:")
            for root, dirs, files in os.walk(temp_dir):
                for file in files:
                    print(f"  {os.path.join(root, file)}")

if __name__ == "__main__":
    main() 