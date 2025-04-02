import argparse
import pysam
from Bio import SeqIO
from tqdm import tqdm
import pandas as pd
import time
import os
from collections import defaultdict
import arrow
import pyensembl
import CAT_thread as ct

global_time1 = time.time()
def star_map(in_fq, genome_lib_dir, output_dir, out_bam_name, out_unmap_fastx=True, soft=False, max_out_filter_multimap=10, sjdb=100,
             max_out_filter_mismatch=1,  min_score_filter=0.66,
             strand=False,  threads=150, alo_mis_match=True):
    if out_bam_name != '':
        bam_dir = os.path.join(os.path.dirname(str(output_dir)), out_bam_name)
    else:
        bam_dir = os.path.join(str(output_dir), 'Aligned.sortedByCoord.out.bam')
    if not os.path.exists(bam_dir):
        fq_str = ' '.join(in_fq)
        out_dir = str(output_dir) if str(output_dir)[-1] == '/' else str(output_dir) + '/'

        os_str = f'STAR --runThreadN {threads} ' + \
                 f'--genomeDir {str(genome_lib_dir)} ' + \
                 f'--readFilesIn {fq_str} ' + \
                 f'--outFileNamePrefix {out_dir} ' + \
                 f'--outSAMtype BAM SortedByCoordinate ' + \
                 f'--outFilterMultimapNmax {max_out_filter_multimap} ' + \
                 f'--alignSJoverhangMin {sjdb} ' + \
                 f'--outBAMsortingThreadN {threads // 4} ' + \
                 f'--quantMode GeneCounts '
        if not soft:
            os_str += f'--alignEndsType EndToEnd '
        if not alo_mis_match:
            os_str += f'--outFilterScoreMin {min_score_filter} ' + f'--outFilterMismatchNmax {max_out_filter_mismatch} '
        if '.gz' in fq_str or '.gzip' in fq_str or '.zip' in fq_str:
            os_str += f'--readFilesCommand zcat '
        if out_unmap_fastx:
            os_str += f'--outReadsUnmapped Fastx '
        if strand:
            os_str += f'--readStrand Reverse'
        os.system(os_str)
        if out_bam_name != '':
            os.system(f'mv {os.path.join(out_dir, "Aligned.sortedByCoord.out.bam")} {bam_dir}')
        print('indexing...')
        os.system(f'samtools index {bam_dir}')
    else:
        print(f"Existing alignment result has been read at {bam_dir}")

def fetch_refer(fas,chr1,locate1,refer_length):
    refer_length = int(refer_length)
    locate1 = int(locate1)
    return fas.fetch(chr1,locate1,locate1+refer_length).upper()

def split_file(file_path, n):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    total_lines = len(lines)
    lines_per_file = total_lines // n

    for i in range(n):
        start_idx = i * lines_per_file
        end_idx = start_idx + lines_per_file if i < n - 1 else total_lines

        base_name = os.path.basename(file_path).split('.')[0]
        dir_name = os.path.join(os.path.dirname(file_path),base_name)
        if not os.path.exists(dir_name):
            os.makedirs(dir_name, exist_ok=True)
        file_name = os.path.join(dir_name, f'part_{i}.txt')

        with open(file_name, 'w') as part_file:
            part_file.writelines(lines[start_idx:end_idx])

    return [f'{dir_name}/part_{i}.txt' for i in range(n)]

# reverse complement if the map to the complement chain of the reference
def reverse_complement(str):
    str = str.upper()
    return str[::-1].replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g').upper()

def parse_sam(samfile_path):
    mapped_reads = defaultdict(list)
    with pysam.AlignmentFile(samfile_path, "r") as samfile,tqdm(desc='parsing sam') as pbar:
        for read in samfile:
            if not read.is_unmapped:
                read_info = {
                    "chromosome": read.reference_name,
                    "position": read.reference_end if not read.is_reverse else read.reference_start,
                    "is_reverse": read.is_reverse
                }
                mapped_reads[read.query_name].append(read_info)
                pbar.update(1)
    return mapped_reads

def parse_sam_paired(samfile_path):
    mapped_reads_m1 = defaultdict(list)
    mapped_reads_m2 = defaultdict(list)
    with pysam.AlignmentFile(samfile_path, "r") as samfile,tqdm(desc='parsing sam') as pbar:
        for read in samfile:
            if not read.is_unmapped:
                read_info = {
                    "chromosome": read.reference_name,
                    "position": read.reference_end if not read.is_reverse else read.reference_start,
                    "is_reverse": read.is_reverse
                }
                if read.is_read1:
                    mapped_reads_m1[read.query_name].append(read_info)
                else:
                    mapped_reads_m2[read.query_name].append(read_info)
                pbar.update(1)
    return mapped_reads_m1,mapped_reads_m2
def get_combined_info(left_mapped_reads, right_mapped_reads):
    combined_info = []
    with tqdm(total=len(left_mapped_reads),desc='searching discordant map') as pbar:
        for read_name in left_mapped_reads:
            if read_name in right_mapped_reads:
                left_info_list = left_mapped_reads[read_name]
                right_info_list = right_mapped_reads[read_name]
                for left_info in left_info_list:
                    for right_info in right_info_list:
                        if left_info["is_reverse"] == right_info["is_reverse"] and right_info["chromosome"] == left_info["chromosome"]:
                            reverse = left_info["is_reverse"]
                        else:
                            continue
                        if not left_info["chromosome"] == "MT" and (reverse == (left_info['position'] < right_info['position']) and abs(left_info['position'] - right_info['position']) < 50000):
                            combined_info.append((read_name, left_info, right_info))
            pbar.update(1)
    return combined_info

def find_fusion_position(seq, left_refer, right_refer,left_flag,right_flag,left_position,right_position,kmis,anchor_len):
    fusion_index = None

    mismatch_count = -1
    i =0
    for i in range(0, len(left_refer)):
        # print(i)
        if seq[i] != left_refer[i]:
            mismatch_count += 1
        if mismatch_count == kmis:
            break
    fusion_index = i

    while (seq[fusion_index - 1] != left_refer[fusion_index - 1]):
        fusion_index -= 1
        mismatch_count -= 1
    remaining_seq = str(seq[fusion_index:])

    len_seq = len(seq)

    loop_len = min(len(remaining_seq), len(right_refer))
    for i in range(0, loop_len):
        if remaining_seq[-i - 1] != right_refer[-i - 1]:
            mismatch_count += 1

    if mismatch_count > kmis:
        return None
    else:
        if (left_flag == True):
            left_fusion_position = left_position - fusion_index +anchor_len
        else:
            left_fusion_position = left_position + fusion_index -anchor_len
        if (right_flag == True):
            right_fusion_position = right_position - fusion_index +len_seq
        else:
            right_fusion_position = right_position + fusion_index -len_seq

        return (left_fusion_position, right_fusion_position)
def process_combined_pair(fusion_position_dir,combined_info,seq_dir,dna_primary_assembly_fa,kmis,anchor_len,mate_no,starnd_type):
    if mate_no == 1 and starnd_type=='F2R1':
        revese_result = True
    elif mate_no == 2 and starnd_type=='F1R2':
        revese_result = True
    else:
        revese_result = False
    fusion_positions = []
    fa = pysam.FastaFile(dna_primary_assembly_fa)
    print("start dump seq to dict...")
    write_mode = "w" if not os.path.exists(fusion_position_dir) else "a"
    # if write_mode =="a":
    #     filee = open(fusion_position_dir, "r")
    #     lines = { i.strip().split('\t')[-1]:i.strip().split('\t') for i in filee.readlines()}
    #     filee.close()
    # else:
    #     lines = False
    idx =0
    record_dict = SeqIO.to_dict(SeqIO.parse(seq_dir,"fastq"))

    with open(fusion_position_dir, write_mode) as f, tqdm(total=len(combined_info), desc='finding junction site') as pbar:
        for read_name, left_info, right_info in combined_info:
            seq = str(record_dict[read_name].seq).upper()  # Replace this with the actual sequence)
            seq_len = len(seq)
            left_flag = left_info["is_reverse"]
            right_flag = right_info["is_reverse"]
            try:
                if left_flag == True:  # reverse if flag==16
                    y_refer = str(fetch_refer(fa,left_info["chromosome"],max(0,left_info["position"]+anchor_len-seq_len),seq_len))
                    y_refer = reverse_complement(y_refer)
                else:
                    y_refer = str(fetch_refer(fa,left_info["chromosome"],max(0,left_info["position"]-anchor_len),seq_len))
                if right_flag == True:  # reverse if flag==16
                    x_refer = str(fetch_refer(fa,right_info["chromosome"],max(0,right_info["position"]),seq_len))
                    x_refer = reverse_complement(x_refer)
                else:
                    x_refer = str(fetch_refer(fa,right_info["chromosome"],max(0,right_info["position"]-seq_len),seq_len))
                left_position = left_info["position"]
                right_position = right_info["position"]
            except KeyError:
                print(f"KeyError: {read_name} {left_info['chromosome']} {right_info['chromosome']}")
                continue
            # if lines and read_name not in lines and :
            #     continue
            try:
                fusion_position = find_fusion_position(seq, str(y_refer).upper(), str(x_refer).upper(),left_flag,right_flag,left_position,right_position,kmis,anchor_len)
            except IndexError:
                print(f"IndexError: {read_name} {str(y_refer).upper()} { str(x_refer).upper()}")
                continue
            # print(seq, y_refer, x_refer)
            if  fusion_position:
                if not (left_info["chromosome"] == right_info["chromosome"] and abs(fusion_position[0] - fusion_position[1]) < 40):
                    strands = "--" if left_flag else "++"
                    l_pos, r_pos = fusion_position[0],fusion_position[1]
                    chrom1, chrom2 = left_info["chromosome"], right_info["chromosome"]
                    if revese_result:
                        strands = "--" if strands == "++" else "++"
                        l_pos, r_pos = r_pos, l_pos
                        chrom1, chrom2 = chrom2, chrom1
                    fusion_positions.append((chrom1, l_pos, chrom2, r_pos, strands))
                    f.write(f"{chrom1}\t{l_pos}\t{chrom2}\t{r_pos}\t{strands}\n")
            pbar.update(1)
            idx+=1
    return fusion_positions

def remove_duplicate_sublists(lst,read_from_file):
    if read_from_file:
        unique_sublists = []
        with open(lst, "r") as f:
            for line in f:
                unique_sublists.append(tuple(line.strip().split("\t")))
        unique_sublists = set(unique_sublists)
    else:
        unique_sublists = set(tuple(sublist) for sublist in lst)
    return list(unique_sublists)
def common_prefix(str1, str2):
    n1, n2 = len(str1), len(str2)
    i, j = 0, 0
    while i <= n1 - 1 and j <= n2 - 1:
        if str1[i] != str2[j]:
            break
        i, j = i + 1, j + 1
    return str2[:j]
def common_suffix(str1, str2):
    n1, n2 = len(str1), len(str2)
    i, j = n1 - 1, n2 - 1
    while i >= 0 and j >= 0:
        if str1[i] != str2[j]:
            break
        i, j = i - 1, j - 1
    return(str1[i+1:])


def remap_ref(chimeric_ref_dir, dna_primary_assembly_fa, finial_list, length_seq):
    fa = pysam.FastaFile(dna_primary_assembly_fa)
    shs_dict = defaultdict(list)
    with tqdm(total=len(finial_list), desc='collecting shs info') as pbar:
        for line in finial_list:
            pos1, pos2 = int(line[1]), int(line[3])
            left_of_leftrefer = fetch_refer(fa, line[0], max(0, pos1 - length_seq), length_seq)
            right_of_rightrefer = fetch_refer(fa, line[2], pos2, length_seq)
            left_of_rightrefer = fetch_refer(fa, line[2], max(0, pos2 - length_seq), length_seq)
            right_of_leftrefer = fetch_refer(fa, line[0], pos1, length_seq)
            reverse_left_of_leftrefer = reverse_complement(right_of_leftrefer)
            reverse_left_of_rightrefer = reverse_complement(right_of_rightrefer)
            reverse_right_of_rightrefer = reverse_complement(left_of_rightrefer)
            reverse_right_of_leftrefer = reverse_complement(left_of_leftrefer)

            if line[4] == "++":
                common_left = common_suffix(left_of_leftrefer, left_of_rightrefer)
                common_right = common_prefix(right_of_leftrefer, right_of_rightrefer)
                shs = common_left + common_right
                motif_20bp = left_of_rightrefer[-10 - len(common_left):] + right_of_rightrefer[:len(
                    common_right)] + shs + left_of_leftrefer[
                                           len(left_of_leftrefer) - (len(common_left)):] + right_of_leftrefer[
                                                                                           :10 + len(common_right)]
                new_pos1_left = pos1 - len(common_left)
                new_pos2_left = pos2 - len(common_left)
                new_left_of_leftrefer = fetch_refer(fa, line[0], max(0, new_pos1_left - length_seq), length_seq)
                new_left_of_rightrefer = fetch_refer(fa, line[2], max(0, new_pos2_left - length_seq), length_seq)
                seq = new_left_of_leftrefer + new_left_of_rightrefer
            elif line[4] == "--":
                common_left = common_suffix(reverse_left_of_leftrefer, reverse_left_of_rightrefer)
                common_right = common_prefix(reverse_right_of_leftrefer, reverse_right_of_rightrefer)
                shs = common_left + common_right
                motif_20bp = reverse_left_of_rightrefer[-10 - len(common_left):] + reverse_right_of_rightrefer[:len(
                    common_right)] + shs + reverse_left_of_leftrefer[
                                           len(left_of_leftrefer) - len(common_left):] + reverse_right_of_leftrefer[
                                                                                         :10 + len(common_right)]
                new_pos1_left = pos1 - len(common_left)
                new_pos2_left = pos2 - len(common_left)
                new_reverse_left_of_leftrefer = reverse_complement(fetch_refer(fa, line[0], pos1, length_seq))
                new_reverse_left_of_rightrefer = reverse_complement(fetch_refer(fa, line[2], pos2, length_seq))
                seq = new_reverse_left_of_leftrefer + new_reverse_left_of_rightrefer
            shs_key = f"{line[0]}:{pos1-len(common_left)}-{pos1+len(common_right)}_{line[2]}:{pos2-len(common_left)}-{pos2+len(common_right)}"
            
            row_res = list(line) + [line[4], shs, motif_20bp, 
                                   pos1 - len(common_left), pos1 + len(common_right),
                                   pos2 - len(common_left), pos2 + len(common_right),
                                   seq]
            shs_dict[shs_key].append(row_res)
            pbar.update(1)

    with tqdm(total=len(shs_dict), desc='writing merged references') as pbar, open(chimeric_ref_dir, "w") as f:
        for shs_key, entries in shs_dict.items():
            str_list = [str(s) for s in entries[0]]
            merged_id = f"{'_'.join(str_list)}"
            fasta_entry = f">{merged_id}\n{str_list[-1]}\n"
            f.write(fasta_entry)
            pbar.update(1)
            
    return



def count_alignments(sam_file,filename,cut_off_count):
    alignment_counts = defaultdict(lambda: [0, 0])
    with tqdm(desc='counting alignments') as pbar, open(sam_file, 'r') as f:
        for line in f:
            pbar.update(1)
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            flag = int(fields[1])
            if not (flag & 0x4) : 
                if not (flag & 0x10) and not (flag & 0x400): 
                    ref_name = fields[2]
                    alignment_counts[ref_name][0] += 1
                    if int(fields[5][:-1]) == len(fields[9]):  
                        alignment_counts[ref_name][1] += 1

    with tqdm(desc='saving csv') as pbar,open(filename, 'w', newline='') as file: 
        for key, value in alignment_counts.items():
            pbar.update(1)
            if value[0] >= cut_off_count and value[1] >= 1:
                file.write(f'{key}_{value[0]}\n')
    return alignment_counts

def multi_rpkm(in_bed, bam_dir, out_depth_txt, threads):
    time1 = time.time()
    print("start count gene expression...", end='')

    if threads > 1:
        bed_list = split_file(in_bed, threads)
        arg_list = []
        res_list = []
        base_name = os.path.basename(out_depth_txt).split('.')[0]
        dir_name = os.path.join(os.path.dirname(out_depth_txt), base_name)
        if not os.path.exists(dir_name):
            os.makedirs(dir_name, exist_ok=True)
        for bed, idx in zip(bed_list, range(len(bed_list))):
            res_name = os.path.join(dir_name, f'part_{idx}_rpkm.txt')
            arg_list.append((bed, bam_dir, res_name,))
            res_list.append(res_name)
        funs = ct.calculate_rpkm

        ct.run_in_parallel([funs] * len(arg_list), arg_list, len(arg_list))

        with open(out_depth_txt, 'w') as outfile:
            for file in res_list:
                with open(file, 'r') as infile:
                    outfile.write(infile.read())
    else:
        ct.calculate_rpkm(in_bed, bam_dir, out_depth_txt)
    time2 = time.time()
    print("rpkm done :", time2 - time1, 's')

def parse_gtf(gtf_file_dir,out_exon_bed):
    exondic = {}
    cdsdic = {}
    gtf_file = open(gtf_file_dir)
    gtf_line = gtf_file.readlines()
    gtf_file.close()
    annotgene = {}
    with tqdm(total=len(gtf_line), desc='parsing gtf') as pbar, open(out_exon_bed, 'w') as bed_file:
        for line in gtf_line:
            pbar.update(1)
            if line.startswith('#'):
                continue
            info = line.split('\t')
            chromo = info[0]
            start = int(info[3])
            end = int(info[4])
            strand = info[6]
            geneinfo = info[8].split('; ')
            genename = ''
            geneid = ''
            transname = ''
            for item in geneinfo:
                if genename != '' and transname != '':
                    break
                if item.startswith('gene_name'):
                    genename = item[11:-1]
                elif item.startswith('transcript_id'):
                    transname = item[15:-1]
                elif item.startswith('gene_id'):
                    geneid = item[9:-1]
            if genename == '':
                genename = geneid
            if info[2] == 'exon':
                if genename != '':
                    if genename not in exondic:
                        exondic[genename] = []
                    exondic[genename].append([start, end, strand, info[2], transname])
            elif info[2] == 'gene':
                biotype = ''
                for item in geneinfo:
                    if item.startswith('gene_biotype'):
                        biotype = str(item[14:-1]).replace('"', '').replace(';', '')
                annotgene[genename] = [[start, end, strand, biotype]]
                bed_file.writelines(chromo + '\t' + str(start) + '\t' + str(end) + '\t' + genename + '\n')
    return exondic, cdsdic, annotgene


def annotate_circRNA(counting_res, gtf_dir, gene_bed, cut_off_count,mid_file,circ_bed):
    with open(counting_res, "r") as infile:
        lines = infile.readlines()
    result_list = []
    annot_obj = pyensembl.Genome(reference_name='my_genome', annotation_name='my_genome_features', gtf_path_or_url=gtf_dir)
    annot_obj.index()
    exondic, cds, annotgene = parse_gtf(gtf_dir,gene_bed)
    with tqdm(desc='classifying',total=len(lines)) as pbar, open(mid_file, 'w', newline='') as outfile, open(circ_bed, 'w', newline='') as circout:
        inter, diff, same = 0, 0, 0
        for line in lines:
            pbar.update(1)
            data = line.strip().split("_")
            readcount = int(data[-1])
            circ_strand = data[4][0] 
            if readcount >= cut_off_count:
                chr1, pos1, chr2, pos2 = data[0], [int(data[8]),int(data[9])], data[2], [int(data[10]),int(data[11])]
                circ_strands, shs, ref_Seq = data[4][0], data[6] if len(data[6])!=0 else 'Na', data[12]
                shs_bod1 = f"{data[8]}-{data[9]}"
                shs_bod2 = f"{data[10]}-{data[11]}"
                if pos1[0] <= pos1[1]:
                    l_cpos1, h_cpos1 = pos1[0], pos1[1]

                else:
                    l_cpos1, h_cpos1 = pos1[1], pos1[0]
                if pos2[0] <= pos2[1]:
                    l_cpos2, h_cpos2 = pos2[0], pos2[1]
                else:
                    l_cpos2, h_cpos2 = pos2[1], pos2[0]
                gene_list1 = annot_obj.gene_names_at_locus(contig=chr1, position=pos1[0])
                gene_list1 += annot_obj.gene_names_at_locus(contig=chr1, position=pos1[1])
                gene_list1 =list(set(gene_list1))
                gene_list2 = annot_obj.gene_names_at_locus(contig=chr2, position=pos2[0])
                gene_list2 += annot_obj.gene_names_at_locus(contig=chr2, position=pos2[1])
                gene_list2 =list(set(gene_list2))
                if '' in gene_list1:
                    gene_list1.remove('')
                if '' in gene_list2:
                    gene_list2.remove('')
                cric_info_list = [chr1, int(data[1]), int(data[3]), circ_strands,shs,shs_bod1,shs_bod2,ref_Seq]
                intersection_gene = [gene for gene in list(set(gene_list1) & set(gene_list2)) if gene != '']

                if not gene_list1 and not gene_list2:
                    gene_name = f"intergenic"
                    gstrands, gtypes = ".", "Na"
                    inter += 1
                elif not intersection_gene:
                    gene_name1 = f"{gene_list1[0]}" if gene_list1  else 'intergenic'
                    gene_name2 = f"{gene_list2[0]}" if gene_list2  else 'intergenic'
                    if not (gene_name1 == 'intergenic' or gene_name2 =='intergenic'):
                        continue
                    try:
                        gene_name = f"{gene_name1}::{gene_name2}"
                        gstrands = annotgene[gene_list1[0]][0][2] if gene_list1 else annotgene[gene_list2[0]][0][2] if gene_list2 else "."
                        gtypes = annotgene[gene_list1[0]][0][3] if gene_list1 else annotgene[gene_list2[0]][0][3] if gene_list2 else "Na"
                    except KeyError:
                        gene_name = f"{gene_name1}::{gene_name2}"
                        gstrands, gtypes = ".", "Na"
                    diff += 1
                else:
                    try:
                        same += 1
                        gene_name = ''
                        for tis_gene in intersection_gene:
                            if circ_strand == annotgene[tis_gene][0][2]:
                                gene_name = tis_gene
                                gstrands, gtypes = annotgene[tis_gene][0][2], annotgene[tis_gene][0][3]
                        if gene_name == '':
                            gene_name = intersection_gene[0]
                            gstrands, gtypes = annotgene[gene_name][0][2], annotgene[gene_name][0][3]
                    except KeyError:
                        gene_name = intersection_gene[0]
                        gstrands, gtypes = ".", "Na"
                anno_res = [gene_name,f"{circ_strand}/{gstrands}",gtypes,'intergenic',[]]

                for this_gene in intersection_gene:
                    try:
                        trs_splice_info = dict() 
                        trs_length = dict()  
                        for item in exondic[this_gene]:
                            this_trs_name = item[4]
                            l_exon, h_exon = item[0], item[1]
                            this_trans_strand = item[2]
                            this_trans_tpye = item[3]
                            if not this_trs_name in trs_splice_info:
                                trs_splice_info[this_trs_name] = ['n/a', 'n/a', this_trans_strand]
                                trs_length[this_trs_name] = [l_exon, h_exon, this_trans_strand]

                            trs_length[this_trs_name][0] = min(trs_length[this_trs_name][0], l_exon)
                            trs_length[this_trs_name][1] = max(trs_length[this_trs_name][1], h_exon)

                            for l_cpos, h_cpos, idx in [(l_cpos1, h_cpos1,0), (l_cpos2, h_cpos2,1)]:
                                if 'splice-site' in trs_splice_info[this_trs_name][idx]:
                                    continue
                                if l_cpos-1 <= h_exon and l_exon <= h_cpos+1: 
                                    rcd_splice_info = this_trans_tpye
                                    if l_cpos-1  <= l_exon <= h_cpos+1  :
                                        if this_trans_strand == '+':
                                            rcd_splice_info += "/5'splice-site"
                                        else:
                                            rcd_splice_info += "/3'splice-site"
                                    elif l_cpos-1  <= h_exon <= h_cpos+1 :
                                        if this_trans_strand == '+':
                                            rcd_splice_info += "/3'splice-site"
                                        else:
                                            rcd_splice_info += "/5'splice-site"
                                    trs_splice_info[this_trs_name][idx] = rcd_splice_info
                    except KeyError:
                        continue

                    nc, lar, sp1, sp2 = [], [], [], []
                    for trs_name, trs_info in trs_splice_info.items():
                        na_count = '_'.join(trs_info[:-1]).count('n/a')
                        sp_count = '_'.join(trs_info[:-1]).count('splice-site')
                        if na_count != 0: 
                            if l_cpos1-1 <= trs_length[trs_name][1] and trs_length[trs_name][0] <= h_cpos1+1 and l_cpos2-1 <= trs_length[trs_name][1] and trs_length[trs_name][0] <= h_cpos2+1:
                                trs_info = ["intron" if x == 'n/a' else x for x in trs_info]
                            else:
                                continue
                        if circ_strand != gstrands:
                            anno_res[3] = 'antisense'
                            if sp_count ==2:
                                l_flg = 'exon' if "3'" in trs_info[0] else 'intron'
                                r_flg = 'exon' if "5'" in trs_info[1] else 'intron'
                                this_types = f'canonical-AS ({l_flg}-{r_flg})'
                                sp2.append([trs_name,*trs_info,this_types])
                            elif sp_count==1: 
                                if trs_info[0] == 'intron' and "3'" in trs_info[1]:
                                    this_types = 'lariat-AS'
                                    lar.append([trs_name,*trs_info,this_types])
                                elif trs_info[1] == 'intron' and "5'" in trs_info[0]:
                                    this_types = 'lariat-AS'
                                    lar.append([trs_name, *trs_info, this_types])
                                else:
                                    this_types = 'partial-AS'
                                    sp1.append([trs_name,*trs_info,this_types])
                            else:
                                this_types = 'interior-AS'
                                nc.append([trs_name,*trs_info,this_types])
                                anno_res[4] = list(sp2 + lar + sp1 + nc)
                        else:
                            if sp_count ==2: 
                                l_flg = 'exon' if "3'" in trs_info[0] else 'intron'
                                r_flg = 'exon' if "5'" in trs_info[1] else 'intron'
                                this_types = f'canonical ({l_flg}-{r_flg})'
                                sp2.append([trs_name,*trs_info,this_types])
                            elif sp_count==1: 
                                if trs_info[0] == 'intron' and "3'" in trs_info[1]:
                                    this_types = 'lariat'
                                    lar.append([trs_name,*trs_info,this_types])
                                elif trs_info[1] == 'intron' and "5'" in trs_info[0]:
                                    this_types = 'lariat'
                                    lar.append([trs_name, *trs_info, this_types])
                                else:
                                    this_types = 'partial'
                                    sp1.append([trs_name,*trs_info,this_types])
                            else:
                                this_types = 'interior'
                                nc.append([trs_name,*trs_info,this_types])
                            anno_res[3] = list(sp2+lar+sp1+nc)[0][-1]
                            anno_res[4] = list(sp2+lar+sp1+nc)
                result_list.append([*anno_res[:-1],readcount,*cric_info_list,anno_res[-1]])
        sorted_data = sorted(result_list, key=lambda x: x[-2], reverse=True)
        with tqdm(desc='writing to file') as pbar:
            for row in sorted_data:
                lines = '\t'.join([str(i) for i in row])
                info_lines = '$'.join([str(i) for i in row])
                chrs = row[5]
                shs1 = row[10].split('-')
                shs2 = row[11].split('-')
                pos1 = int(row[6])
                pos2 = int(row[7])
                if pos1 > pos2:
                    shs1, shs2 = shs2, shs1
                else:
                    shs1, shs2 = shs1, shs2
                pbar.update(1)
                try:
                    circout_line = f"{chrs}\t{int(shs1[0])-5}\t{int(shs1[1])+5}\t{info_lines}\n{chrs}\t{int(shs2[0])-5}\t{int(shs2[1])+5}\t{info_lines}\n"
                    circout.write(f"{circout_line}")
                    outfile.write(f"{lines}\n")
                except ValueError:
                    print(row)
                    continue

def main(project_dir, seq_dir, star_ref, genome_fa, bt2_dna_ref, gtf_dir,threads = 150,length_seq = 150,strand_specific = 'F2R1'):
    print("start********************")
    # todo check below
    kmis = length_seq//100   
    anchor_len = length_seq//5  
    if anchor_len > 40:  
        anchor_len = 40
    kmap = "1" 
    cut_off_count = 2
    min_score = str(kmis*-6)
    if not strand_specific in ['F2R1','F1R2','']:
        print("strand_specific should be F2R1 or F1R2")
        return

    print(seq_dir)

    map_dir = project_dir+"/MappingResults"
    bam_dir = os.path.join(str(map_dir), 'Aligned.sortedByCoord.out.bam')
    mate1_dir = map_dir + "/Mate1Anchor"
    mate2_dir = map_dir + "/Mate2Anchor"


    processed_dir = project_dir+"/ProcessedResults"
    raw_res = project_dir+'/circ_raw_result.tsv'
    fusion_position_dir = processed_dir+"/fusion_pos_output.txt"
    fusion_position_checkpoint = processed_dir+"/fusion_pos.ok"
    gene_bed = processed_dir+"/gene_pos.bed"
    circ_bed = processed_dir+"/circ_pos.bed"
    output_rpkm_txt = processed_dir+"/gene_exp_rpkm.txt"
    circ_rpkm_txt = processed_dir+"/circ_exp_rpkm.txt"

    remap_dir = processed_dir+ "/remap"
    remap_reference_fa = remap_dir + "/chimeric_refer.fasta"
    remap_res_dir = remap_dir+'/remap.sam'
    remap_idx_path = remap_dir + "/chimeric_refer_index"

    counting_res = processed_dir+"/counting_output.txt"

    mid_file = processed_dir+"/classify_output_gene_cut.tsv"
    output_file = project_dir+"/circleRNA_predict_result.tsv"

    if not os.path.exists(project_dir):
        os.makedirs(project_dir)
    if not os.path.exists(map_dir):
        os.makedirs(map_dir)
    if not os.path.exists(processed_dir):
        os.makedirs(processed_dir)
    if not os.path.exists(mate1_dir):
        os.makedirs(mate1_dir)
        os.makedirs(mate2_dir)
    if not os.path.exists(remap_dir):
        os.makedirs(remap_dir)
    time_total = time.time()
    """
    ___________________________________________________
                      [0]first map
    ___________________________________________________
    """
    star_map(seq_dir, star_ref, map_dir, "", out_unmap_fastx=True, threads=150 if threads > 150 else threads,alo_mis_match=False,max_out_filter_mismatch=kmis)
    """
    ___________________________________________________
            [1]anchor split 
    ___________________________________________________
    """
    unmap_path = [map_dir+"/Unmapped.out.mate1",map_dir+"/Unmapped.out.mate2"]
    anch_path = [[map_dir+"/mate1_anchor1.fa",map_dir+"/mate1_anchor2.fa"],[map_dir+"/mate2_anchor1.fa",map_dir+"/mate2_anchor2.fa"]]
    out_anch_sam = [mate1_dir+"/mate1_anchor.sam",mate2_dir+"/mate2_anchor.sam"]
    is_single = False
    if not os.path.exists(map_dir+"/Unmapped.out.mate2"):
        print("Unmapped.out.mate2 not exist")
        unmap_path = unmap_path[:1]
        anch_path = anch_path[:1]
        out_anch_sam = out_anch_sam[:1]
        is_single = True
    """
    ___________________________________________________
       [2]second map & build chimeric reference
    ___________________________________________________
    """

    if (not os.path.exists(fusion_position_dir)) or (os.path.getsize(fusion_position_dir) == 0):
        for unmap_mate_dir,sam_dir, anch_dir in zip(unmap_path,out_anch_sam,anch_path):
            cmd1 = f"seqkit fq2fa {unmap_mate_dir} | seqkit subseq -r 1:{str(anchor_len)} > {anch_dir[0]}"
            cmd2 = f"seqkit fq2fa {unmap_mate_dir} | seqkit subseq -r -{str(anchor_len)}:-1 > {anch_dir[1]}"
            print("start cut unmapped reads to get anchor")
            time1 = time.time()
            if not os.path.exists(anch_dir[0]) or os.path.getsize(anch_dir[0]) == 0 :
                os.system(cmd1)
            if not os.path.exists(anch_dir[1]) or os.path.getsize(anch_dir[1]) == 0 :
                os.system(cmd2)
            time2 = time.time()
            print("finish cut unmapped reads,cost:",arrow.get(time2-time1).format("HH:mm:ss"))

            bowtie_secondmap = f"bowtie2 -f -p {threads} -x {bt2_dna_ref} -k {kmap} -1 {anch_dir[0]} -2 {anch_dir[1]} -S {sam_dir} --score-min C,-6,0 > {sam_dir}.log 2>&1"

            if not os.path.exists(sam_dir) or os.path.getsize(sam_dir) == 0:
                print(f"start second map {sam_dir} to the DNA reference")
                time1 = time.time()
                os.system(bowtie_secondmap)
                time2 = time.time()
                print("finish first map left half to DNA,cost:", arrow.get(time2 - time1).format("HH:mm:ss"))
            else:
                print("Existing second map left half has been read at",sam_dir)
            print("parsing both sides sam file")
            left_mapped_reads,right_mapped_reads = parse_sam_paired(sam_dir)
            mates = 1 if "mate1" in sam_dir else 2 if "mate2" in sam_dir else 0
            combined_info = get_combined_info(left_mapped_reads, right_mapped_reads)
            process_combined_pair(fusion_position_dir, combined_info, unmap_mate_dir, genome_fa, kmis, anchor_len,mates,strand_specific)

    else:
        print("Existing fusion position has been read at",fusion_position_dir)
        
    if not os.path.exists(remap_reference_fa):
        finial_list =  remove_duplicate_sublists(fusion_position_dir,True)
        remap_ref(remap_reference_fa, genome_fa, finial_list, length_seq)
    else:
        print("Existing remap reference has been read at",remap_reference_fa)


    """
    ___________________________________________________
          [3]remap reads to junction site
    ___________________________________________________
    """

    if not os.path.exists(remap_idx_path+'.1.bt2'):
        print("start building remap idx...")
        time_1 = time.time()
        os.system("bowtie2-build --threads 100 "+remap_reference_fa + " " + remap_idx_path+" > /dev/null 2>&1")
        time_2 = time.time()
        print("remap bulid refer time cost:",arrow.get(time_2-time_1).format("HH:mm:ss"))
    else:
        print("Existing remap index has been read at",remap_idx_path)

    if is_single:
        bt2_str = f"bowtie2 -p {threads} --seed 42 -x {remap_idx_path} -U {unmap_path[0]} -S {remap_res_dir} --score-min C,{min_score},0 --rfg 6,6 > {remap_res_dir}.log 2>&1"
    else:
        bt2_str = f"bowtie2 -p {threads} --seed 42 -x {remap_idx_path} -1 {unmap_path[0]} -2 {unmap_path[1]} -S {remap_res_dir} --score-min C,{min_score},0 --rfg 6,6 > {remap_res_dir}.log 2>&1"

    if not os.path.exists(remap_res_dir) or os.path.getsize(remap_res_dir) == 0:
        print("start remap...")
        time1 = time.time()
        os.system(bt2_str)
        time2 = time.time()
        print("remap time cost:",arrow.get(time2-time1).format("HH:mm:ss"))
    else:
        print("Existing remap result has been read at",remap_res_dir)
    if (not os.path.exists(counting_res)) or (os.path.getsize(counting_res) == 0):
        count_alignments(remap_res_dir,counting_res,cut_off_count)
        global_time2 = time.time()
        print("total running time cost:", arrow.get(global_time2 - global_time1).format("HH:mm:ss"))
    else:
        print("Existing counting result has been read at",counting_res)


    """
    ___________________________________________________
          [4]classify circRNA and counting gene exp
    ___________________________________________________
    """

    if (not os.path.exists(mid_file) or os.path.getsize(mid_file) == 0) or (not os.path.exists(gene_bed) or os.path.getsize(gene_bed) == 0):
        annotate_circRNA(counting_res, gtf_dir, gene_bed, cut_off_count, mid_file,circ_bed)
    else:
        print("Existing classify result has been read at",mid_file)

    if os.path.exists(output_rpkm_txt) and os.path.getsize(output_rpkm_txt) > 0:
        print("Existing rpkm result has been read at",output_rpkm_txt)
    else:
        multi_rpkm(gene_bed, bam_dir, output_rpkm_txt, threads)

    if os.path.exists(circ_rpkm_txt) and os.path.getsize(circ_rpkm_txt) > 0:
        print("Existing rpkm result has been read at",circ_rpkm_txt)
    else:
        multi_rpkm(circ_bed, bam_dir, circ_rpkm_txt, threads)
    """
    ___________________________________________________
          [5] process classify result to excel
    ___________________________________________________
    """

    heads = ['genename', 'strand(circ/gene)' , 'genge_type', 'circ_type','circ_exp','chr', 'pos1', 'pos2', 'circ_strands','shsseq','shs_bod1','shs_bod2','ref_seq','other_transcripts']

    c_info = pd.read_csv(circ_rpkm_txt,sep='\t',names=['info','cov1','cov2'])

    res_dic = dict()
    filter_list = []
    with tqdm(desc='filtering',total=len(c_info)) as pbar:
        for idx, row in c_info.iterrows():
            if row['info'].startswith('total_reads'):
                continue
            try:
                if row['info'] not in res_dic:
                    res_dic[row['info']] = [int(row['cov1'])]
                else:
                    res_dic[row['info']].append(int(row['cov1']))
                    filter_list.append(row['info']+"$"+"$".join([str(i) for i in res_dic[row['info']]]))
            except ValueError:
                print(f'except row:{row["info"]}')
            pbar.update(1)
    df = pd.DataFrame(filter_list,columns=['info'],dtype=str)['info'].str.split('$',expand=True).iloc[:, :len(heads + ['cov1','cov2'])]
    df.columns = heads + ['cov1','cov2']
    df.to_csv(raw_res,sep='\t',index=False)
    df = df[df['cov1'].astype(int) > 0]
    df = df[df['cov2'].astype(int) > 0]
    df['circ_exp'].astype(int)
    df = df.sort_values(by='circ_exp',ascending=False)
    df['pos1'],df['pos2'] = df['pos1'].astype(int),df['pos2'].astype(int)
    filtered_df = pd.DataFrame()
    filtered_df['CircID'] = df['genename'].astype(str) + '/'+ df['chr'].astype(str) + ':' + (df['shs_bod1'] ).astype(str) + '-' + (df['shs_bod2'] ).astype(str) + '-' + df['strand(circ/gene)'] # // 4
    filtered_df['Gene'] = df['genename']
    filtered_df['GeneType'] = df['genge_type']
    filtered_df['Strand(circ/gene)'] = df['strand(circ/gene)']
    filtered_df['ReadCount'] = df['circ_exp'].astype(int)
    filtered_df['CircType'] = df['circ_type']
    filtered_df['CircPos'] = df['chr'].astype(str) + ':' + (df['pos1']).astype(str) + '-' + (df['pos2']).astype(str)
    filtered_df['CircLen'] = abs(df['pos1'].astype(int) - df['pos2'].astype(int))
    filtered_df['ShsSeq'] = df['shsseq']
    filtered_df['LenShs'] = filtered_df['ShsSeq'].apply(len)
    filtered_df.loc[filtered_df['ShsSeq'] == 'Na', 'LenShs'] = 0
    filtered_df['ShsArea1'] = df['chr'].astype(str) + ':' +df['shs_bod1']
    filtered_df['ShsArea2'] = df['chr'].astype(str) + ':' +df['shs_bod2']
    filtered_df["SplicOnOtherTrans"] = df['other_transcripts'].str[1:-1]
    filtered_df["Reference"] = df["ref_seq"]

    filtered_df['ReadCount'] = filtered_df['ReadCount'].astype(int)
    filtered_df = filtered_df.sort_values(by='ReadCount', ascending=False)

    filtered_df['chr'] = df['chr'].astype(str)
    filtered_df['left boundary of shs area1'] = df['shs_bod1'].apply(lambda x: int(x.split('-')[0]))
    filtered_df['right boundary of shs area1'] = df['shs_bod1'].apply(lambda x: int(x.split('-')[1]))
    filtered_df['left boundary of shs area2'] = df['shs_bod2'].apply(lambda x: int(x.split('-')[0]))
    filtered_df['right boundary of shs area2'] = df['shs_bod2'].apply(lambda x: int(x.split('-')[1]))
    filtered_df['brkpt1Cover'] = df['cov1']
    filtered_df['brkpt2Cover'] = df['cov2']

    filtered_df = filtered_df.drop_duplicates(subset=['left boundary of shs area1', 'right boundary of shs area1', 'left boundary of shs area2', 'right boundary of shs area2'], keep='first')
    filtered_df = filtered_df.drop_duplicates(subset=['left boundary of shs area2', 'right boundary of shs area2', 'left boundary of shs area1', 'right boundary of shs area1'], keep='first')


    filtered_df['CircLen'].astype(int)
    filtered_df = filtered_df[filtered_df['CircLen'] <= 50000]
    filtered_df['LenShs'].astype(int)
    filtered_df = filtered_df[filtered_df['LenShs'] <= 20]
    filtered_df.to_csv(output_file,sep='\t',index=False)
    print(f"CAT process complete total detected circRNA{len(filtered_df)}:", arrow.get(time.time() - time_total).format("HH:mm:ss"))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Process some directories.')
    parser.add_argument('--project_dir', type=str, required=True, help='Path to the project directory')
    parser.add_argument('--seq_dir', type=str, required=True, help='Path to the sequence directory')
    parser.add_argument( "--length", type=int, default=150, help="avg length of read")
    parser.add_argument( "--gtf", type=str, required=True, help="annotation GTF file")
    parser.add_argument( "--bt2", type=str, required=True, help="bowtie2 ref file")
    parser.add_argument( "--star", type=str, required=True, help="STAR ref file")
    parser.add_argument( "--strand", type=str, default='F2R1', required=True, help="F2R1 or F1R2")
    parser.add_argument( "--threads, type=int, default=4, required=True, help="Number of threads running cat")

    args = parser.parse_args()
    main(args.project_dir, args.seq_dir, args.star, args.bt2, args.dna, args.gtf, threads = args.threads,length_seq = args.length,strand_specific = args.strand)
