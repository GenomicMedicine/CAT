import multiprocessing
import os
import time
from tqdm import tqdm
import subprocess


def run_in_parallel(func_list, args_list, num_processes):
    pool = multiprocessing.get_context("spawn").Pool(processes=num_processes)
    results = []
    with tqdm(total=len(args_list)) as pbar:
        pbar.set_description("submitting processes")
        for idx, (func, args) in enumerate(zip(func_list, args_list)):
            result = pool.apply_async(func, args=args)
            # result = pool.map_async(func, args)
            results.append(result)
            pbar.update()
    pool.close()
    pool.join()
    output = []
    with tqdm(total=len(results)) as pbar:
        pbar.set_description("Collecting results processes")
        for result in results:
            output.append(result.get())
            pbar.update()
    return output

def calculate_rpkm(bed_file, bam_file, output_rpkm_txt, calc_total_reads=False):
    if calc_total_reads:
        output = subprocess.getoutput(f"samtools idxstats {bam_file}")
        total_reads = 0
        for line in output.splitlines():
            fields = line.split('\t')
            if len(fields) >= 3:
                read_count = int(fields[2])
                total_reads += read_count
    else:
        total_reads=1000000000
    time1 = time.time()
    with open(output_rpkm_txt, 'w') as f:
        lines = subprocess.getoutput(f"bedtools multicov -bams {bam_file} -bed {bed_file}").splitlines()
        if calc_total_reads:
            f.write(f"total_reads\t{total_reads}\t0\n")
        for line in lines:
            fields = line.split('\t')
            length = int(fields[2]) - int(fields[1])
            if length < 1:
                f.write(f"{fields[3]}\t{fields[4]}\t0\n")
            else:
                rpkm = (1000000000 * int(fields[4])) / (length * total_reads)
                f.write(f"{fields[3]}\t{fields[4]}\t{rpkm:.2f}\n")
