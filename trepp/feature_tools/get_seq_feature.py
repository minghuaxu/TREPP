import os
import sys
from read_fasta import read_fasta
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

def GC_content(seq):
    GC = 0
    AT = 0
    seq_len = len(seq)
    for base in seq:
        if base in ['G', 'C']:
            GC += 1
        else:
            AT += 1
    return GC, AT

def process_row(row, fasta, reference_name):
    contig = row[0]
    query_contig = contig
    if reference_name == "hg19":
        query_contig = contig.replace('chr', '')
    start = int(row[1])
    end = int(row[2])
    motif = row[3]
    motif = motif.upper()
    motif = motif.replace(' ', '')
    
    # 获取序列片段
    seq_1000 = fasta[1][query_contig][start-500:end+500]
    gc_1000, at_1000 = GC_content(seq_1000)
    seq_10000 = fasta[1][query_contig][start-5000:end+5000]
    gc_10000, at_10000 = GC_content(seq_10000)
    
    return [contig, start, end, motif, gc_1000, at_1000, gc_10000, at_10000]

def main(input_bed, output_bed, fasta_file, reference_name, threads):
    fasta = read_fasta(fasta_file)
    input_df = pd.read_csv(input_bed, sep='\t', header=None)
    # 使用ThreadPoolExecutor来并行处理每一行
    with ThreadPoolExecutor(max_workers=(min(threads, os.cpu_count()))) as executor:
        datalist = list(executor.map(lambda row: process_row(row, fasta, reference_name), input_df.itertuples(index=False, name=None)))
    
    # 将结果转换为DataFrame并保存到输出文件
    df = pd.DataFrame(datalist, columns=['chrom', 'start', 'end', 'motif', 'gc_1000', 'at_1000', 'gc_10000', 'at_10000'])
    df.to_csv(output_bed, sep='\t', index=False)

if __name__ == '__main__':
    input_bed = sys.argv[1]
    output_bed = sys.argv[2]
    fasta_file = sys.argv[3]
    reference_name = sys.argv[4]
    threads = int(sys.argv[5])
    
    main(input_bed, output_bed, fasta_file, reference_name, threads)
