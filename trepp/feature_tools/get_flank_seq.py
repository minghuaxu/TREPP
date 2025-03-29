import os
import sys
from read_fasta import read_fasta
from multiprocessing import Pool

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

def process_bed_line(line, fasta, reference_name):
    line = line.strip().split('\t')
    contig = line[0]
    query_contig = contig
    if reference_name == "hg19":
        query_contig = contig.replace('chr', '')
    start = int(line[1])
    end = int(line[2])
    motif = line[3]
    motif = motif.upper()
    motif = motif.replace(' ', '')
    
    # 获取序列片段
    seq_1000 = fasta[1][query_contig][start-500:start+500]
    seq_10000 = fasta[1][query_contig][start-5000:start+5000]
    
    # 构建输出数据
    result = []
    result.append('>' + contig + '_' + str(start) + '_' + str(end) + '_' + motif + '_seq1000\n' + seq_1000)
    result.append('>' + contig + '_' + str(start) + '_' + str(end) + '_' + motif + '_seq10000\n' + seq_10000)
    
    return result

def write_output(output_bed, results):
    with open(output_bed, 'w') as fw:
        for result in results:
            for line in result:
                fw.write(line + '\n')

def main(input_bed, output_bed, fasta_file, reference_name, threads):
    fasta = read_fasta(fasta_file)
    
    with open(input_bed, 'r') as fr:
        lines = fr.readlines()

    # 使用Pool进行并行处理
    with Pool(processes=min(threads, os.cpu_count())) as pool:
        results = pool.starmap(process_bed_line, [(line, fasta, reference_name) for line in lines])
    
    # 写入输出文件
    write_output(output_bed, results)

if __name__ == '__main__':
    input_bed = sys.argv[1]
    output_bed = sys.argv[2]
    fasta_file = sys.argv[3]
    reference_name = sys.argv[4]
    threads = int(sys.argv[5])
    
    main(input_bed, output_bed, fasta_file, reference_name, threads)
