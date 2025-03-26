import os
import sys
from read_fasta import read_fasta
import pandas as pd


def GC_content(seq):
    GC = 0
    AT = 0
    seq_len = len(seq)
    for base in seq:
        if base in ['G', 'C']:
            GC += 1
        else:
            AT += 1
    # return GC / seq_len, AT / seq_len
    return GC, AT

# 主函数
if __name__ == '__main__':
    input_bed = sys.argv[1]
    output_bed = sys.argv[2]
    fasta_file = sys.argv[3]
    fasta = read_fasta(fasta_file)
    input_df = pd.read_csv(input_bed, sep='\t', header=None)
    datalist = []

    for i in range(len(input_df)):
        row = input_df.iloc[i]
        contig = row[0]
        start = int(row[1])
        end = int(row[2])
        motif = row[3]
        seq_1000 = fasta[1][contig][start-500:end+500]
        gc_1000, at_1000 = GC_content(seq_1000)
        seq_10000 = fasta[1][contig][start-5000:end+5000]
        gc_10000, at_10000 = GC_content(seq_10000)
        datalist.append([contig, start, end, motif, gc_1000, at_1000, gc_10000, at_10000])
    df = pd.DataFrame(datalist, columns=['chrom', 'start', 'end', 'motif', 'gc_1000', 'at_1000', 'gc_10000', 'at_10000'])
    df.to_csv(output_bed, sep='\t', index=False)
