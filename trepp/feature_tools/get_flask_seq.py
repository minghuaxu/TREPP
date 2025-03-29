import os
import sys
from read_fasta import read_fasta


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


if __name__ == '__main__':
    input_bed = sys.argv[1]
    output_bed = sys.argv[2]
    fasta_file = sys.argv[3]
    fasta = read_fasta(fasta_file)

    fw = open(output_bed, 'w')
    with open(input_bed , 'r') as fr:
        for line in fr:
            line = line.strip().split('\t')
            contig = line[0]
            start = int(line[1])
            end = int(line[2])
            motif = line[3]
            seq_1000 = fasta[1][contig][start-500:start+500]
            seq_10000 = fasta[1][contig][start-5000:start+5000]
            fw.write('>'+contig+'_'+str(start)+'_'+str(end)+ '_' + motif + '_seq1000\n'+seq_1000 + '\n')
            fw.write('>'+contig+'_'+str(start)+'_'+str(end)+ '_' + motif + '_seq10000\n'+seq_10000 + '\n')
        fw.close()

