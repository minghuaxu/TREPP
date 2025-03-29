import os
import sys
import itertools
import numpy as np
from read_fasta import read_fasta


def generate_kmer_dic(repeat_num):
    ##initiate a dic to store the kmer dic
    ##kmer_dic = {'ATC':0,'TTC':1,...}
    kmer_dic = {}
    bases = ['A','G','C','T', '-']
    kmer_list = list(itertools.product(bases, repeat=int(repeat_num)))
    reduce_mer = '-' * repeat_num # 去掉全空的kmer
    for eachitem in kmer_list:
        each_kmer = ''.join(eachitem)
        if each_kmer == reduce_mer:
            continue
        kmer_dic[each_kmer] = 0
    return (kmer_dic)

def get_freq_feature(seq, Ks):
    freq_features = []
    for k in Ks:
        kmer_dic = generate_kmer_dic(k)
        kmer_list_len = len(seq) - k + 1
        mers = []
        for i in range(kmer_list_len):
            mers.append(seq[i:i+k])
        reduce_mer = '-' * k
        for mer in mers:
            if mer == reduce_mer:
                continue
            kmer_dic[mer] += 1
        num_list = []  ##this dic stores num_dic = [0,1,1,0,3,4,5,8,2...]
        for eachkmer in kmer_dic:
            num_list.append(kmer_dic[eachkmer])
        freq_features += num_list
    freq_features = np.array(freq_features)
    return freq_features

if __name__ == '__main__':
    input_bed = sys.argv[1]
    output_bed = sys.argv[2]
    fasta_file = sys.argv[3]
    reference_name = sys.argv[4]
    fasta = read_fasta(fasta_file)

    fw = open(output_bed, 'w')
    with open(input_bed , 'r') as fr:
        for line in fr:
            line = line.strip().split('\t')
            contig = line[0]
            start = int(line[1])
            end = int(line[2])
            motif = line[3]
            motif = motif.upper()
            motif = motif.replace(' ', '')
            seq_1000 = fasta[1][contig][start-500:end+500]
            
            Ks = [1, 2, 3, 4]
            kmer_feature = get_freq_feature(seq_1000, Ks)

            str_kmer = [str(item) for item in kmer_feature]
            fw.write('\t'.join([contig, str(start), str(end), motif, line[4], '\t'.join(str_kmer)]) + '\n')
        fw.close()

