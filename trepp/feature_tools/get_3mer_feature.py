import os
import sys
import itertools
import numpy as np
from read_fasta import read_fasta
import pandas as pd

class KmerFeature:
    def __init__(self, input_bed_path, output_file, fasta_path, k=3):
        self.k = k
        self.input_bed_path = input_bed_path
        self.output_file = output_file
        self.fasta_path = fasta_path
        self.fasta = read_fasta(fasta_path)
        self.kmer_list = self._generate_3mer_list()
        self.kmer_ids_list = None
        self.str_ids = None
        self.labels = None

    def _generate_3mer_list(self):
        bases = ['A','G','C','T', '-']
        kmer_tuple_list = list(itertools.product(bases, repeat=int(self.k)))
        kmer_list = [''.join(kmer_tuple) for kmer_tuple in kmer_tuple_list]
        return kmer_list

    def _get_seq(self, contig, start, end, seq_len = 500):
        left_seq = self.fasta[1][contig][start-500:start]
        if len(left_seq) < seq_len:
            left_seq = '-'*(seq_len-len(left_seq)) + left_seq
        right_seq = self.fasta[1][contig][end:end+500]
        if len(right_seq) < seq_len:
            right_seq = right_seq + '-'*(seq_len-len(right_seq))
        return left_seq, right_seq


    def _parse_bed(self):
        left_kmer_ids_list = []
        right_kmer_ids_list = []
        str_ids = []
        with open(self.input_bed_path, 'r') as fr:
            for line in fr:
                line = line.strip().split('\t')
                contig = line[0]
                contig = contig.replace('chr', '') # remove 'chr' prefix, 只有HG19有chr
                start = int(line[1])
                end = int(line[2])
                motif = line[3]
                motif = motif.upper()
                motif = motif.replace(' ', '')
                label = line[4]
                left_seq, right_seq = self._get_seq(contig, start, end)
                left_kmer_feature = self.get_3mer_feature(left_seq)
                left_kmer_ids_list.append(left_kmer_feature)
                right_kmer_feature = self.get_3mer_feature(right_seq)
                right_kmer_ids_list.append(right_kmer_feature)
                str_ids.append(f'chr{contig}_{start}_{end}_{motif}_{label}')
        left_kmer_ids_list = np.array(left_kmer_ids_list)
        right_kmer_ids_list = np.array(right_kmer_ids_list)
        kmer_ids_list = np.hstack((left_kmer_ids_list, right_kmer_ids_list))
        return kmer_ids_list.tolist() , str_ids

    def _write_output(self):
        contigs = []
        starts = []
        ends = []
        motifs = []
        labels = []
        for str_id in self.str_ids:
            contig, start, end, motif, label = str_id.split('_')
            contigs.append(contig)
            starts.append(start)
            ends.append(end)
            motifs.append(motif)
            labels.append(label)
        self.labels = labels
        df = pd.DataFrame({'chrom': contigs, 'start': starts, 'end': ends, 'motif': motifs, 'label': labels, '3mer_ids': self.kmer_ids_list})
        df.to_csv(self.output_file, sep='\t', index=False)

    def get_3mer_feature(self, seq):
        seq_3mers = [seq[i:i+3] for i in range(len(seq)-3+1)]
        seq_3mer_ids = [self.kmer_list.index(mer) for mer in seq_3mers]
        return seq_3mer_ids

    def main(self):
        self.kmer_ids_list, self.str_ids = self._parse_bed()
        self._write_output()

if __name__ == '__main__':
    input_bed_path = sys.argv[1]
    output_file = sys.argv[2]
    fasta_path = sys.argv[3]

    kmer_feature = KmerFeature(input_bed_path, output_file, fasta_path)
    kmer_feature.main()

