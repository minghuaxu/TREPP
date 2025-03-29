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
        self.seq_feature_list = None
        self.str_ids = None
        self.labels = None

    def _get_seq(self, contig, start, end, seq_len = 500):
        left_seq = self.fasta[1][contig][start-500:start]
        if len(left_seq) < seq_len:
            left_seq = '-'*(seq_len-len(left_seq)) + left_seq
        right_seq = self.fasta[1][contig][end:end+500]
        if len(right_seq) < seq_len:
            right_seq = right_seq + '-'*(seq_len-len(right_seq))
        return left_seq, right_seq


    def _parse_bed(self):
        left_feature_list = []
        right_feature_list = []
        str_ids = []
        base_dict = {
            '-': 0,
            'A': 1,
            'G': 2,
            'C': 3,
            'T': 4
        }
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
                left_seq_feature = [base_dict[base] for base in left_seq]
                left_feature_list.append(left_seq_feature)
                right_seq_feature = [base_dict[base] for base in right_seq]
                right_feature_list.append(right_seq_feature)
                str_ids.append(f'chr{contig}_{start}_{end}_{motif}_{label}')
        left_feature_list = np.array(left_feature_list)
        right_feature_list = np.array(right_feature_list)
        seq_feature_list = np.hstack((left_feature_list, right_feature_list))
        return seq_feature_list.tolist() , str_ids

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
        df = pd.DataFrame({'chrom': contigs, 'start': starts, 'end': ends, 'motif': motifs, 'label': labels, '3mer_ids': self.seq_feature_list})
        df.to_csv(self.output_file, sep='\t', index=False)

    def main(self):
        self.seq_feature_list, self.str_ids = self._parse_bed()
        self._write_output()

if __name__ == '__main__':
    input_bed_path = sys.argv[1]
    output_file = sys.argv[2]
    fasta_path = sys.argv[3]

    kmer_feature = KmerFeature(input_bed_path, output_file, fasta_path)
    kmer_feature.main()

