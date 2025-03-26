import sys
import pandas as pd


if __name__ == '__main__':
    input_file = sys.argv[1]
    output_bed = sys.argv[2]
    fw = open(output_bed, 'w')
    seq1000_dict = {}
    seq10000_dict = {}
    start_flag = False
    repeats_array = [0 for i in range(10001)]
    name = ""
    with open(input_file, 'r') as rf:
        for line in rf:
            if line.startswith('Sequence:'):
                if start_flag:
                    repeats_len = sum(repeats_array)
                    name = name.replace(" ", "")
                    chrom, start, end, motif, seq_type = name.split('_')
                    keys = f'{chrom}_{start}_{end}_{motif}'
                    if seq_type == 'seq1000':
                        seq1000_dict[keys] = repeats_len
                    else:
                        seq10000_dict[keys] = repeats_len
                    repeats_array = [0 for i in range(10001)]
                start_flag = True
                name = line.strip().split(':')[1]
                
            elif line.startswith('Parameters'):
                continue
            elif line == '\n':
                continue
            else:
                line_info = line.strip().split()
                start = int(line_info[0])
                end = int(line_info[1])
                for i in range(start, end):
                    repeats_array[i] = 1
        repeats_len = sum(repeats_array)
        name = name.replace(" ", "")
        chrom, start, end, motif, seq_type = name.split('_')
        keys = f'{chrom}_{start}_{end}_{motif}'
        if seq_type == 'seq1000':
            seq1000_dict[keys] = repeats_len
        else:
            seq10000_dict[keys] = repeats_len
    merges = {k: [seq1000_dict[k], seq10000_dict[k]] for k in seq1000_dict}
    df = pd.DataFrame(merges).T
    df.columns = ['strc_1000', 'strc_10000']
    df['chrom'] = df.index.str.split('_').str[0]
    df['start'] = df.index.str.split('_').str[1]
    df['end'] = df.index.str.split('_').str[2]
    df['motif'] = df.index.str.split('_').str[3]
    df = df[['chrom', 'start', 'end', 'motif', 'strc_1000', 'strc_10000']]
    df.to_csv(output_bed, sep='\t', index=False)
