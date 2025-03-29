import os
import pandas as pd
import sys


def format_motif(motif):
    motif = motif.upper()
    motif = motif.replace(' ', '')
    reversed_motif = ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}[base] for base in motif])
    # 得到最小字典序motif单元
    motif_mer = []
    for i in range(len(motif)):
        motif_mer.append(motif[i:] + motif[:i])
    for i in range(len(reversed_motif)):
        motif_mer.append(reversed_motif[i:] + reversed_motif[:i])
    motif_mer = list(set(motif_mer))
    motif = min(motif_mer) # 选择字典序最小的motif
    return motif


if __name__ == '__main__':
    input_path = sys.argv[1]
    input_file = sys.argv[2]
    output_filename = sys.argv[3]

    input_df = pd.read_csv(input_file, sep='\t', header=None)
    if len(input_df.columns) == 4:
        input_df.columns = ['chrom', 'start', 'end', 'motif']
    else:
        input_df.columns = ['chrom', 'start', 'end', 'motif', 'label']
    input_df['motif'] = input_df['motif'].apply(lambda x: x.upper())
    input_df['motif'] = input_df['motif'].apply(lambda x: x.replace(' ', ''))
    input_df['id'] = input_df['chrom'] + '_' + input_df['start'].astype(str) + '_' + input_df['end'].astype(str) + '_' + input_df['motif']

    input_df['motif_len'] = input_df['motif'].apply(lambda x: len(x))
    input_df['format_motif'] = input_df['motif'].apply(format_motif)
    input_df['g_value'] = input_df['format_motif'].apply(lambda x: x.count('G') + x.count('N'))
    input_df['c_value'] = input_df['format_motif'].apply(lambda x: x.count('C') + x.count('N'))
    input_df['t_value'] = input_df['format_motif'].apply(lambda x: x.count('T') + x.count('N'))
    input_df['a_value'] = input_df['format_motif'].apply(lambda x: x.count('A') + x.count('N'))
    input_df.drop(['format_motif'], axis=1, inplace=True)
    input_df['trf_len'] = abs(input_df['end'] - input_df['start'])

    input_df.drop(['chrom', 'start', 'end', 'motif'], axis=1, inplace=True)
    
    geneanno_df = pd.read_csv(input_path + '/annotation_feature.csv', sep='\t')
    geneanno_df['id'] = geneanno_df['chrom'] + '_' + geneanno_df['start'].astype(str) + '_' + geneanno_df['end'].astype(str) + '_' + geneanno_df['motif']
    geneanno_df.drop(['chrom', 'start', 'end', 'motif'], axis=1, inplace=True)
    assert len(geneanno_df) == len(input_df)
    # 修改时间 2025-2-23
    # 两个去重
    geneanno_df = geneanno_df.drop_duplicates(subset=['id'])
    input_df = input_df.drop_duplicates(subset=['id'])

    merged_df = pd.merge(input_df, geneanno_df, on='id', how='left')
    assert len(merged_df) == len(input_df)

    # orf_df = pd.read_csv(input_path + '/orf.bed', sep='\t')
    # #  修改时间 2025-2-23
    # # 去重
    # orf_df = orf_df.drop_duplicates(subset=['chrom', 'start', 'end', 'motif'])
    
    # # orf_df.columns = ['chrom', 'start', 'end', 'motif', 'orf_len_100', 'orf_len_reverse_100', 'orf_number_100', 'orf_len_1000', 'orf_len_reverse_1000', 'orf_number_1000']
    # orf_df['id'] = orf_df['chrom'] + '_' + orf_df['start'].astype(str) + '_' + orf_df['end'].astype(str) + '_' + orf_df['motif']
    # orf_df.drop(['chrom', 'start', 'end', 'motif'], axis=1, inplace=True)
    # assert len(orf_df) == len(merged_df)
    # merged_df_temp = pd.merge(merged_df, orf_df, on='id', how='left')
    # assert len(merged_df_temp) == len(merged_df)
    # merged_df = merged_df_temp

    gc_df = pd.read_csv(input_path + '/gc.bed', sep='\t')
    #  修改时间 2025-2-23
    # 去重
    gc_df = gc_df.drop_duplicates(subset=['chrom', 'start', 'end', 'motif'])
    # gc_df.columns = ['chrom', 'start', 'end', 'motif', 'gc_1000', 'at_1000', 'gc_10000', 'at_10000']
    gc_df['id'] = gc_df['chrom'] + '_' + gc_df['start'].astype(str) + '_' + gc_df['end'].astype(str) + '_' + gc_df['motif']
    gc_df.drop(['chrom', 'start', 'end', 'motif'], axis=1, inplace=True)
    assert len(gc_df) == len(merged_df)
    merged_df_temp = pd.merge(merged_df, gc_df, on='id', how='left')
    assert len(merged_df_temp) == len(merged_df)
    merged_df = merged_df_temp

    str_content_df = pd.read_csv(input_path + '/str_content.bed', sep='\t')
    #  修改时间 2025-2-23
    # 去重
    str_content_df = str_content_df.drop_duplicates(subset=['chrom', 'start', 'end', 'motif'])
    
    # str_content_df.columns = ['chrom', 'start', 'end', 'motif', 'strc_1000', 'strc_10000']
    str_content_df['id'] = str_content_df['chrom'] + '_' + str_content_df['start'].astype(str) + '_' + str_content_df['end'].astype(str) + '_' + str_content_df['motif']
    str_content_df.drop(['chrom', 'start', 'end', 'motif'], axis=1, inplace=True)
    assert len(str_content_df) == len(merged_df)
    merged_df_temp = pd.merge(merged_df, str_content_df, on='id', how='left')
    assert len(merged_df_temp) == len(merged_df)
    merged_df = merged_df_temp

    merged_df.to_csv(output_filename + '.csv', sep='\t', index=False)

    # 将除label，id外的列归一化
    columns = merged_df.columns
    for column in columns:
        if column not in ['label', 'id']:
            if merged_df[column].max() == merged_df[column].min():
                merged_df[column] = 0
            elif merged_df[column].max() == 0 or merged_df[column].min() == 0:
                merged_df[column] = 0
            else:
                merged_df[column] = (merged_df[column] - merged_df[column].min()) / (merged_df[column].max() - merged_df[column].min())
                # 6位小数
                merged_df[column] = merged_df[column].apply(lambda x: round(x, 6))
    merged_df.to_csv(output_filename + ".normalized", sep='\t', index=False)

