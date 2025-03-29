import re
import sys
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import os

# 解析单条序列的重复区长度
def parse_sequence(sequence_block):
    seq_id_pattern = r"(chr[^\s]+_[\d]+_[\d]+_[A-Za-z]+_seq\d+)"
    seq_id_match = re.search(seq_id_pattern, sequence_block)
    if not seq_id_match:
        return None
    seq_id = seq_id_match.group(1)
    chrom, start, end, motif, group_tag = seq_id.split('_')
    seq_id = f"{chrom}_{start}_{end}_{motif}"
    repeat_pattern = r"^\s*(\d+)\s+(\d+)"
    repeats = re.findall(repeat_pattern, sequence_block, re.MULTILINE) # MULTILINE：多行模式
    # 计算每个重复区的长度
    repeat_lengths = []
    for start, end in repeats:
        start_pos = int(start)
        end_pos = int(end)
        repeat_length = end_pos - start_pos + 1
        repeat_lengths.append(repeat_length)
    return seq_id, sum(repeat_lengths), group_tag

# 解析文件内容
def parse_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    sequence_blocks = re.split(r"(Sequence:)", content)[1:]  # split at "Sequence:"并去掉首个空块
    sequence_data = []
    for i in range(0, len(sequence_blocks), 2):
        sequence_block = sequence_blocks[i + 1]
        sequence_block = re.sub(r"Parameters:[^\n]*", "", sequence_block) # 排除 Parameters 行，并获取重复区数据
        result = parse_sequence(sequence_block)
        if result:
            sequence_data.append(result)
    return sequence_data

def process_all_files(input_dir):
    # 获取所有.dat文件路径
    dat_files = [f for f in Path(input_dir).glob("*.dat")]

    # 使用线程池并行处理所有文件
    with ThreadPoolExecutor(max_workers=min(len(dat_files), os.cpu_count())) as executor:
        results = list(executor.map(parse_file, dat_files))

    # 汇总所有线程的结果
    all_data = []
    for file_data in results:
        print(len(file_data))
        all_data.extend(file_data)

    # 将结果转化为DataFrame
    df = pd.DataFrame(all_data, columns=['seq_id', 'strc', 'group_tag'])

    df_seq1000 = df[df['group_tag'] == 'seq1000']
    df_seq1000 = df_seq1000.drop(columns=['group_tag'])

    df_seq10000 = df[df['group_tag'] == 'seq10000']
    df_seq10000 = df_seq10000.drop(columns=['group_tag'])

    df_merged = pd.merge(df_seq1000, df_seq10000, on='seq_id', suffixes=('_1000', '_10000'))
    df_merged['chrom'] = df_merged['seq_id'].apply(lambda x: x.split('_')[0])
    df_merged['start'] = df_merged['seq_id'].apply(lambda x: int(x.split('_')[1]))
    df_merged['end'] = df_merged['seq_id'].apply(lambda x: int(x.split('_')[2]))
    df_merged['motif'] = df_merged['seq_id'].apply(lambda x: x.split('_')[3])
    df_merged = df_merged[['chrom', 'start', 'end', 'motif', 'strc_1000', 'strc_10000']]

    return df_merged


# 主函数
if __name__ == "__main__":
    input_dir = sys.argv[1] # .dat文件路径
    output_bed = sys.argv[2]
    df = process_all_files(input_dir)
    df.to_csv(output_bed, sep='\t', index=False)
