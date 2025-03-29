import os
import threading
import subprocess
from Bio import SeqIO
import shutil
import sys

# 分割FASTA文件为多个子文件
def split_fasta(fasta_file, output_dir, threads):
    # 自定义方法来处理每两个连续的行
    sequences = []
    with open(fasta_file, "r") as f:
        lines = f.readlines()
        for i in range(0, len(lines), 2):
            seq_id = lines[i].strip()
            seq = lines[i + 1].strip()
            sequences.append((seq_id, seq))  # 存储序列名和序列为元组
    
    total_seqs = len(sequences)
    seqs_per_thread = total_seqs // threads
    temp_dir = os.path.join(output_dir, "tmp_splited_seqs")
    
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.makedirs(temp_dir)
    
    split_files = []
    for i in range(threads):
        start_idx = i * seqs_per_thread
        if i == threads - 1:  # 最后一部分包括剩余的序列
            end_idx = total_seqs
        else:
            end_idx = (i + 1) * seqs_per_thread
        
        split_file = os.path.join(temp_dir, f"trf_{i}.fasta")
        split_files.append(split_file)
        
        with open(split_file, "w") as f:
            for seq_id, seq in sequences[start_idx:end_idx]:
                f.write(f"{seq_id}\n{seq}\n")  # 逐行写入序列名和序列
    
    return split_files, temp_dir

# # 分割FASTA文件为多个子文件
# def split_fasta(fasta_file, output_dir, threads):
#     sequences = list(SeqIO.parse(fasta_file, "fasta"))
#     total_seqs = len(sequences)
#     seqs_per_thread = total_seqs // threads
#     temp_dir = os.path.join(output_dir, "tmp_splited_seqs")
    
#     if os.path.exists(temp_dir):
#         shutil.rmtree(temp_dir)
#     os.makedirs(temp_dir)
    
#     split_files = []
#     for i in range(threads):
#         start_idx = i * seqs_per_thread
#         if i == threads - 1:  # 最后一部分包括剩余的序列
#             end_idx = total_seqs
#         else:
#             end_idx = (i + 1) * seqs_per_thread
        
#         split_file = os.path.join(temp_dir, f"trf_{i}.fasta")
#         split_files.append(split_file)
#         SeqIO.write(sequences[start_idx:end_idx], split_file, "fasta")
    
#     return split_files, temp_dir

# 调用TRF工具处理FASTA文件
def run_trf(fasta_file, trf_path, output_dir):
    os.chdir(os.path.join(output_dir, "tmp_splited_seqs"))
    print(f"当前目录: {os.getcwd()}")
    print(f"运行TRF: {fasta_file}")
    trf_command = [
        trf_path, fasta_file, '2', '5', '7', '80', '10', '30', '2000', '-l', '2', '-h', '>', '/dev/null', '2>&1'
    ]
    print(f"{' '.join(trf_command)}")
    subprocess.run(trf_command, capture_output=True, text=True)

# 去除TRF结果文件前9行
def remove_first_n_lines(file_path, n=7):
    with open(file_path, 'r') as f:
        lines = f.readlines()[n:]
    
    with open(file_path, 'w') as f:
        f.writelines(lines)

# 合并所有结果文件
def merge_results(temp_dir, output_dir):
    all_files = [f for f in os.listdir(temp_dir) if f.endswith('.dat')]
    all_files.sort()  # 可以根据需要调整合并顺序
    filename = "merged_trf_output.dat"
    output_file = os.path.join(output_dir, filename)
    with open(output_file, 'w') as outfile:
        for filename in all_files:
            file_path = os.path.join(temp_dir, filename)
            with open(file_path, 'r') as infile:
                outfile.write(infile.read())

# 主执行逻辑
def process_fasta(fasta_file, trf_path, output_dir, threads):
    split_files, temp_dir = split_fasta(fasta_file, output_dir, threads)
    
    def thread_task(fasta_file):
        run_trf(fasta_file, trf_path, output_dir)
        remove_first_n_lines(fasta_file + '.2.5.7.80.10.30.2000.dat')
    
    threads_list = []
    for split_file in split_files:
        t = threading.Thread(target=thread_task, args=(split_file,))
        threads_list.append(t)
        t.start()
    
    for t in threads_list:
        t.join()
    merge_results(temp_dir, output_dir)

    # 清理临时文件夹
    # shutil.rmtree(temp_dir)

# 调用主函数
if __name__ == "__main__":
    fasta_file = sys.argv[1]
    trf_path = sys.argv[2]
    output_dir = sys.argv[3]
    threads = int(sys.argv[4])  # 从命令行传入线程数

    process_fasta(fasta_file, trf_path, output_dir, threads=threads)
