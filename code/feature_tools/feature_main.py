import os
import subprocess


def get_feature(script_dir, target_file, db_file, fasta_file, out_dir, outname, logger, threads=20):
    # 获取绝对路径
    script_dir = os.path.abspath(script_dir)

    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir, exist_ok=True)

    logger.debug(f"当前目录: {os.getcwd()}")
    
    gene_annotation_script = os.path.join(script_dir, 'get_gene_annotation.py')
    command = ['python', gene_annotation_script, target_file, db_file, temp_dir, str(threads)]
    logger.debug(f"{' '.join(command)}")
    subprocess.run(command, check=True)
    logger.info("基因注释特征提取完成")
    # Run feature extraction scripts
    feature_scripts = [
        ("get_seq_feature.py", "gc.bed", "GC含量特征提取完成"),
        ("get_orf_feature.py", "orf.bed", "ORF特征提取完成"),
        ("get_flask_seq.py", "trf.fasta", "序列提取完成, 运行TRF")
    ]
    for script, output_name, message in feature_scripts:
        logger.debug(f"{' '.join(command)}")
        command = ['python', os.path.join(script_dir, script), target_file, os.path.join(temp_dir, output_name), fasta_file]
        subprocess.run(command, check=True)
        logger.info(message)

    work_dir = os.getcwd()

    # Run TRF
    os.chdir(script_dir)
    logger.debug(f"当前目录: {os.getcwd()}")
    print(os.getcwd())
    # Run TRF
    # trf_out = os.path.join(temp_dir, 'trf.dat')
    trf_path = os.path.join(script_dir, 'trf409.macosx')
    trf_command = [trf_path, os.path.join(temp_dir, 'trf.fasta'), '2', '5', '7', '80', '10', '30', '2000', '-l', '2', '-h', '2> /dev/null']
    logger.debug(f"{' '.join(trf_command)}")
    subprocess.run(trf_command, capture_output=True, text=True)
    
    mv_command = ['mv', os.path.join(script_dir, 'trf.fasta.2.5.7.80.10.30.2000.dat'), os.path.join(temp_dir, 'trf.fasta.2.5.7.80.10.30.2000.dat')]
    logger.debug(f"{' '.join(mv_command)}")
    subprocess.run(mv_command, check=True)

    # subprocess.run(f"{' '.join(trf_command)}", shell=True, check=True)
    os.chdir(work_dir)
    logger.debug(f"当前目录: {os.getcwd()}")
    # Clean up TRF output
    tail_command = ['tail', '-n', '+9', os.path.join(temp_dir, 'trf.fasta.2.5.7.80.10.30.2000.dat'), '>', os.path.join(temp_dir, 'temp_file'), '&&', 'mv', os.path.join(temp_dir, 'temp_file'), os.path.join(temp_dir, 'trf.fasta.2.5.7.80.10.30.2000.dat')]
    logger.debug(f"{' '.join(tail_command)}")
    subprocess.run(f"{' '.join(tail_command)}", shell=True, check=True)
    
    logger.info("TRF运行完成")

    logger.debug(f"当前目录: {os.getcwd()}")
    # Run final feature extraction
    tf_command = ['python', os.path.join(script_dir, 'get_strc_feature.py'), os.path.join(temp_dir, 'trf.fasta.2.5.7.80.10.30.2000.dat'), os.path.join(temp_dir, 'str_content.bed')]
    logger.debug(f"{' '.join(tf_command)}")
    subprocess.run(tf_command, check=True)
    logger.info("侧翼序列串联重复含量提取完成")
    # Merge features
    merge_command = ['python', os.path.join(script_dir, 'merge_feature.py'), temp_dir, target_file, os.path.join(out_dir, outname)]
    logger.debug(f"{' '.join(merge_command)}")
    subprocess.run(merge_command, check=True)
    logger.info("特征合并完成")


if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__name__)
    target_file = "/public/home/hpc214707002/code_space/strp3/data/train.new.sorted.bed"
    db_file = "/public/home/hpc214707002/code_space/strp3_2/feature_tools/genome_annotation.db"
    fasta_file = '/public/home/hpc214707002/data/GRCH38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
    out_dir = '/public/home/hpc214707002/code_space/strp3_2/out'
    script_dir = "/public/home/hpc214707002/code_space/strp3_2/feature_tools"
    outname = "strp3"
    threads=20

    get_feature(script_dir, target_file, db_file, fasta_file, out_dir, outname, logger, threads=threads)
