import os
import subprocess


def get_feature(script_dir, target_file, db_file, fasta_file, out_dir, outname, logger, reference_name, resume, threads=20):
    # 获取绝对路径
    script_dir = os.path.abspath(script_dir)

    out_dir = os.path.abspath(out_dir)

    temp_dir = os.path.join(out_dir, 'temp_' + outname)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir, exist_ok=True)

    logger.debug(f"当前目录: {os.getcwd()}")
    
    gene_annotation_script = os.path.join(script_dir, 'get_gene_annotation.py')
    if resume and os.path.exists(os.path.join(temp_dir, 'annotation_feature.csv')):
        logger.info("基因注释特征提取完成 (跳过)")
    else:
        command = ['python', gene_annotation_script, target_file, db_file, temp_dir, str(threads)]
        logger.debug(f"{' '.join(command)}")
        subprocess.run(command, check=True)
        logger.info("基因注释特征提取完成")
    # Run feature extraction scripts
    feature_scripts = [
        ("get_seq_feature.py", "gc.bed", "GC含量特征提取完成"),
        ("get_flank_seq.py", "trf.fasta", "序列提取完成, 运行TRF")
    ]
    for script, output_name, message in feature_scripts:
        if resume and os.path.exists(os.path.join(temp_dir, output_name)):
            logger.info(f"{message} (跳过)")
            continue
        else:
            logger.debug(f"{' '.join(command)}")
            command = ['python', os.path.join(script_dir, script), target_file, os.path.join(temp_dir, output_name), fasta_file, reference_name, str(threads)]
            subprocess.run(command, check=True)
            logger.info(message)

    # work_dir = os.getcwd()

    # Run TRF
    # os.chdir(script_dir)
    # logger.debug(f"当前目录: {os.getcwd()}")
    # print(os.getcwd())
    # Run TRF
    # trf_out = os.path.join(temp_dir, 'trf.dat')
    # linux
    trf_path = os.path.join(script_dir, 'trf-4.10.0')
    # maxos
    # trf_path = os.path.join(script_dir, 'trf409.macosx')
    run_trf_command = ['python', os.path.join(script_dir, 'run_trf.py'), os.path.join(temp_dir, 'trf.fasta'), trf_path, temp_dir, str(threads)]
    if resume and os.path.exists(os.path.join(temp_dir, 'str_content.bed')):
        logger.info("TRF运行完成 (跳过)")
    else:
        logger.debug(f"{' '.join(run_trf_command)}")
        subprocess.run(run_trf_command, capture_output=True, text=True)
        logger.info("TRF运行完成")
    # trf_command = [trf_path, os.path.join(temp_dir, 'trf.fasta'), '2', '5', '7', '80', '10', '30', '2000', '-l', '2', '-h', '2> /dev/null']
    # logger.debug(f"{' '.join(trf_command)}")
    # subprocess.run(trf_command, capture_output=True, text=True)
    
    # mv_command = ['mv', os.path.join(script_dir, 'trf.fasta.2.5.7.80.10.30.2000.dat'), os.path.join(temp_dir, 'trf.fasta.2.5.7.80.10.30.2000.dat')]
    # logger.debug(f"{' '.join(mv_command)}")
    # subprocess.run(mv_command, check=True)

    # # subprocess.run(f"{' '.join(trf_command)}", shell=True, check=True)
    # os.chdir(work_dir)
    # logger.debug(f"当前目录: {os.getcwd()}")
    # # Clean up TRF output
    # tail_command = ['tail', '-n', '+9', os.path.join(temp_dir, 'trf.fasta.2.5.7.80.10.30.2000.dat'), '>', os.path.join(temp_dir, 'temp_file'), '&&', 'mv', os.path.join(temp_dir, 'temp_file'), os.path.join(temp_dir, 'trf.fasta.2.5.7.80.10.30.2000.dat')]
    # logger.debug(f"{' '.join(tail_command)}")
    # subprocess.run(f"{' '.join(tail_command)}", shell=True, check=True)
    
    

    logger.debug(f"当前目录: {os.getcwd()}")
    # Run final feature extraction
    tf_command = ['python', os.path.join(script_dir, 'get_strc_feature.py'), os.path.join(temp_dir, 'tmp_splited_seqs'), os.path.join(temp_dir, 'str_content.bed')]
    if resume and os.path.exists(os.path.join(temp_dir, 'strc_content.bed')):
        logger.info("侧翼序列串联重复含量提取完成 (跳过)")
    else:
        logger.debug(f"{' '.join(tf_command)}")
        subprocess.run(tf_command, check=True)
        logger.info("侧翼序列串联重复含量提取完成")
    # Merge features
    merge_command = ['python', os.path.join(script_dir, 'merge_feature.py'), temp_dir, target_file, os.path.join(out_dir, outname)]
    if resume and os.path.exists(os.path.join(out_dir, outname + '.csv')):
        logger.info("特征合并完成 (跳过)")
    else:
        logger.debug(f"{' '.join(merge_command)}")
        subprocess.run(merge_command, check=True)
        logger.info("特征合并完成")


if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__name__)
    target_file = "/home/xuminghua/STRP3/code/strp3_2_hg19/data/rexprt_input.sorted.bed"
    db_file = "/home/xuminghua/STRP3/hg19_genome/genome_annotation.db"
    fasta_file = '/home/xuminghua/STRP3/code/strp3_2_hg19/data/hs37d5.fa'
    out_dir = '/home/xuminghua/STRP3/code/strp3_2_hg19/feature0307'
    reference_name = "hg19"
    script_dir = "/home/xuminghua/STRP3/code/strp3_2_hg19/feature_tools"
    outname = "strp3_eval"
    threads=20

    get_feature(script_dir, target_file, db_file, fasta_file, out_dir, outname, logger, reference_name, resume=False, threads=threads)
