import pandas as pd
import sys
import sqlite3
import os
from multiprocessing import Pool, Manager
from functools import partial

def parse_chrom(chrom):
        chrom = chrom.replace("chr", "")
        if chrom == "X":
            return 23
        if chrom == "Y":
            return 24
        if chrom == "M":
            return 25
        return  int(chrom)

def parse_chrom_reverse(chrom):
    if chrom == 23:
        return "chrX"
    if chrom == 24:
        return "chrY"
    if chrom == 25:
        return "chrM"
    return f"chr{chrom}"

querys = {
    "gene": "SELECT * FROM gene WHERE chromosome = ? AND ((start > ? AND start <= ?) or (end >= ? AND end < ?) or (start <= ? AND end >= ?));",
    "gene_flank5000": "SELECT * FROM gene WHERE chromosome = ? AND  ((start-5000 > ? AND start-5000 <= ?) or (end+5000 >= ? AND end+5000 < ?) or (start-5000 <= ? AND end+5000 >= ?)) LIMIT 1;",
    "transcript": "SELECT * FROM transcript WHERE chromosome = ? AND ((start > ? AND start <= ?) or (end >= ? AND end < ?) or (start <= ? AND end >= ?));",
    "trainscript_flank5000": "SELECT * FROM transcript WHERE chromosome = ? AND  ((start-5000 > ? AND start-5000 <= ?) or (end+5000 >= ? AND end+5000 < ?) or (start-5000 <= ? AND end+5000 >= ?)) LIMIT 1;",
    "UTR3": "SELECT * FROM UTR3 WHERE chromosome = ? AND ((start > ? AND start <= ?) or (end >= ? AND end < ?) or (start <= ? AND end >= ?));",
    "UTR3_flank500": "SELECT * FROM UTR3 WHERE chromosome = ? AND  ((start-500 > ? AND start-500 <= ?) or (end+500 >= ? AND end+500 < ?) or (start-500 <= ? AND end+500 >= ?)) LIMIT 1;",
    "UTR5": "SELECT * FROM UTR5 WHERE chromosome = ? AND  ((start > ? AND start <= ?) or (end >= ? AND end < ?) or (start <= ? AND end >= ?));",
    "UTR5_flank500": "SELECT * FROM UTR5 WHERE chromosome = ? AND  ((start-500 > ? AND start-500 <= ?) or (end+500 >= ? AND end+500 < ?) or (start-500 <= ? AND end+500 >= ?)) LIMIT 1;",
    "intron": "SELECT * FROM Intron WHERE chromosome = ? AND  ((start > ? AND start <= ?) or (end >= ? AND end < ?) or (start <= ? AND end >= ?));",
    "intron_flank500": "SELECT * FROM Intron WHERE chromosome = ? AND  ((start-500 > ? AND start-500 <= ?) or (end+500 >= ? AND end+500 < ?) or (start-500 <= ? AND end+500 >= ?)) LIMIT 1;",
    "exon": "SELECT * FROM exon WHERE chromosome = ? AND  ((start > ? AND start <= ?) or (end >= ? AND end < ?) or (start <= ? AND end >= ?));",
    "exon_flank500": "SELECT * FROM exon WHERE chromosome = ? AND  ((start-500 > ? AND start-500 <= ?) or (end+500 >= ? AND end+500 < ?) or (start-500 <= ? AND end+500 >= ?)) LIMIT 1;",
    "CDS": "SELECT * FROM CDS WHERE chromosome = ? AND  ((start > ? AND start <= ?) or (end >= ? AND end < ?) or (start <= ? AND end >= ?));",
    "CDS_flank500": "SELECT * FROM CDS WHERE chromosome = ? AND  ((start-500 > ? AND start-500 <= ?) or (end+500 >= ? AND end+500 < ?) or (start-500 <= ? AND end+500 >= ?)) LIMIT 1;",
    "Promoter": "SELECT * FROM promoter WHERE chromosome = ? AND  ((start > ? AND start <= ?) or (end >= ? AND end < ?) or (start <= ? AND end >= ?));",
    "Promoter_flank500": "SELECT * FROM promoter WHERE chromosome = ? AND  ((start-500 > ? AND start-500 <= ?) or (end+500 >= ? AND end+500 < ?) or (start-500 <= ? AND end+500 >= ?)) LIMIT 1;",
    "Enhancer": "SELECT * FROM enhancer WHERE chromosome = ? AND  ((start > ? AND start <= ?) or (end >= ? AND end < ?) or (start <= ? AND end >= ?));",
    "Enhancer_flank500": "SELECT * FROM enhancer WHERE chromosome = ? AND  ((start-500 > ? AND start-500 <= ?) or (end+500 >= ? AND end+500 < ?) or (start-500 <= ? AND end+500 >= ?)) LIMIT 1;",
    "CTCF": "SELECT * FROM CTCF_binding_site WHERE chromosome = ? AND  ((start > ? AND start <= ?) or (end >= ? AND end < ?) or (start <= ? AND end >= ?));",
    "CTCF_flank500": "SELECT * FROM CTCF_binding_site WHERE chromosome = ? AND  ((start-500 > ? AND start-500 <= ?) or (end+500 >= ? AND end+500 < ?) or (start-500 <= ? AND end+500 >= ?)) LIMIT 1;",
    "stop_codon": "SELECT * FROM stop_codon WHERE chromosome = ? AND  ((start > ? AND start <= ?) or (end >= ? AND end < ?) or (start <= ? AND end >= ?));",
    "stop_codon_flank100": "SELECT * FROM stop_codon WHERE chromosome = ? AND  ((start-100 > ? AND start-100 <= ?) or (end+100 >= ? AND end+100 < ?) or (start-100 <= ? AND end+100 >= ?)) LIMIT 1;",
    "start_codon": "SELECT * FROM start_codon WHERE chromosome = ? AND  ((start > ? AND start <= ?) or (end >= ? AND end < ?) or (start <= ? AND end >= ?));",
    "start_codon_flank100": "SELECT * FROM start_codon WHERE chromosome = ? AND  ((start-100 > ? AND start-100 <= ?) or (end+100 >= ? AND end+100 < ?) or (start-100 <= ? AND end+100 >= ?)) LIMIT 1;"
}

columns = ["chrom", "start", "end", "motif", "gene", "gene_flank5000", "transcript", "trainscript_flank5000", "UTR3", "UTR3_flank500", "UTR5", "UTR5_flank500", "exon", "exon_flank500", "intron", "intron_flank500", "promoter", "promoter_flank500", "enhancer", "enhancer_flank500", "CTCF", "CTCF_flank500", "stop_codon", "start_codon", "start_codon_flank100", "stop_codon_flank100"]

def query_database_for_row(cursor, row, querys):
    chrom, start, end, motif = parse_chrom(row[0]), int(row[1]), int(row[2]), row[3]

    gene, gene_flank5000, transcript, trainscript_flank5000, UTR3, UTR3_flank500, UTR5, UTR5_flank500, exon, exon_flank500, intron, intron_flank500, promoter, promoter_flank500, enhancer, enhancer_flank500, CTCF, CTCF_flank500, start_codon, start_codon_flank100,stop_codon , stop_codon_flank100 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 
    
    for key in querys:
        cursor.execute(querys[key], (chrom, start, end, start, end, start, end))
        result = cursor.fetchall()
        if result:
            if key == "gene":
                gene = 1
            elif key == "gene_flank5000":
                gene_flank5000 = 1
            elif key == "transcript":
                transcript = 1
            elif key == "trainscript_flank5000":
                trainscript_flank5000 = 1
            elif key == "UTR3":
                UTR3 = 1    
            elif key == "UTR3_flank500":
                UTR3_flank500 = 1
            elif key == "UTR5":
                UTR5 = 1
            elif key == "UTR5_flank500":
                UTR5_flank500 = 1
            elif key == "exon":
                exon = 1
            elif key == "exon_flank500":
                exon_flank500 = 1
            elif key == "intron":
                intron = 1
            elif key == "intron_flank500":
                intron_flank500 = 1
            elif key == "Promoter":
                promoter = 1
            elif key == "Promoter_flank500":
                promoter_flank500 = 1
            elif key == "Enhancer":
                enhancer = 1
            elif key == "Enhancer_flank500":
                enhancer_flank500 = 1
            elif key == "CTCF":
                CTCF = 1
            elif key == "CTCF_flank500":
                CTCF_flank500 = 1
            elif key == "stop_codon":
                stop_codon = 1
            elif key == "start_codon":
                start_codon = 1
            elif key == "start_codon_flank100":
                start_codon_flank100 = 1
            elif key == "stop_codon_flank100":
                stop_codon_flank100 = 1

    return [parse_chrom_reverse(chrom), start, end, motif, gene, gene_flank5000, transcript, trainscript_flank5000, UTR3, UTR3_flank500, UTR5, UTR5_flank500, exon, exon_flank500, intron, intron_flank500, promoter, promoter_flank500, enhancer, enhancer_flank500, CTCF, CTCF_flank500, start_codon, start_codon_flank100,stop_codon , stop_codon_flank100]

def process_chunk(chunk_df, db_file, querys):
        conn = sqlite3.connect(db_file)
        cursor = conn.cursor()
        results = []
        for _, row in chunk_df.iterrows():
            result = query_database_for_row(cursor, row.iloc, querys)
            results.append(result)
        conn.close()
        return results

def main(bed_file, db_file, output_dir, threads):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    df = pd.read_csv(bed_file, sep="\t", header=None)
    
    # Open database connection once per process

    # Set up multiprocessing
    chunk_size = 100  # Adjust based on your machine
    num_chunks = len(df) // chunk_size + 1
    chunks = [df[i*chunk_size:(i+1)*chunk_size] for i in range(num_chunks)]
    print(f"threads: {threads}, chunks: {num_chunks}")

    with Pool(processes=min(threads, os.cpu_count())) as pool:
        process_chunk_partial = partial(process_chunk, db_file=db_file, querys=querys)
        result_chunks = pool.map(process_chunk_partial, chunks)


    # Combine all results
    all_results = [item for sublist in result_chunks for item in sublist]

    df_result = pd.DataFrame(all_results, columns=columns)
    df_result.to_csv(f"{output_dir}/annotation_feature.csv", index=False, header=True, sep="\t")

if __name__ == "__main__":
    bed_file = sys.argv[1]
    db_file = sys.argv[2]
    output_dir = sys.argv[3]
    threads = int(sys.argv[4])
    main(bed_file, db_file, output_dir, threads)
