import pandas as pd
import sqlite3
import argparse
import os
import sys
import gzip

def parse_gff3_attributes(attr_str):
    """解析 GFF3 属性列"""
    attrs = {}
    for item in attr_str.strip().split(';'):
        if '=' in item:
            key, val = item.split('=', 1)
            attrs[key] = val
    return attrs

def build_database(gff_path, gnomad_path, db_path):
    print(f"Building robust database from {os.path.basename(gff_path)}...")
    
    # 准备数据容器
    data = {
        'genes': [],       # 基因主体
        'exons': [],       # 外显子
        'cds': [],         # 编码区
        'utrs': [],        # UTRs (5' and 3')
        'promoters': []    # 计算出的 TSS 位点
    }
    
    # 打开文件 (支持 .gz 或普通文本)
    open_func = gzip.open if gff_path.endswith('.gz') else open
    
    with open_func(gff_path, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            chrom = parts[0]
            if not chrom.startswith('chr'): chrom = 'chr' + chrom
            
            ftype = parts[2]
            start, end = int(parts[3]), int(parts[4])
            strand = parts[6]
            attr_str = parts[8]
            
            # 提取 gene_name 用于关联
            attrs = parse_gff3_attributes(attr_str)
            gene_name = attrs.get('gene_name', attrs.get('Name', '.'))
            
            row = (chrom, start, end, strand, gene_name)
            
            if ftype == 'gene':
                data['genes'].append(row)
                # 计算 TSS (Promoter center)
                tss = start if strand == '+' else end
                # 定义 Promoter 为 TSS 前后 500bp 核心区，用于快速索引
                data['promoters'].append((chrom, max(0, tss-500), tss+500, strand, gene_name))
                
            elif ftype == 'exon':
                data['exons'].append(row)
            elif ftype == 'CDS':
                data['cds'].append(row)
            elif ftype in ['five_prime_UTR', 'three_prime_UTR']:
                data['utrs'].append(row)

    # 存入 SQLite
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    
    print("Writing tables to SQLite...")
    
    # 定义表结构和数据
    tables = {
        'anno_gene': data['genes'],
        'anno_exon': data['exons'],
        'anno_cds': data['cds'],
        'anno_utr': data['utrs'],
        'anno_promoter': data['promoters']
    }
    
    cols = ['chrom', 'start', 'end', 'strand', 'gene_name']
    
    for table, rows in tables.items():
        if not rows:
            print(f"Warning: No data found for {table}")
            continue
        
        df = pd.DataFrame(rows, columns=cols)
        df.to_sql(table, conn, if_exists='replace', index=False)
        # 创建强力索引
        conn.execute(f"CREATE INDEX idx_{table}_loc ON {table} (chrom, start, end)")
        print(f"  -> {table}: {len(df)} records")

    # 加载 GnomAD
    if os.path.exists(gnomad_path):
        print("Loading GnomAD...")
        try:
            df_gnomad = pd.read_csv(gnomad_path, sep='\t')
            # 只要 pLI 和 oe_lof
            valid_cols = [c for c in ['gene', 'pLI', 'oe_lof'] if c in df_gnomad.columns]
            df_gnomad[valid_cols].to_sql('gene_metrics', conn, if_exists='replace', index=False)
            conn.execute("CREATE INDEX idx_metrics_gene ON gene_metrics (gene)")
        except Exception as e:
            print(f"GnomAD Error: {e}")
            
    conn.close()
    print("Database ready.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff', required=True)
    parser.add_argument('--gnomad', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    
    build_database(args.gff, args.gnomad, args.out)