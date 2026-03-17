import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import os
from tqdm import tqdm

plt.rcParams.update({'font.size': 16})

# ================= 配置路径 =================
# RExPRT 带有分数的原始文件
REXPRT_SCORE_FILE = "reference_RExPRTscores.txt"
# TREPP 的结果
TREPP_CSV = "trepp_hg19_production.csv"

# 映射与数据库
GFF3_FILE = "gencode.v49lift37.annotation.gff3"
OMIM_FILE = "mim2gene.txt"
GENETREK_FILE = "genetrek-data-2024-04-26.tsv"
TOP_K = 10000
# 输出文件
OUTPUT_OMIM_IMG = f"venn_omim_top{TOP_K}_ranked.png"
OUTPUT_GENETREK_IMG = f"venn_genetrek_top{TOP_K}_ranked.png"


# ================= 1. 基因映射工具 (保持公平，统一使用坐标映射) =================

def parse_gff3_to_genes(gff_path):
    print(f"Loading Genes from GFF3: {gff_path} ...")
    genes_by_chr = {}
    with open(gff_path, 'r') as f:
        for line in tqdm(f, desc="Parsing GFF3"):
            if line.startswith("#"): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            if parts[2] != 'gene': continue
            
            chrom = parts[0]
            if not chrom.startswith('chr'): chrom = 'chr' + chrom
            start, end = int(parts[3]), int(parts[4])
            
            # 解析 Ensembl ID
            gene_id = None
            for attr in parts[8].split(';'):
                if attr.startswith('gene_id='):
                    gene_id = attr.split('=')[1].split('.')[0].split('_')[0] # 去版本号
                    break
            
            if gene_id:
                if chrom not in genes_by_chr: genes_by_chr[chrom] = []
                genes_by_chr[chrom].append((start, end, gene_id))
    return genes_by_chr

def map_loci_to_best_gene(loci_list, genes_db):
    mapped = set()
    for chrom, l_start, l_end in loci_list:
        if chrom not in genes_db:
            continue
        best_gene = None
        best_olap = 0
        for g_start, g_end, g_id in genes_db[chrom]:
            olap = min(l_end, g_end) - max(l_start, g_start)
            if olap > best_olap:   # 只接受真正 overlap（>0）
                best_olap = olap
                best_gene = g_id
        if best_gene is not None:
            mapped.add(best_gene)
    return mapped


def map_loci_to_genes(loci_list, genes_db):
    """
    loci_list: [(chr, start, end), ...]
    返回: set(gene_ids)
    """
    mapped_genes = set()
    loci_by_chr = {}
    for chrom, start, end in loci_list:
        if chrom not in loci_by_chr: loci_by_chr[chrom] = []
        loci_by_chr[chrom].append((start, end))
        
    for chrom, loci in loci_by_chr.items():
        if chrom not in genes_db: continue
        genes = genes_db[chrom]
        # 简单的重叠检查
        for g_start, g_end, g_id in genes:
            for l_start, l_end in loci:
                if max(l_start, g_start) < min(l_end, g_end):
                    mapped_genes.add(g_id)
    return mapped_genes

# ================= 2. 数据加载与排序 =================

def get_trepp_top_k(csv_path, genes_db, k):
    print(f"\nProcessing TREPP (Top {k})...")
    # 1. 读取
    try:
        df = pd.read_csv(csv_path, sep='\t')
        if len(df.columns) < 2: df = pd.read_csv(csv_path, sep=',')
    except:
        df = pd.read_csv(csv_path, delim_whitespace=True)
    
    # 2. 排序 (按 prob 降序)
    df_sorted = df.sort_values(by='prob', ascending=False)
    
    # 3. 取前 K
    df_top = df_sorted.head(k)
    print(f"  -> TREPP Score Range: {df_top['prob'].max():.4f} to {df_top['prob'].min():.4f}")
    
    # 4. 提取坐标
    loci = []
    for _, row in df_top.iterrows():
        if isinstance(row['id'], str):
            p = row['id'].split('_')
            if len(p) >= 3: loci.append((p[0], int(p[1]), int(p[2])))
            
    # 5. 映射基因
    # genes = map_loci_to_genes(loci, genes_db)
    genes = map_loci_to_best_gene(loci, genes_db)

    print(f"  -> Mapped to {len(genes)} genes.")
    return genes

def get_rexprt_top_k(txt_path, genes_db, k):
    print(f"\nProcessing RExPRT (Top {k})...")
    df = pd.read_csv(txt_path, sep='\t')
    
    # ensembleMax 代表集成的最大概率/分数，通常是最佳排序指标
    if 'ensembleMax' in df.columns:
        sort_col = 'ensembleMax'
    elif 'ensembleScore' in df.columns:
        sort_col = 'ensembleScore'
    else:
        raise KeyError("Could not find 'ensembleMax' or 'ensembleScore' in RExPRT file")
    
    df_sorted = df.sort_values(by=sort_col, ascending=False)
    
    # 取前 K
    df_top = df_sorted.head(k)
    print(f"  -> RExPRT Score Range ({sort_col}): {df_top[sort_col].max():.4f} to {df_top[sort_col].min():.4f}")
    
    # 提取坐标
    loci = []
    for _, row in df_top.iterrows():
        # 确保 chr 格式统一
        c = str(row['chr'])
        if not c.startswith('chr'): c = 'chr' + c
        loci.append((c, int(row['start']), int(row['end'])))
        
    # 映射基因，为了和 TREPP 公平对比，重新映射到 Ensembl ID
    # genes = map_loci_to_genes(loci, genes_db)
    genes = map_loci_to_best_gene(loci, genes_db)
    print(f"  -> Mapped to {len(genes)} genes.")
    return genes

# ================= 3. 数据库加载 =================

def get_omim_genes(path):
    print("Loading OMIM...")
    genes = set()
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            p = line.strip().split('\t')
            if len(p) < 5: continue
            if 'gene' in p[1] and p[4].strip():
                genes.add(p[4].strip()) # Ensembl ID
    return genes

def get_genetrek_genes(path):
    print("Loading GeneTrek...")
    df = pd.read_csv(path, sep='\t')
    if 'Ensembl ID' in df.columns:
        return set(df['Ensembl ID'].dropna().unique())
    return set()

# ================= 4. 绘图与统计 =================

def analyze_and_plot(trepp_set, rexprt_set, db_set, db_name, filename):
    plt.figure(figsize=(8, 8))
    
    # Venn Diagram
    venn = venn3([trepp_set, rexprt_set, db_set], set_labels=('TREPP', 'RExPRT', db_name))
    
    # 设置标签字体大小
    for label in venn.set_labels:
        if label is not None:
            label.set_fontsize(16)
    
    # 设置数字字体大小
    for text in venn.subset_labels:
        if text is not None:
            text.set_fontsize(14)

    # Colors
    colors = {'100': '#88C9BF', '010': '#7FADCA', '001': '#E67C71', 
              '110': '#F9B29C', '101': '#C6DCB9', '011': '#A2C4F1', '111': '#B8BBD9'}


    for pid, color in colors.items():
        p = venn.get_patch_by_id(pid)
        if p: p.set_color(color); p.set_alpha(0.8)
            
    plt.title(f"Top-{TOP_K} Loci Overlap with {db_name}", fontsize=16)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    # --- Statistics ---
    trepp_hits = len(trepp_set & db_set)
    rexprt_hits = len(rexprt_set & db_set)
    
    trepp_rate = trepp_hits / len(trepp_set) * 100 if len(trepp_set) > 0 else 0
    rexprt_rate = rexprt_hits / len(rexprt_set) * 100 if len(rexprt_set) > 0 else 0
    
    print(f"\n>>> Results for {db_name} (Top {TOP_K}) <<<")
    print(f"TREPP Genes mapped: {len(trepp_set)}")
    print(f"RExPRT Genes mapped: {len(rexprt_set)}")
    print(f"------------------------------------------------")
    print(f"TREPP Hit Count:  {trepp_hits} ({trepp_rate:.2f}%)")
    print(f"RExPRT Hit Count: {rexprt_hits} ({rexprt_rate:.2f}%)")
    
    if trepp_rate > rexprt_rate:
        print(f"SUCCESS: TREPP has higher precision (+{trepp_rate - rexprt_rate:.2f}%)")
    else:
        print(f"Note: RExPRT has higher precision.")

# ================= 主程序 =================

def main():
    # 1. 准备映射库
    genes_db = parse_gff3_to_genes(GFF3_FILE)
    
    # 2. 获取各自的 Top K 基因集合
    trepp_genes = get_trepp_top_k(TREPP_CSV, genes_db, TOP_K)
    rexprt_genes = get_rexprt_top_k(REXPRT_SCORE_FILE, genes_db, TOP_K)
    
    # 3. 加载金标准
    omim = get_omim_genes(OMIM_FILE)
    trek = get_genetrek_genes(GENETREK_FILE)
    
    # 4. 分析
    analyze_and_plot(trepp_genes, rexprt_genes, omim, "OMIM", OUTPUT_OMIM_IMG)
    analyze_and_plot(trepp_genes, rexprt_genes, trek, "GeneTrek", OUTPUT_GENETREK_IMG)

if __name__ == "__main__":
    main()