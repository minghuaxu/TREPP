import pandas as pd
import numpy as np
import sqlite3
import pysam
import argparse
import os
import subprocess
import tempfile
import shutil
from tqdm import tqdm

def add_engineered_features(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    # -----------------------------
    # 1. 基本重复结构
    # -----------------------------
    motif_len_safe = df["motif_len"].replace(0, np.nan)
    df["repeat_units"] = df["tr_len"] / motif_len_safe
    df["repeat_units"] = df["repeat_units"].fillna(0.0)
    df["log_tr_len"] = np.log10(df["tr_len"].values + 1.0)
    df["log_repeat_units"] = np.log10(df["repeat_units"].values + 1e-3)

    # -----------------------------
    # 2. motif 碱基组成特征
    # -----------------------------
    def _motif_stats(m):
        m = str(m).upper()
        # 处理空 motif 的情况
        if not m or m == 'NAN': 
            return pd.Series({
                "motif_gc": 0.0, "motif_at": 0.0, 
                "motif_purine_frac": 0.0, "motif_pyrimidine_frac": 0.0, 
                "motif_has_cpg": 0
            })
            
        bases = [b for b in m if b in "ACGT"]
        L = len(bases) if len(bases) > 0 else (len(m) if len(m) > 0 else 1)
        A = bases.count("A")
        C = bases.count("C")
        G = bases.count("G")
        T = bases.count("T")
        gc = (G + C) / L
        at = (A + T) / L
        pur = (A + G) / L
        pyr = (C + T) / L
        has_cpg = 1 if "CG" in m else 0
        return pd.Series(
            {
                "motif_gc": gc,
                "motif_at": at,
                "motif_purine_frac": pur,
                "motif_pyrimidine_frac": pyr,
                "motif_has_cpg": has_cpg,
            }
        )

    motif_feats = df["motif"].apply(_motif_stats)
    df = pd.concat([df, motif_feats], axis=1)

    # -----------------------------
    # 3. motif 类型/长度类别
    # -----------------------------
    df["is_trinuc"] = (df["motif_len"] == 3).astype(int)
    df["is_pentanuc"] = (df["motif_len"] == 5).astype(int)

    df["log_tr_gc"] = np.log10(df["tr_gc"].values + 1e-3)
    df["log_motif_gc"] = np.log10(df["motif_gc"].values + 1e-3)
    df["delta_gc_tr_motif"] = df["log_tr_gc"] - df["log_motif_gc"]

    # -----------------------------
    # 4. STR 密度对比
    # -----------------------------
    eps = 1e-6
    if "str_density_1000" in df.columns and "str_density_10000" in df.columns:
        df["str_density_ratio"] = df["str_density_1000"] / (df["str_density_10000"] + eps)
    else:
        df["str_density_ratio"] = 0.0
        
    df["log_str_density_ratio"] = np.log10(df["str_density_ratio"].values + 1e-3)

    # -----------------------------
    # 5. GC 相对差异 & STR 密度对比
    # -----------------------------
    df["delta_gc_tr_flank"] = df["tr_gc"] - df["flank_gc_1k"]
    df["abs_delta_gc_tr_flank"] = df["delta_gc_tr_flank"].abs()

    # -----------------------------
    # 6. 基因/位置相关特征
    # -----------------------------
    gene_len_safe = df["gene_len"].replace(0, 1)
    df["norm_dist_tss"] = df["dist_to_tss"] / gene_len_safe

    df["in_gene"] = (df["gene_name"] != ".").astype(int)
    df["in_exon"] = (df["dist_to_exon"] == 0).astype(int)
    df["in_cds"] = (df["dist_to_cds"] == 0).astype(int)
    df["near_tss_1k"] = (df["dist_to_tss"] <= 1000).astype(int)
    df["near_tss_5k"] = (df["dist_to_tss"] <= 5000).astype(int)

    return df

# ==========================================
# 2. 原始特征提取逻辑 (DB & TRF)
# ==========================================

class DBQuery:
    def __init__(self, db_path):
        self.conn = sqlite3.connect(db_path)
        self.cursor = self.conn.cursor()
    
    def close(self):
        self.conn.close()

def get_min_dist(chrom, start, end, table, db, max_dist=50000):
    s_search = max(0, start - max_dist)
    e_search = end + max_dist
    
    def _query(c):
        db.cursor.execute(f"SELECT start, end FROM {table} WHERE chrom=? AND end >= ? AND start <= ?", (c, s_search, e_search))
        return db.cursor.fetchall()

    rows = _query(chrom)
    if not rows:
        alt = chrom.replace('chr', '') if str(chrom).startswith('chr') else 'chr'+str(chrom)
        rows = _query(alt)

    if not rows: return max_dist

    min_d = max_dist
    for fs, fe in rows:
        if max(start, fs) < min(end, fe):
            return 0 
        d = min(abs(start - fe), abs(end - fs))
        if d < min_d: min_d = d
    return min_d

def get_gene_advanced(chrom, start, end, db):
    search_limit = 50000
    
    def _find(c):
        db.cursor.execute("SELECT start, end, strand, gene_name FROM anno_gene WHERE chrom=? AND end >= ? AND start <= ? LIMIT 1", (c, start, end))
        res = db.cursor.fetchone()
        if res: return res
        db.cursor.execute("SELECT start, end, strand, gene_name FROM anno_gene WHERE chrom=? AND end >= ? AND start <= ?", (c, start-search_limit, end+search_limit))
        rows = db.cursor.fetchall()
        best, min_d = None, search_limit
        for r in rows:
            d = min(abs(start-r[1]), abs(end-r[0]))
            if d < min_d: min_d = d; best = r
        return best

    gene_info = _find(chrom)
    if not gene_info:
        alt = chrom.replace('chr', '') if str(chrom).startswith('chr') else 'chr'+str(chrom)
        gene_info = _find(alt)
    
    dist_tss = 50000
    rel_pos = -1.0
    g_len = 0
    pLI, oe_lof = 0.0, 1.0
    gene_name = "."
    
    if gene_info:
        gs, ge, strand, gname = gene_info
        gene_name = gname
        g_len = ge - gs
        tr_center = (start + end) / 2
        
        if strand == '+':
            tss = gs
            dist_tss = abs(tr_center - tss)
            rel_pos = (tr_center - gs) / g_len if g_len > 0 else 0
        else:
            tss = ge
            dist_tss = abs(tr_center - tss)
            rel_pos = (ge - tr_center) / g_len if g_len > 0 else 0
            
        rel_pos = max(0.0, min(1.0, rel_pos))
        
        try:
            db.cursor.execute("SELECT pLI, oe_lof FROM gene_metrics WHERE gene=?", (gname,))
            m = db.cursor.fetchone()
            if m: pLI, oe_lof = float(m[0]), float(m[1])
        except: pass

    return dist_tss, rel_pos, g_len, pLI, oe_lof, gene_name

def compute_trf_density(df, fasta, trf_path):
    import os
    print("Calculating TRF Density...")
    radii = [1000, 10000]
    max_r = 10000
    cwd = os.getcwd()
    tmp_dir = tempfile.mkdtemp()
    
    try:
        os.chdir(tmp_dir)
        with open("batch.fa", 'w') as f:
            for i, row in df.iterrows():
                try:
                    seq = fasta.fetch(str(row['chr']), max(0, row['start']-max_r), row['end']+max_r).upper()
                    f.write(f">{i}\n{seq}\n")
                except: pass
        
        subprocess.run([trf_path, "batch.fa", "2", "5", "7", "80", "10", "30", "2000", "-l", "2", "-h", "-d"], 
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        dat_files = [x for x in os.listdir('.') if x.endswith('.dat')]
        intervals = {}
        if dat_files:
            curr = None
            with open(dat_files[0], 'r') as f:
                for line in f:
                    if line.startswith("Sequence:"): 
                        try: curr = int(line.split()[1]); intervals[curr] = []
                        except: pass
                    elif curr is not None:
                        p = line.split()
                        if len(p)>5 and p[0].isdigit(): intervals[curr].append((int(p[0]), int(p[1])))
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp_dir)

    res = {r: [0.0]*len(df) for r in radii}
    for i in df.index:
        tr_len = df.loc[i, 'end'] - df.loc[i, 'start']
        ints = intervals.get(i, [])
        merged = []
        if ints:
            ints.sort()
            cs, ce = ints[0]
            for ns, ne in ints[1:]:
                if ns < ce: ce = max(ce, ne)
                else: merged.append((cs, ce)); cs, ce = ns, ne
            merged.append((cs, ce))
            
        center = max_r
        for r in radii:
            ws, we = center-r, center+tr_len+r
            wlen = we-ws
            tot = 0
            for s, e in merged:
                os_ov, oe_ov = max(ws, s), min(we, e)
                if os_ov < oe_ov: tot += (oe_ov - os_ov)
            res[r][i] = tot/wlen if wlen>0 else 0.0
            
    return pd.DataFrame({f'str_density_{r}': res[r] for r in radii})

def calculate_gc(fasta, chrom, start, end):
    try:
        tr_seq = fasta.fetch(chrom, start, end).upper()
        tr_gc = (tr_seq.count('G') + tr_seq.count('C')) / len(tr_seq) if tr_seq else 0
        
        seq_1k = fasta.fetch(chrom, max(0, start-1000), end+1000).upper()
        flank_gc = (seq_1k.count('G') + seq_1k.count('C')) / len(seq_1k) if seq_1k else 0
        return tr_gc, flank_gc
    except: return 0.0, 0.0

# ==========================================
# 主流程
# ==========================================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_bed', required=True)
    parser.add_argument('--ref_fasta', required=True)
    parser.add_argument('--db_path', required=True)
    parser.add_argument('--trf_path', required=True)
    parser.add_argument('--output_csv', required=True)
    args = parser.parse_args()

    # 1. Load Data
    raw_df = pd.read_csv(args.input_bed, sep='\t', header=None)
    raw_df = raw_df.iloc[:, :4]
    raw_df.columns = ['chr', 'start', 'end', 'motif']

    raw_df = raw_df.reset_index(drop=True)

    # 2. TRF Density
    fasta = pysam.FastaFile(args.ref_fasta)
    dens_df = compute_trf_density(raw_df, fasta, args.trf_path)

    # 3. DB Features
    print("Extracting Raw Features...")
    db = DBQuery(args.db_path)
    feat_rows = []
    
    for idx, row in tqdm(raw_df.iterrows(), total=len(raw_df)):
        c, s, e = str(row['chr']), int(row['start']), int(row['end'])
        
        # Raw Calcs
        d_exon = get_min_dist(c, s, e, 'anno_exon', db)
        d_cds = get_min_dist(c, s, e, 'anno_cds', db)
        d_utr = get_min_dist(c, s, e, 'anno_utr', db)
        dist_tss, rel_pos, g_len, pLI, oe_lof, gname = get_gene_advanced(c, s, e, db)
        tr_gc, flank_gc = calculate_gc(fasta, c, s, e)
        
        item = {
            'chr': c, 'start': s, 'end': e, 
            'motif': str(row.get('motif', '')),
            'gene_name': gname,
            
            # Raw Numerics
            'motif_len': len(str(row.get('motif', ''))),
            'tr_len': e - s,
            'tr_gc': tr_gc,
            'flank_gc_1k': flank_gc,
            'dist_to_exon': d_exon,
            'dist_to_cds': d_cds,
            'dist_to_utr': d_utr,
            'dist_to_tss': dist_tss,
            'gene_len': g_len,
            'gene_rel_pos': rel_pos,
            'pLI': pLI,
            'oe_lof': oe_lof
        }
        feat_rows.append(item)
    
    db.close()
    
    # 4. Merge & Clip
    df_raw = pd.DataFrame(feat_rows)
    df_combined = pd.concat([df_raw, dens_df], axis=1)
    
    for c in ['dist_to_exon', 'dist_to_cds', 'dist_to_utr', 'dist_to_tss']:
        df_combined[c] = df_combined[c].clip(upper=50000)

    # 5. Apply Feature Engineering HERE
    df_final = add_engineered_features(df_combined)
    
    # 6. Save
    df_final.to_csv(args.output_csv, index=False)
    print(f"Saved to {args.output_csv}")

if __name__ == "__main__":
    main()