import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import argparse
from utils_load import load_processed_data

# --- 统一配色方案 ---
COLOR_RED = '#E67C71'   # Enriched (OR > 1)
COLOR_BLUE = '#7FADCA'  # Depleted (OR < 1)
COLOR_GREEN = '#88C9BF' # Reference (OR = 1)

def perform_matching(df, threshold=0.5):
    print("Performing stratified matching...")
    high_risk = df[df['prob'] > threshold].copy()
    pool = df[df['prob'] <= threshold].copy()
    
    try:
        ret = pd.qcut(df['tr_len'], q=5, retbins=True, duplicates='drop')
        len_bins = ret[1]
    except Exception as e:
        print(f"Warning: qcut failed ({e}), using simple linspace bins.")
        len_bins = np.linspace(df['tr_len'].min(), df['tr_len'].max(), 6)

    gc_bins = np.linspace(0, 1, 6)
    
    for d in [high_risk, pool]:
        d['bin_id'] = (
            d['motif_len'].astype(str) + "_" + 
            pd.cut(d['tr_gc'], bins=gc_bins, labels=False, include_lowest=True).astype(str) + "_" +
            pd.cut(d['tr_len'], bins=len_bins, labels=False, include_lowest=True).astype(str)
        )

    matched_controls = []
    pool_groups = pool.groupby('bin_id')
    
    matched_count = 0
    total_needed = 0
    
    for bin_id, group in high_risk.groupby('bin_id'):
        n_needed = len(group)
        total_needed += n_needed
        
        if bin_id in pool_groups.groups:
            candidates = pool_groups.get_group(bin_id)
            n_sample = min(len(candidates), n_needed * 1) 
            sampled = candidates.sample(n=n_sample, replace=(len(candidates) < n_needed), random_state=42)
            matched_controls.append(sampled)
            matched_count += len(sampled)
            
    if not matched_controls:
        raise ValueError("Matching failed! No control samples found.")
        
    matched_df = pd.concat(matched_controls)
    print(f"Matching stats: High Risk={len(high_risk)}, Matched Control={len(matched_df)} (Coverage: {matched_count/total_needed:.1%})")
    
    return high_risk, matched_df

def calculate_enrichment(case_df, ctrl_df, feature_col):
    a = case_df[feature_col].sum()
    b = len(case_df) - a
    c = ctrl_df[feature_col].sum()
    d = len(ctrl_df) - c
    
    # Haldane-Anscombe correction
    if a == 0 or b == 0 or c == 0 or d == 0:
        a_ = a + 0.5
        b_ = b + 0.5
        c_ = c + 0.5
        d_ = d + 0.5
    else:
        a_, b_, c_, d_ = a, b, c, d
        
    odds_ratio = (a_ * d_) / (b_ * c_)
    
    # Fisher test uses original integers
    _, p_value = stats.fisher_exact([[a, b], [c, d]])
    
    log_or = np.log(odds_ratio)
    se = np.sqrt(1/a_ + 1/b_ + 1/c_ + 1/d_)
    
    lower_ci = np.exp(log_or - 1.96 * se)
    upper_ci = np.exp(log_or + 1.96 * se)
    
    return odds_ratio, lower_ci, upper_ci, p_value

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prob_file', required=True)
    parser.add_argument('--feat_file', required=True)
    parser.add_argument('--output', default='fig_enrichment_matched.png')
    args = parser.parse_args()
    
    df = load_processed_data(args.prob_file, args.feat_file)
    cases, controls = perform_matching(df, threshold=0.5)
    
    features = {
        'Promoter': 'is_promoter',
        'Exon': 'is_exon',
        'UTR': 'is_utr',
        'Intron': 'is_intron',
        'Intergenic': 'is_intergenic'
    }
    
    results = []
    for label, col in features.items():
        if col in cases.columns:
            orr, low, up, p = calculate_enrichment(cases, controls, col)
            results.append({'Label': label, 'OR': orr, 'Lower': low, 'Upper': up, 'P': p})
            
    res_df = pd.DataFrame(results).iloc[::-1]
    
    # Plot
    plt.figure(figsize=(10, 6)) # 稍微加宽一点
    
    colors = [COLOR_RED if row.OR > 1 else COLOR_BLUE for row in res_df.itertuples()]
    y_pos = np.arange(len(res_df))
    
    # 获取X轴的最大值，用于动态计算文字偏移量
    max_or = res_df['Upper'].max()
    # 动态偏移量：X轴总长度的 2%
    text_offset = max_or * 0.02
    
    for i, row in enumerate(res_df.itertuples()):
        err = [[row.OR - row.Lower], [row.Upper - row.OR]]
        
        plt.errorbar(
            x=row.OR, 
            y=i, 
            xerr=err, 
            fmt='o', 
            color=colors[i], 
            ecolor=colors[i], 
            capsize=5, 
            markersize=10,
            elinewidth=3,
            capthick=2
        )
    
    plt.axvline(x=1, color=COLOR_GREEN, linestyle='--', linewidth=2, alpha=0.8)
    plt.yticks(y_pos, res_df['Label'], fontsize=16)
    plt.xticks(fontsize=16)
    plt.xlabel("Odds Ratio (vs. Matched Background)", fontsize=16)
    # plt.title("Genomic Region Enrichment (Matched)", fontsize=16)
    
    # Add text
    for i, row in enumerate(res_df.itertuples()):
        sig = "***" if row.P < 0.001 else "**" if row.P < 0.01 else "*" if row.P < 0.05 else "ns"
        
        # 无论红蓝，文字都放在点的右侧 这样蓝色文字就不会撞到 Y 轴了
        if row.OR > 1:
            pos_x = row.Upper + text_offset
            txt_color = COLOR_RED
        else:
            # 对于蓝色 (OR < 1)，Error bar 通常很短或接近0，放在 OR 值 (接近0) 的右侧，加上偏移量，会显示在点和参考线(1.0)之间
            pos_x = max(row.Upper, row.OR) + text_offset
            txt_color = COLOR_BLUE
            
        plt.text(pos_x, i, f"OR={row.OR:.2f} {sig}", 
                 va='center', ha='left',  # 统一左对齐 (向右延伸)
                 fontsize=12, color=txt_color, fontweight='bold')

    plt.grid(axis='x', linestyle='--', alpha=0.3)
    
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Saved to {args.output}")

if __name__ == "__main__":
    main()