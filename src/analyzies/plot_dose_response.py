import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from scipy.stats import spearmanr
from utils_load import load_processed_data

# --- 统一配色 ---
COLOR_TREND = '#E67C71' # Red for pathogenic trend

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prob_file', required=True)
    parser.add_argument('--feat_file', required=True)
    parser.add_argument('--output', default='fig_dose_response.png')
    args = parser.parse_args()
    
    df = load_processed_data(args.prob_file, args.feat_file)
    df['is_functional'] = df[['is_promoter', 'is_exon', 'is_utr']].max(axis=1)
    df['score_bin'] = pd.cut(df['prob'], bins=np.linspace(0, 1, 11), include_lowest=True)
    
    agg = df.groupby('score_bin').agg(
        prop_functional=('is_functional', 'mean'),
        count=('is_functional', 'size'),
        mean_score=('prob', 'mean')
    ).reset_index()
    agg['prop_pct'] = agg['prop_functional'] * 100
    
    # 去除NaN
    agg_clean = agg.dropna(subset=['mean_score', 'prop_pct'])
    
    if len(agg_clean) < 2:
        corr, pval = 0, 1
    else:
        corr, pval = spearmanr(agg_clean['mean_score'], agg_clean['prop_pct'])
    
    # Plot
    plt.figure(figsize=(8, 5))
    
    # 使用统一的红色
    sns.pointplot(data=agg_clean, x='score_bin', y='prop_pct', color=COLOR_TREND, scale=1.2)
    
    plt.xticks(rotation=45, fontsize=14)
    plt.xlabel("TREPP Pathogenicity Score (Binned)", fontsize=12)
    plt.ylabel("% of TRs in Functional Regions \n (Promoter, UTR, Exon)", fontsize=12)
    # plt.title(f"Dose-Response: Score vs. Genomic Context\n(Spearman r={corr:.2f}, P={pval:.2e})", fontsize=14)
    
    plt.grid(True, axis='y', linestyle='--', alpha=0.5, color='#CCCCCC')
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Saved to {args.output}")

if __name__ == "__main__":
    main()