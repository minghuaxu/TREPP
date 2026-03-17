import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
from matplotlib.colors import LinearSegmentedColormap
from utils_load import load_processed_data
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import FuncFormatter  # 用于自定义刻度格式

# --- 统一配色 ---
custom_cmap = LinearSegmentedColormap.from_list("CustomRed", ["#F7F7F7", "#E67C71", "#C0392B"])

HG19_LENS = {
    'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276,
    'chr5': 180915260, 'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022,
    'chr9': 141213431, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
    'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753,
    'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520,
    'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566
}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prob_file', required=True)
    parser.add_argument('--feat_file', required=True)
    parser.add_argument('--output', default='fig_chromosome_norm.png')
    parser.add_argument('--bin_size', type=int, default=2000000) 
    args = parser.parse_args()
    
    df = load_processed_data(args.prob_file, args.feat_file)
    heatmap_data = []
    chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    
    for chrom in chroms:
        if chrom not in HG19_LENS: continue
        length = HG19_LENS[chrom]
        n_bins = int(np.ceil(length / args.bin_size))
        c_data = df[df['chr'] == chrom]
        
        for i in range(n_bins):
            start = i * args.bin_size
            end = (i + 1) * args.bin_size
            in_bin = c_data[(c_data['start'] >= start) & (c_data['start'] < end)]
            n_total = len(in_bin)
            n_pathogenic = len(in_bin[in_bin['prob'] > 0.5])
            
            ratio = n_pathogenic / n_total if n_total > 5 else np.nan 
            heatmap_data.append({'chr': chrom, 'bin_idx': i, 'ratio': ratio, 'start': start})
            
    hm_df = pd.DataFrame(heatmap_data)
    
    fig, ax = plt.subplots(figsize=(12, 12))
    chr_map = {c: i for i, c in enumerate(chroms)}
    
    for chrom in chroms:
        y = chr_map[chrom]
        length = HG19_LENS[chrom]
        rect = mpatches.Rectangle((0, y-0.25), length, 0.5, color='#F0F0F0', zorder=0)
        ax.add_patch(rect)
        
    sc = ax.scatter(
        hm_df['start'], 
        hm_df['chr'].map(chr_map), 
        c=hm_df['ratio'], 
        cmap=custom_cmap, 
        s=45, 
        marker='s', 
        edgecolor='none',
        vmin=0, vmax=0.3 
    )
    
    ax.set_yticks(range(len(chroms)))
    ax.set_yticklabels(chroms, fontsize=18)
    ax.invert_yaxis()
    
    # 将 X 轴单位转换为 Mb (百万碱基)
    def mb_formatter(x, pos):
        return f'{int(x / 1e6)}' # 将数值除以 1,000,000 并取整

    ax.xaxis.set_major_formatter(FuncFormatter(mb_formatter))
    
    # 更新标签，注明单位是 Mb
    ax.set_xlabel('Genomic Position (Mb)', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=18)
    
    # 使用用科学计数法且只想改字体大小，注释掉上面三行，使用下面这行：
    # ax.xaxis.get_offset_text().set_fontsize(18) 

    # --- 颜色条设置 ---
    cax = inset_axes(ax, 
                width="3%",    
                height="35%",  
                loc='lower right',  
                borderpad=10)  

    cbar = plt.colorbar(sc, cax=cax)
    cbar.set_label('Pathogenic Ratio', fontsize=16, labelpad=10)
    cbar.ax.tick_params(labelsize=14)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Saved to {args.output}")

if __name__ == "__main__":
    main()