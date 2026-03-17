import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
import os

# 设置绘图风格
sns.set(style="whitegrid", context="paper", font_scale=1.2)

def load_data(csv_path):
    print(f"Loading data from {csv_path}...")
    df = pd.read_csv(csv_path)
    
    # 排除非数值列和元数据
    meta_cols = ['chr', 'start', 'end', 'motif', 'label', 'gene_name', 'id']
    feature_cols = [c for c in df.columns if c not in meta_cols]
    
    X = df[feature_cols]
    y = df['label']
    
    # 填充缺失值
    imputer = SimpleImputer(strategy='median')
    X_filled = pd.DataFrame(imputer.fit_transform(X), columns=feature_cols)
    
    return X_filled, y, feature_cols

def plot_dim_reduction(X, y, out_dir):
    print("Running Dimensionality Reduction (PCA & t-SNE)...")
    
    # 标准化
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # 1. PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    var_explained = pca.explained_variance_ratio_
    
    # 2. t-SNE 
    perp = min(30, len(X) - 1)
    tsne = TSNE(n_components=2, perplexity=perp, random_state=42, init='pca', learning_rate='auto')
    X_tsne = tsne.fit_transform(X_scaled)
    
    # 绘图
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # 定义颜色: 0->Blue(Benign), 1->Red(Pathogenic)
    palette = {0: "#3498db", 1: "#e74c3c"}
    
    # Plot PCA
    sns.scatterplot(
        x=X_pca[:, 0], y=X_pca[:, 1], hue=y, palette=palette, 
        alpha=0.7, s=60, edgecolor='k', ax=axes[0]
    )
    axes[0].set_title(f'PCA (Explained Var: {sum(var_explained):.2%})')
    axes[0].set_xlabel('PC1')
    axes[0].set_ylabel('PC2')
    
    # Plot t-SNE
    sns.scatterplot(
        x=X_tsne[:, 0], y=X_tsne[:, 1], hue=y, palette=palette, 
        alpha=0.7, s=60, edgecolor='k', ax=axes[1]
    )
    axes[1].set_title('t-SNE Visualization')
    
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "plot_clustering.png"), dpi=300)
    print("Saved plot_clustering.png")

def plot_feature_importance(X, y, feature_names, out_dir):
    print("Calculating Feature Importance...")
    rf = RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced')
    rf.fit(X, y)
    
    importances = rf.feature_importances_
    indices = np.argsort(importances)[::-1]
    
    # 取前 15 个重要特征
    top_n = min(15, len(feature_names))
    top_indices = indices[:top_n]
    
    plt.figure(figsize=(10, 8))
    sns.barplot(x=importances[top_indices], y=np.array(feature_names)[top_indices], palette="viridis")
    plt.title(f'Top {top_n} Feature Importance (Random Forest)')
    plt.xlabel('Importance')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "plot_importance.png"), dpi=300)
    print("Saved plot_importance.png")
    
    return np.array(feature_names)[top_indices]

def plot_top_distributions(X, y, top_features, out_dir):
    print("Plotting distributions for top features...")
    # 取前 6 个特征画小提琴图
    top_6 = top_features[:6]
    
    data = X[top_6].copy()
    data['Label'] = y.map({0: 'Benign', 1: 'Pathogenic'})
    
    # 融化数据以便绘图
    data_melt = data.melt(id_vars='Label', var_name='Feature', value_name='Value')
    
    # 由于不同特征尺度差异巨大（距离 vs GC含量），分图画
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    
    for i, feature in enumerate(top_6):
        sns.violinplot(
            data=data, x='Label', y=feature, 
            palette={"Benign": "#3498db", "Pathogenic": "#e74c3c"},
            ax=axes[i], split=False, inner="quartile"
        )
        axes[i].set_title(feature)
        axes[i].set_xlabel('')
        
        # 如果是距离特征，需要对数坐标看才清楚
        if 'dist' in feature or 'len' in feature:
            # 检查是否有0，如果有0不能直接log，加1
            if (data[feature] > 0).all():
                axes[i].set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "plot_distributions.png"), dpi=300)
    print("Saved plot_distributions.png")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help="Path to train_features.csv")
    parser.add_argument('--out_dir', required=True, help="Output directory for plots")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    
    X, y, feat_names = load_data(args.input)
    
    print(f"Data shape: {X.shape}, Positive samples: {y.sum()}")
    
    # 1. 降维可视化
    plot_dim_reduction(X, y, args.out_dir)
    
    # 2. 特征重要性
    top_feats = plot_feature_importance(X, y, feat_names, args.out_dir)
    
    # 3. 核心特征分布对比
    plot_top_distributions(X, y, top_feats, args.out_dir)

if __name__ == "__main__":
    main()