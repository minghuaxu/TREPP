#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TREPP SHAP Analysis Tool
========================
原理：计算所有基模型(Base Models)的 SHAP 值并取平均，以解释集成模型的决策逻辑。
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import joblib
import shap
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib.colors import LinearSegmentedColormap

# 必须导入，否则 pickle 加载会报错
from catboost import CatBoostClassifier

# 设置全局绘图参数 (作为基础配置)
# 注意：这通常需要在 plt.figure() 之前设置
PLOT_CONFIG = {
    'font.size': 18,
    'axes.titlesize': 20,
    'axes.labelsize': 18,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16,
    'figure.titlesize': 22,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans', 'sans-serif'] # 优先使用 Arial
}
plt.rcParams.update(PLOT_CONFIG)


# ---------------------------------------------------------
# SHAP 分析主逻辑
# ---------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="TREPP Ensemble SHAP Analysis")
    parser.add_argument('--input', required=True, help="Input CSV (Evaluation dataset)")
    parser.add_argument('--model', required=True, help="Path to trepp_final.pkl")
    parser.add_argument('--out_dir', required=True, help="Directory to save plots")
    parser.add_argument('--samples', type=int, default=1000, help="Number of samples to explain")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # 1. 加载模型
    print(f"Loading model: {args.model}")
    artifact = joblib.load(args.model)
    base_models = artifact["base_models"]
    feat_cols = artifact["feat_cols"]
    
    print(f"Loaded {len(base_models)} base models.")
    print(f"Feature set size: {len(feat_cols)}")

    # 2. 加载数据并预处理
    print(f"Loading data: {args.input}")
    df = pd.read_csv(args.input)
    
    # 3. 采样数据
    if len(df) > args.samples:
        print(f"Sampling {args.samples} random samples for SHAP analysis...")
        df_sample = df.sample(n=args.samples, random_state=42)
    else:
        df_sample = df
        
    X = df_sample[feat_cols]
    
    # 4. 计算 Ensemble SHAP Values
    print("Calculating SHAP values (Averaging across base models)...")
    ensemble_shap_values = np.zeros(X.shape)
    
    for i, model in enumerate(tqdm(base_models, desc="Explaining Models")):
        # 处理 CalibratedClassifierCV 包装的情况
        if hasattr(model, 'calibrated_classifiers_'):
            # 使用第一个 fold 的 estimator 来解释 (通常足够代表性)
            model = model.calibrated_classifiers_[0].estimator

        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X)
        
        # 累加
        ensemble_shap_values += shap_values
        
    # 取平均
    ensemble_shap_values /= len(base_models)
    
    print("SHAP calculation done.")

    # 5. 绘图 (增强版)
    
    custom_colors = ["#7FADCA", "#E67C71"] # 蓝红配色
    custom_cmap = LinearSegmentedColormap.from_list("custom", custom_colors)
    
    # --- A. Summary Plot (Beeswarm) ---
    print("Generating Summary Plot...")
    
    # 创建画布，稍微加宽一点以容纳长特征名
    fig = plt.figure(figsize=(9, 10)) 
    
    # 绘制 SHAP
    shap.summary_plot(
        ensemble_shap_values, 
        X, 
        cmap=custom_cmap, 
        max_display=15, 
        show=False, 
        plot_size=None
    )
    
    # [关键修复] 强制修改当前 Axes 的属性
    ax = plt.gca()
    
    # 1. 修改特征名 (Y轴) 字体大小
    ax.tick_params(axis='y', labelsize=20) 
    # 2. 修改 X 轴 (SHAP Value) 字体大小
    ax.tick_params(axis='x', labelsize=20)
    # 3. 修改 X 轴标签
    ax.set_xlabel("SHAP value", fontsize=20, labelpad=10)
    
    # 4. 修改 Colorbar 字体 (SHAP 的 Colorbar 是独立的 Axes)
    # 通常它是 figure 中的最后一个 axes
    if len(fig.axes) > 1:
        cbar_ax = fig.axes[-1]
        cbar_ax.tick_params(labelsize=20)
        cbar_ax.set_ylabel("Feature value", fontsize=20, labelpad=10)

    # plt.title("TREPP Feature Importance (Ensemble Average)", fontsize=22, pad=20)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, "shap_summary_beeswarm.png"), dpi=600, bbox_inches='tight')
    plt.close()
    

    # --- B. Bar Plot (绝对重要性) ---
    print("Generating Bar Plot...")
    fig = plt.figure(figsize=(9, 10))
    
    shap.summary_plot(
        ensemble_shap_values, 
        X, 
        plot_type="bar", 
        show=False, 
        plot_size=None, 
        max_display=15, 
        color="#7FADCA"
    )
    
    # [关键修复] 强制修改属性
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=20) # 特征名
    ax.tick_params(axis='x', labelsize=20) # SHAP值
    ax.set_xlabel("mean(|SHAP value|)", fontsize=20, labelpad=10)
    
    # 可以在条形图上增加数值标签
    # for p in ax.patches:
    #     ax.annotate(f"{p.get_width():.2f}", (p.get_width() * 1.01, p.get_y() + p.get_height()/2), fontsize=12, va='center')

    # plt.title("TREPP Mean Absolute SHAP Values", fontsize=22, pad=20)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, "shap_importance_bar.png"), dpi=600, bbox_inches='tight')
    plt.close()

    print(f"Plots saved to {args.out_dir}")

if __name__ == "__main__":
    main()

# python scripts/analyze_shap.py \
#     --input data/processed/eval_features.csv \
#     --model data/models_2/trepp_final.pkl \
#     --out_dir results/shap_plots \
#     --samples 2000