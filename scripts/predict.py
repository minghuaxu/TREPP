#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TREPP Prediction Script (Compatible with V7 Training)
=====================================================
功能：
1. 加载 V7 训练生成的 trepp_final.pkl
2. 对新数据应用完全相同的特征工程 (Feature Engineering)
3. 对齐特征列 (Feature Alignment)
4. 执行 Stacked Inference (Base Models -> Meta Model)
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import joblib
from tqdm import tqdm

# 即使只加载模型，也需要导入 CatBoost，否则 pickle 反序列化会找不到类定义
from catboost import CatBoostClassifier
from sklearn.calibration import CalibratedClassifierCV

# ---------------------------------------------------------
# 1. 特征工程 (必须与 train_final_v7.py 完全一致)
# ---------------------------------------------------------
def add_engineered_features(df: pd.DataFrame) -> pd.DataFrame:
    """
    复现训练时的特征工程逻辑
    """
    df = df.copy()

    # 1. 基本重复结构
    motif_len_safe = df["motif_len"].replace(0, np.nan)
    df["repeat_units"] = df["tr_len"] / motif_len_safe
    df["repeat_units"] = df["repeat_units"].fillna(0.0)
    df["log_tr_len"] = np.log10(df["tr_len"].values + 1.0)
    df["log_repeat_units"] = np.log10(df["repeat_units"].values + 1e-3)

    # 2. motif 碱基组成特征
    def _motif_stats(m):
        m = str(m).upper()
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

    # 避免空 dataframe 报错
    if len(df) > 0:
        motif_feats = df["motif"].apply(_motif_stats)
        df = pd.concat([df, motif_feats], axis=1)
    else:
        # 如果是空的，手动补列防止后续报错
        for c in ["motif_gc", "motif_at", "motif_purine_frac", "motif_pyrimidine_frac", "motif_has_cpg"]:
            df[c] = 0.0

    # 3. motif 类型/长度类别
    df["is_trinuc"] = (df["motif_len"] == 3).astype(int)
    df["is_pentanuc"] = (df["motif_len"] == 5).astype(int)

    df["log_tr_gc"] = np.log10(df["tr_gc"].values + 1e-3)
    # 处理可能存在的缺失列 (如果输入数据不完整)
    if "motif_gc" not in df.columns: df["motif_gc"] = 0
    df["log_motif_gc"] = np.log10(df["motif_gc"].values + 1e-3)
    df["delta_gc_tr_motif"] = df["log_tr_gc"] - df["log_motif_gc"]

    # 4. STR 密度对比
    eps = 1e-6
    # 确保列存在
    if "str_density_1000" in df.columns and "str_density_10000" in df.columns:
        df["str_density_ratio"] = df["str_density_1000"] / (df["str_density_10000"] + eps)
    else:
        df["str_density_ratio"] = 0.0
    
    df["log_str_density_ratio"] = np.log10(df["str_density_ratio"].values + 1e-3)

    # 5. GC 相对差异
    if "flank_gc_1k" in df.columns:
        df["delta_gc_tr_flank"] = df["tr_gc"] - df["flank_gc_1k"]
    else:
        df["delta_gc_tr_flank"] = 0.0
    df["abs_delta_gc_tr_flank"] = df["delta_gc_tr_flank"].abs()

    # 6. 基因/位置相关
    if "gene_len" in df.columns:
        gene_len_safe = df["gene_len"].replace(0, 1)
        df["norm_dist_tss"] = df["dist_to_tss"] / gene_len_safe
    else:
        df["norm_dist_tss"] = 0.0

    df["in_gene"] = (df["gene_name"] != ".").astype(int)
    df["in_exon"] = (df["dist_to_exon"] == 0).astype(int)
    df["in_cds"] = (df["dist_to_cds"] == 0).astype(int)
    df["near_tss_1k"] = (df["dist_to_tss"] <= 1000).astype(int)
    df["near_tss_5k"] = (df["dist_to_tss"] <= 5000).astype(int)

    return df

# ---------------------------------------------------------
# 2. 模型加载与预测逻辑
# ---------------------------------------------------------
def load_model(model_path):
    print(f"Loading model artifact from: {model_path}")
    try:
        artifact = joblib.load(model_path)
        required_keys = ["base_models", "meta_info", "feat_cols"]
        for key in required_keys:
            if key not in artifact:
                raise ValueError(f"Model file is missing key: {key}. Are you using a V7 model?")
        return artifact
    except Exception as e:
        print(f"Error loading model: {e}")
        sys.exit(1)

def predict(df, artifact):
    """
    执行预测流程：
    1. 特征工程
    2. 特征筛选与对齐
    3. 基模型预测
    4. 元模型预测
    """
    base_models = artifact["base_models"]
    meta_info = artifact["meta_info"]
    
    # 训练时最终保留的特征列名 (经过了 SelectFeatures 筛选的)
    train_feat_cols = artifact["feat_cols"] 
    
    # 1. 应用特征工程
    # print("Applying feature engineering...")
    # df_eng = add_engineered_features(df)
    
    # 2. 特征对齐 (Feature Alignment)
    print(f"Aligning features to match training set ({len(train_feat_cols)} features)...")
    
    # 检查缺失列并补0
    # missing_cols = set(train_feat_cols) - set(df_eng.columns)
    # if missing_cols:
    #     print(f"Warning: {len(missing_cols)} features missing, filling with 0.")
    #     for c in missing_cols:
    #         df_eng[c] = 0.0
            
    # 严格按照训练时的顺序提取特征，并转为 float
    # 这一步至关重要，CatBoost 对列顺序敏感
    try:
        X = df[train_feat_cols].values.astype(float)
    except KeyError as e:
        print(f"Error during feature alignment: {e}")
        sys.exit(1)
        
    # 3. 基模型预测 (Base Models Inference)
    n_samples = len(X)
    n_models = len(base_models)
    base_preds = np.zeros((n_samples, n_models))
    
    print(f"Running inference with {n_models} base models...")
    for i, model in enumerate(tqdm(base_models, desc="Base Models")):
        # model 是 CalibratedClassifierCV
        base_preds[:, i] = model.predict_proba(X)[:, 1]
        
    # 4. 元模型预测 (Meta Model Inference)
    print("Running Meta-Model (Stacking)...")
    
    # 这里的 meta_model 是训练好的 LogisticRegression (Calibrated)
    # 它需要的输入是基模型的预测概率
    # 注意：V7 代码中，meta_model 是用 oof[:, selected_idx] 训练的
    
    selected_idx = meta_info['selected_idx']
    meta_model = meta_info['meta_model']
    
    # 筛选列
    X_meta = base_preds[:, selected_idx]
    
    # 预测最终概率
    final_probs = meta_model.predict_proba(X_meta)[:, 1]
    
    return final_probs

def main():
    parser = argparse.ArgumentParser(description="TREPP V7 Prediction Tool")
    parser.add_argument('--input', required=True, help="Input CSV file (raw features)")
    parser.add_argument('--model', required=True, help="Path to trepp_final.pkl")
    parser.add_argument('--output', required=True, help="Output CSV path")
    parser.add_argument('--threshold', type=float, default=None, help="Custom threshold (optional)")
    
    args = parser.parse_args()
    
    # 1. 加载模型
    artifact = load_model(args.model)
    
    # 确定阈值
    if args.threshold is not None:
        threshold = args.threshold
        print(f"Using user-defined threshold: {threshold}")
    else:
        # 优先使用 Macro F1 最佳阈值
        threshold = 0.5

    # 2. 加载数据
    print(f"Reading input data: {args.input}")
    try:
        df = pd.read_csv(args.input)
    except Exception as e:
        print(f"Error reading input CSV: {e}")
        sys.exit(1)
        
    # 3. 构造 ID 列 (如果不存在)
    if 'id' not in df.columns:
        # 尝试根据 chr, start, end, motif 构造
        required_id_cols = ['chr', 'start', 'end', 'motif']
        if all(c in df.columns for c in required_id_cols):
            df['id'] = df['chr'].astype(str) + '_' + df['start'].astype(str) + '_' + df['end'].astype(str) + '_' + df['motif'].astype(str)
        else:
            print("Warning: Could not construct 'id' column. Using row index.")
            df['id'] = df.index.astype(str)

    # 4. 执行预测
    probs = predict(df, artifact)
    
    # 5. 构造输出
    output_df = pd.DataFrame()
    output_df['id'] = df['id']
    
    # 如果原始数据有 label，保留它以便比对
    if 'label' in df.columns:
        output_df['true_label'] = df['label']
        
    output_df['prob'] = probs
    output_df['pred_label'] = (probs >= threshold).astype(int)
    
    output_df = output_df[['id','prob']]
    # 6. 保存
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    output_df.to_csv(args.output, index=False, sep='\t')
    
    print(f"\nPredictions saved to: {args.output}")
    print(f"Total samples: {len(output_df)}")
    # print(f"Positive predictions: {output_df['pred_label'].sum()}")

if __name__ == "__main__":
    main()