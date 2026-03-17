#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Strict Inference & Evaluation Script
====================================
目标: 直接加载训练好的物理模型文件，在 Eval 数据集上进行单向预测，并计算核心评估指标。
"""

import argparse
import os
import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score, 
    roc_auc_score, 
    f1_score, 
    precision_score, 
    recall_score
)

def evaluate(model_path, eval_path):
    print(f"[*] 正在加载模型文件: {model_path}")
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"找不到模型文件: {model_path}")
        
    # 1. 加载模型 Artifact
    model_artifact = joblib.load(model_path)
    
    # 提取组件（假设你的 pkl 保存的是包含特征列和两层模型的字典）
    try:
        feat_cols = model_artifact["feat_cols"]
        base_models = model_artifact["base_models"]
        meta_info = model_artifact["meta_info"]
        meta_model = meta_info["meta_model"]
        # threshold = meta_info.get("threshold", 0.5) # 如果没有保存阈值，默认0.5
        threshold = 0.5
    except KeyError as e:
        raise KeyError(f"模型文件结构不匹配，缺失关键键值: {e}。请确认 pkl 文件的保存格式。")

    print(f"[*] 正在加载 Eval 数据集: {eval_path}")
    df_eval = pd.read_csv(eval_path)
    
    missing_cols = [c for c in feat_cols if c not in df_eval.columns]
    if missing_cols:
        raise ValueError(f"Eval 数据集缺失必要的特征列: {missing_cols}")
        
    X_eval = df_eval[feat_cols].values.astype(float)
    y_eval = df_eval["label"].values.astype(int)
    
    print(f"[*] Eval 数据集形状: {X_eval.shape}, 正样本: {sum(y_eval)}, 负样本: {len(y_eval)-sum(y_eval)}")
    print("[*] 正在执行前向预测...")

    # 2. Level-1: 基模型预测
    level_1_eval = np.zeros((len(X_eval), len(base_models)))
    for i, model in enumerate(base_models):
        level_1_eval[:, i] = model.predict_proba(X_eval)[:, 1]
        
    # 3. Level-2: Meta 模型预测最终概率
    final_probs = meta_model.predict_proba(level_1_eval)[:, 1]
    
    # 4. 应用阈值生成硬标签
    final_preds = (final_probs >= threshold).astype(int)
    
    # 5. 计算指标
    print("[*] 正在计算评估指标...")
    auprc = average_precision_score(y_eval, final_probs)
    auroc = roc_auc_score(y_eval, final_probs)
    f1_mac = f1_score(y_eval, final_preds, average="macro")
    prec_mac = precision_score(y_eval, final_preds, average="macro", zero_division=0)
    rec_mac = recall_score(y_eval, final_preds, average="macro", zero_division=0)
    
    # 输出结果
    print("\n" + "="*45)
    print(" 最终评估结果 (Eval Dataset)")
    print("="*45)
    print(f" 使用的分类阈值 (Threshold): {threshold:.4f}")
    print("-" * 45)
    print(f" AUPRC:             {auprc:.4f}")
    print(f" AUROC:             {auroc:.4f}")
    print(f" F1 Macro:          {f1_mac:.4f}")
    print(f" Precision Macro:   {prec_mac:.4f}")
    print(f" Recall Macro:      {rec_mac:.4f}")
    print("="*45)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="评估已保存的 TREPP 模型")
    parser.add_argument("--model", type=str, default="trepp_experiments.pkl", help="pkl模型路径")
    parser.add_argument("--eval_data", type=str, required=True, help="eval特征csv文件路径")
    args = parser.parse_args()
    
    evaluate(args.model, args.eval_data)