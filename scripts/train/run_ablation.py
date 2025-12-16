#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TREPP Ablation Study (V7 Based)
===============================
功能：
1. 定义特征组 (Gene vs Seq)。
2. 依次移除特定特征组，重新训练 TREPP V7 模型。
3. 输出每次消融的性能指标和预测结果 (id, prob)。
"""

import argparse
import os
import sys
import logging
import numpy as np
import pandas as pd
import joblib
from sklearn.metrics import average_precision_score, roc_auc_score, f1_score

# 复用 train_final_v7.py 中的核心逻辑
# 假设 train_final_v7.py 在 scripts/ 目录下，且名为 train_final_v7.py
# 如果文件名不同，请修改下面的 import
try:
    from trepp_v7_macro_2 import (
        load_data, 
        setup_logger, 
        select_features, 
        train_base_models, 
        train_meta_learner,
        SEED
    )
except ImportError:
    print("Error: Could not import trepp_v7_macro_2.py. Please ensure it is in the same directory or PYTHONPATH.")
    sys.exit(1)

# ==========================================
# 1. 定义特征分组逻辑
# ==========================================
def get_feature_groups(all_features: list) -> dict:
    """
    根据特征名定义分组。
    需要根据你实际的列名逻辑进行调整。
    """
    groups = {
        "All": [], # 占位，表示所有特征
        "No_Gene": [],
        "No_Seq": []
    }
    
    # 定义关键词规则
    # Gene Features: 距离、位置、基因属性、GnomAD评分
    gene_keywords = [
        'dist_', 'gene_', 'tss', 'exon', 'cds', 'utr', 'promoter', 
        'pLI', 'oe_lof', 'in_gene', 'in_exon', 'in_cds', 'delta_'
    ]
    
    # Seq Features: 序列长度、GC含量、Motif属性、重复单元、TRF密度
    seq_keywords = [
        'tr_len', 'motif', 'gc', 'at', 'repeat', 'str_density', 
        'trinuc', 'pentanuc', 'delta_'
    ]
    
    # 生成 No_Gene 列表 (只保留 Seq 特征)
    no_gene_feats = []
    for f in all_features:
        # 如果特征名包含任何 gene 关键词，剔除
        is_gene = any(k in f for k in gene_keywords)
        if not is_gene:
            no_gene_feats.append(f)
    groups["No_Gene"] = no_gene_feats
    
    # 生成 No_Seq 列表 (只保留 Gene 特征)
    no_seq_feats = []
    for f in all_features:
        is_seq = any(k in f for k in seq_keywords)
        if not is_seq:
            no_seq_feats.append(f)
    groups["No_Seq"] = no_seq_feats
    
    groups["All"] = all_features
    
    return groups

# ==========================================
# 2. 核心消融流程
# ==========================================
def run_ablation(args):
    os.makedirs(args.out_dir, exist_ok=True)
    logger = setup_logger(args.out_dir)
    logger.info("=== Starting Ablation Study ===")
    
    # 1. 加载全量数据
    # 注意：load_data 内部已经做了 add_engineered_features
    # 返回的 X 是 numpy 数组，feat_cols 是列名列表
    # 我们这里需要稍微修改一下调用逻辑，因为我们需要根据列名筛选
    
    logger.info(f"Loading data from {args.train} ...")
    df_train = pd.read_csv(args.train)
    
    df_eval = pd.read_csv(args.eval)
    
    # 构造 ID 列 (用于输出)
    for df in [df_train, df_eval]:
        if 'id' not in df.columns:
            df['id'] = df['chr'].astype(str) + '_' + df['start'].astype(str) + '_' + df['end'].astype(str) + '_' + df['motif'].astype(str)
            
    # 提取标签
    y_train = df_train['label'].values.astype(int)
    y_eval = df_eval['label'].values.astype(int)
    
    # 提取特征名
    meta_cols = ["chr", "start", "end", "motif", "gene_name", "label", "id", "oe_lof", "pLI", "tr_gc", "motif_gc"]
    all_feat_cols = sorted([c for c in df_train.columns if c not in meta_cols])
    
    # 定义分组
    feature_groups = get_feature_groups(all_feat_cols)
    
    summary_results = []
    
    # 2. 循环实验
    for group_name, feats_to_use in feature_groups.items():
        logger.info(f"\n{'='*40}")
        logger.info(f"Running Ablation: {group_name}")
        logger.info(f"Features count: {len(feats_to_use)}")
        logger.info(f"{'='*40}")
        
        if len(feats_to_use) == 0:
            logger.warning(f"Skipping {group_name}: No features selected.")
            continue
            
        # 准备数据矩阵
        X_tr = df_train[feats_to_use].values.astype(float)
        X_ev = df_eval[feats_to_use].values.astype(float)
        
        # 特征筛选 (可选，但建议保留以保持与主模型一致)
        # 注意：这里 select_features 返回的是索引，是相对于 feats_to_use 的索引
        if args.use_feature_selection:
            keep_indices, selected_names = select_features(
                X_tr, y_train, feats_to_use, logger, top_k=min(args.top_k, len(feats_to_use))
            )
            X_tr = X_tr[:, keep_indices]
            X_ev = X_ev[:, keep_indices]
            logger.info(f"Features after selection: {X_tr.shape[1]}")
        
        # 训练基模型
        models, oof, scores = train_base_models(
            X_tr, y_train, args.n_models, args.n_trials, args.out_dir, logger
        )
        
        # 训练元模型
        meta_info = train_meta_learner(
            oof, y_train, scores, args.out_dir, logger
        )
        
        # 预测评估集
        # 1. 基模型预测
        n_samples = len(X_ev)
        base_preds = np.zeros((n_samples, len(models)))
        for i, m in enumerate(models):
            base_preds[:, i] = m.predict_proba(X_ev)[:, 1]
            
        # 2. 元模型预测
        sel_idx = meta_info["selected_idx"]
        meta_model = meta_info["meta_model"]
        train_thresh = meta_info["threshold"]
        
        final_probs = meta_model.predict_proba(base_preds[:, sel_idx])[:, 1]
        
        # 计算指标
        auprc = average_precision_score(y_eval, final_probs)
        auroc = roc_auc_score(y_eval, final_probs)
        
        # 使用训练集阈值
        pred_labels = (final_probs >= 0.5).astype(int)
        f1 = f1_score(y_eval, pred_labels, average='macro')
        
        logger.info(f"[{group_name}] Result: AUPRC={auprc:.4f}, AUROC={auroc:.4f}, MacroF1={f1:.4f}")
        
        summary_results.append({
            "Experiment": group_name,
            "AUPRC": auprc,
            "AUROC": auroc,
            "MacroF1": f1,
            "Threshold": train_thresh,
            "Num_Features": X_tr.shape[1]
        })
        
        # 输出预测结果 (id \t prob)
        out_file = os.path.join(args.out_dir, f"preds_{group_name}.txt")
        result_df = pd.DataFrame({
            "id": df_eval['id'],
            "prob": final_probs
        })
        # 保存为制表符分隔，无 header (或者根据需求加 header)
        result_df.to_csv(out_file, sep='\t', index=False, header=True)
        logger.info(f"Predictions saved to {out_file}")

    # 保存汇总结果
    summary_df = pd.DataFrame(summary_results)
    summary_file = os.path.join(args.out_dir, "ablation_summary.csv")
    summary_df.to_csv(summary_file, index=False)
    logger.info(f"\nAll ablation experiments completed. Summary saved to {summary_file}")
    print(summary_df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", required=True, help="Training CSV")
    parser.add_argument("--eval", required=True, help="Evaluation CSV")
    parser.add_argument("--out_dir", required=True, help="Output directory")
    
    # 复用主模型的参数
    parser.add_argument("--n_models", type=int, default=8)
    parser.add_argument("--n_trials", type=int, default=150) # 可以设小一点加速消融实验
    parser.add_argument("--top_k", type=int, default=30)
    parser.add_argument("--use_feature_selection", action="store_true", default=True)
    
    args = parser.parse_args()
    
    run_ablation(args)