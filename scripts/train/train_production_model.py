#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 2: Production Training with JSON Config
============================================
1. 读取 best_params.json 获取特征列表和模型参数。
2. 读取 Train + Eval 数据并合并。
3. 使用固定参数在全量数据上训练 CatBoost。
4. 使用 OOF 预测训练 Meta Learner 并确定阈值。
5. 输出最终模型 trepp_production.pkl。
"""

import argparse
import os
import logging
import time
import json
import joblib
import numpy as np
import pandas as pd

from catboost import CatBoostClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import f1_score, average_precision_score
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression

# 环境变量设置
os.environ['PYTHONHASHSEED'] = '0'
SEED = 2025
np.random.seed(SEED)

def setup_logger(out_dir: str) -> logging.Logger:
    os.makedirs(out_dir, exist_ok=True)
    log_file = os.path.join(out_dir, f"prod_train_{time.strftime('%m%d_%H%M')}.log")
    logger = logging.getLogger("TREPP_PROD")
    logger.setLevel(logging.INFO)
    if logger.handlers: logger.handlers = []
    fh = logging.FileHandler(log_file)
    sh = logging.StreamHandler()
    formatter = logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S')
    fh.setFormatter(formatter); sh.setFormatter(formatter)
    logger.addHandler(fh); logger.addHandler(sh)
    return logger

def load_full_data(train_path, eval_path, feat_cols, logger):
    logger.info("Loading and merging Train + Eval datasets...")
    df_train = pd.read_csv(train_path)
    df_eval = pd.read_csv(eval_path)
    
    # 合并
    df = pd.concat([df_train, df_eval], ignore_index=True)
    y = df["label"].values.astype(int)
    
    # 严格按照 JSON 中的特征顺序提取
    missing = [c for c in feat_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing features in dataset: {missing}")
        
    X = df[feat_cols].values.astype(float)
    logger.info(f"Full Dataset: {X.shape}, Pos: {sum(y)}, Neg: {len(y)-sum(y)}")
    return X, y

def scan_thresholds_macro(y_true, y_probs, logger, prefix=""):
    # thresholds = np.arange(0.1, 0.95, 0.05)
    thresholds = [0.5]
    best_f1 = 0.0
    best_th = 0.5
    for th in thresholds:
        preds = (y_probs >= th).astype(int)
        f1 = f1_score(y_true, preds, average="macro")
        if f1 > best_f1:
            best_f1 = f1
            best_th = th
    logger.info(f"{prefix} Best Threshold: {best_th:.2f} -> Macro F1: {best_f1:.4f}")
    return best_th, best_f1

def train_base_models_fixed(X, y, params_list, logger):
    logger.info(f"\n=== Training {len(params_list)} Base Models (Full Data) ===")
    
    neg_idx = np.where(y == 0)[0]
    pos_idx = np.where(y == 1)[0]
    pos_n = len(pos_idx)
    
    models = []
    oof_preds = np.zeros((len(X), len(params_list)))
    
    for i, params in enumerate(params_list):
        # 确保随机种子的一致性 (使用 params 里自带的 random_state 或外部控制)
        # 这里的 current_seed 主要用于 Bagging 和 CV 的划分
        current_seed = SEED + i 
        
        # --- Bagging 逻辑 (保持和实验一致) ---
        np.random.seed(current_seed)
        np.random.shuffle(neg_idx)
        curr_neg = np.random.choice(neg_idx, size=min(len(neg_idx), pos_n * 2), replace=False)
        curr_idx = np.concatenate([pos_idx, curr_neg])
        # -----------------------------------
        
        # 初始化模型
        base_est = CatBoostClassifier(**params)
        
        # 初始化 CV (用于 CalibratedClassifier 和 OOF 生成)
        cv_strategy = StratifiedKFold(n_splits=5, shuffle=True, random_state=current_seed)
        
        # 1. 训练校准模型 (全量训练)
        calibrated = CalibratedClassifierCV(
            estimator=base_est, method="sigmoid", cv=cv_strategy, n_jobs=1
        )
        calibrated.fit(X[curr_idx], y[curr_idx])
        models.append(calibrated)
        
        # 2. 生成 OOF 预测 (用于 Meta Learner)
        # 在 Bagging 的子集上做 CV 预测
        subset_oof = cross_val_predict(
            calibrated, X[curr_idx], y[curr_idx], 
            cv=cv_strategy, method="predict_proba", n_jobs=1
        )[:, 1]
        
        # 填充 OOF 矩阵
        full_pred = calibrated.predict_proba(X)[:, 1] # 先填满
        full_pred[curr_idx] = subset_oof              # 用 CV 预测覆盖 Bagging 部分
        oof_preds[:, i] = full_pred
        
        logger.info(f"Model {i+1} trained. Full Data AUPRC (approx): {average_precision_score(y, full_pred):.4f}")
        
    return models, oof_preds

def train_meta_learner(oof, y, logger):
    logger.info("\n=== Training Meta Learner (Logistic Regression) ===")
    
    # 1. 使用 CV 确定最佳阈值 (模拟未知数据)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    meta_oof = np.zeros(len(y))
    
    for tr, val in skf.split(oof, y):
        lr = LogisticRegression(penalty="l1", solver="liblinear", C=0.5, class_weight="balanced", random_state=SEED)
        lr.fit(oof[tr], y[tr])
        meta_oof[val] = lr.predict_proba(oof[val])[:, 1]
        
    best_th, best_f1 = scan_thresholds_macro(y, meta_oof, logger, prefix="[CV Scan]")
    
    # 2. 全量训练 Meta Model
    final_meta_model = LogisticRegression(
        penalty="l1", solver="liblinear", C=0.5, class_weight="balanced",
        max_iter=200, random_state=SEED
    )
    # 同样使用校准
    calibrated_meta = CalibratedClassifierCV(estimator=final_meta_model, method="sigmoid", cv=5)
    calibrated_meta.fit(oof, y)
    
    meta_info = {
        "meta_model": calibrated_meta,
        "threshold": float(best_th),
        "cv_f1_macro": float(best_f1)
    }
    return meta_info

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json_config", default="best_params.json", help="包含参数的 JSON 文件")
    parser.add_argument("--train", default="data/processed/train_features.csv")
    parser.add_argument("--eval", default="data/processed/eval_features.csv")
    parser.add_argument("--out_dir", default="production_model_json")
    args = parser.parse_args()

    logger = setup_logger(args.out_dir)

    # 1. 读取配置
    logger.info(f"Reading configuration from {args.json_config}")
    with open(args.json_config, "r") as f:
        config = json.load(f)
    
    feat_cols = config["feat_cols"]
    params_list = config["models_params"]
    
    logger.info(f"Features: {len(feat_cols)}")
    logger.info(f"Models to train: {len(params_list)}")

    # 2. 加载全量数据
    X_full, y_full = load_full_data(args.train, args.eval, feat_cols, logger)

    # 3. 训练基模型
    models, oof = train_base_models_fixed(X_full, y_full, params_list, logger)

    # 4. 训练 Meta Learner
    meta_info = train_meta_learner(oof, y_full, logger)

    # 5. 保存
    final_artifact = {
        "base_models": models,
        "meta_info": meta_info,
        "feat_cols": feat_cols
    }
    
    save_path = os.path.join(args.out_dir, "trepp_production.pkl")
    joblib.dump(final_artifact, save_path)
    
    logger.info("-" * 50)
    logger.info(f"Production Model Saved: {save_path}")
    logger.info(f"Threshold: {meta_info['threshold']}")
    logger.info(f"Expected Macro F1: {meta_info['cv_f1_macro']:.4f}")
    logger.info("-" * 50)

if __name__ == "__main__":
    main()