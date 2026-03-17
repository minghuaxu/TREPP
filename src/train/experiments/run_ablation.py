#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import os
import sys
import logging
from optuna.samplers import TPESampler  # [修改点1] 引入采样器
import numpy as np
import pandas as pd
import joblib
from sklearn.metrics import average_precision_score, roc_auc_score, f1_score
from catboost import CatBoostClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import (
    average_precision_score,
    roc_auc_score,
    f1_score,
    precision_score,
    recall_score,
)

import optuna
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression
import time

SEED =2025

def setup_logger(out_dir: str) -> logging.Logger:
    os.makedirs(out_dir, exist_ok=True)
    log_file = os.path.join(out_dir, f"training_v7_macro_{time.strftime('%m%d_%H%M')}.log")
    logger = logging.getLogger("TREPP_V7")
    logger.setLevel(logging.INFO)

    if logger.handlers:
        logger.handlers = []

    fh = logging.FileHandler(log_file)
    sh = logging.StreamHandler()
    formatter = logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S')
    fh.setFormatter(formatter)
    sh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger

def load_data(csv_path: str, logger: logging.Logger):
    logger.info(f"Loading {csv_path}...")
    df = pd.read_csv(csv_path)
    y = df["label"].values.astype(int)

    df['id'] = df['chr'].astype(str) + '_' + df['start'].astype(str) + '_' + df['end'].astype(str) + '_' + df['motif'].astype(str)

    meta_cols = ["chr", "start", "end", "motif", "gene_name", "label", "id", "oe_lof", "pLI", "tr_gc", "motif_gc"]
    # [修改点3] 排序特征列，防止因读取顺序不同导致的特征错位
    feat_cols = sorted([c for c in df.columns if c not in meta_cols])

    X = df[feat_cols].values.astype(float)
    return X, y, feat_cols

def select_features(X: np.ndarray, y: np.ndarray, feature_names, logger: logging.Logger, top_k: int = 50):
    """
    使用一个快速的 CatBoost 模型筛选特征，保留重要性最高的 Top-K 个特征。
    """
    logger.info(f"\n=== Feature Selection (Keeping Top {top_k}) ===")
    
    # 训练一个快速模型用于评估重要性
    # 这里的参数可以设得简单一点，速度优先
    quick_model = CatBoostClassifier(
        iterations=500,
        learning_rate=0.1,
        depth=6,
        random_seed=SEED,
        verbose=False,
        allow_writing_files=False,
        auto_class_weights="Balanced"
    )
    
    quick_model.fit(X, y)
    
    # 获取特征重要性
    importances = quick_model.get_feature_importance()
    
    # 创建 DataFrame 方便排序
    df_imp = pd.DataFrame({"feature": feature_names, "importance": importances})
    df_imp = df_imp.sort_values("importance", ascending=False)
    
    # 打印前 10 个最重要的特征
    logger.info("Top 10 Important Features:")
    for i, row in df_imp.head(20).iterrows():
        logger.info(f"  {row['feature']:<25} : {row['importance']:.4f}")
        
    # 选出 Top K 的特征名称
    selected_feat_names = set(df_imp.head(top_k)["feature"].values)
    
    # 找到这些特征对应的原始索引列号
    selected_indices = [i for i, f in enumerate(feature_names) if f in selected_feat_names]
    selected_names_sorted = [feature_names[i] for i in selected_indices]
    
    logger.info(f"Dropped {len(feature_names) - len(selected_indices)} features. Remaining: {len(selected_indices)}")
    
    return np.array(selected_indices), selected_names_sorted


def tune_catboost(X, y, idx: int, n_trials: int, random_seed, logger: logging.Logger):
    logger.info(f"Tuning Base Model {idx + 1} with seed {random_seed}...")

    pos_idx = np.where(y == 1)[0]
    neg_idx = np.where(y == 0)[0]

    np.random.seed(random_seed)
    np.random.shuffle(neg_idx)
    sample_neg = neg_idx[: len(pos_idx) * 2]

    curr_X = np.vstack([X[pos_idx], X[sample_neg]])
    curr_y = np.hstack([y[pos_idx], y[sample_neg]])

    def obj(trial):
        params = {
            "iterations": 2000,
            "learning_rate": trial.suggest_float("learning_rate", 0.001, 0.5, log=True),
            "depth": trial.suggest_int("depth", 2, 10),
            "l2_leaf_reg": trial.suggest_float("l2_leaf_reg", 2, 200, log=True),
            "verbose": False,
            "allow_writing_files": False,
            "rsm": trial.suggest_float("rsm", 0.5, 1.0),
            "grow_policy": "Lossguide",
            "eval_metric": "Precision",
            "loss_function": "Logloss",
            "auto_class_weights": "Balanced",
            "random_state": random_seed, # CatBoost 内部种子
            "thread_count": -1, # [可选] 限制线程数有助于减少并行计算浮点误差带来的微小差异
        }

        # 这里的 SKF 需要固定随机种子
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=random_seed)
        scores = []
        for tr, val in skf.split(curr_X, curr_y):
            model = CatBoostClassifier(**params)
            model.fit(
                curr_X[tr],
                curr_y[tr],
                eval_set=(curr_X[val], curr_y[val]),
                early_stopping_rounds=20,
                verbose=False,
            )
            preds = model.predict_proba(curr_X[val])[:, 1]
            if np.isnan(preds).any():
                return 0.0
            scores.append(average_precision_score(curr_y[val], preds))
        return float(np.mean(scores))

    optuna.logging.set_verbosity(optuna.logging.WARNING)
    
    # [修改点1] 关键：Optuna 必须使用带 Seed 的 Sampler 才能复现搜索结果
    sampler = TPESampler(seed=random_seed) 
    study = optuna.create_study(direction="maximize", sampler=sampler)
    
    study.optimize(obj, n_trials=n_trials, show_progress_bar=False)

    bp = study.best_params
    bp.update(
        {
            "loss_function": "Logloss",
            "verbose": False,
            "allow_writing_files": False,
            "random_state": random_seed,
        }
    )
    return bp


def train_base_models(
    X: np.ndarray,
    y: np.ndarray,
    n_models: int,
    n_trials: int,
    model_dir: str,
    logger: logging.Logger,
):
    logger.info(f"\n=== Training {n_models} Base Models (CatBoost) ===")

    neg_idx = np.where(y == 0)[0]
    pos_idx = np.where(y == 1)[0]
    pos_n = len(pos_idx)

    models = []
    oof_preds = np.zeros((len(X), n_models))
    scores = []

    for i in range(n_models):
        current_seed = SEED + i
        np.random.seed(current_seed)
        np.random.shuffle(neg_idx)
        curr_neg = np.random.choice(
            neg_idx, size=min(len(neg_idx), pos_n * 2), replace=False
        )
        curr_idx = np.concatenate([pos_idx, curr_neg])

        params = tune_catboost(X, y, i, n_trials, current_seed, logger)

        base_est = CatBoostClassifier(**params)

        # [修改点2] 创建固定的 CV 对象，而不是只传数字 cv=5
        # 这样能保证 Platt Scaling 内部 fold 的划分每次都完全一样
        cv_strategy = StratifiedKFold(n_splits=5, shuffle=True, random_state=current_seed)

        calibrated = CalibratedClassifierCV(
            estimator=base_est, method="sigmoid", cv=cv_strategy, n_jobs=1
        )
        calibrated.fit(X[curr_idx], y[curr_idx])

        subset_oof = cross_val_predict(
            calibrated,
            X[curr_idx],
            y[curr_idx],
            cv=cv_strategy, # 同样使用固定的 CV 对象
            method="predict_proba",
            n_jobs=1,
        )[:, 1]

        full_pred = calibrated.predict_proba(X)[:, 1]
        full_pred[curr_idx] = subset_oof

        oof_preds[:, i] = full_pred
        models.append(calibrated)

        score = average_precision_score(y, full_pred)
        scores.append(score)
        logger.info(f"Model {i + 1} OOF AUPRC: {score:.4f}")

        joblib.dump(calibrated, os.path.join(model_dir, f"base_{i}.pkl"))

    return models, oof_preds, np.array(scores)

def scan_thresholds_macro(
    y_true: np.ndarray,
    y_probs: np.ndarray,
    logger: logging.Logger = None,
    prefix: str = "",
    thresholds: np.ndarray = None,
):
    if thresholds is None:
        thresholds = np.arange(0.1, 0.95, 0.05)

    best_f1_macro = 0.0
    best_th = 0.5

    if logger:
        logger.info(f"\n--- {prefix} Threshold Scan (Macro F1) ---")
        logger.info(f"{'Thresh':<10} | {'F1_pos':<10} | {'F1_neg':<10} | {'F1_macro':<10} | {'Prec':<10} | {'Recall':<10}")
        logger.info("-" * 80)

    for th in thresholds:
        preds = (y_probs >= th).astype(int)
        f1_pos = f1_score(y_true, preds, zero_division=0, average="binary")
        f1_neg = f1_score(1 - y_true, 1 - preds, zero_division=0, average="binary")
        f1_macro = 0.5 * (f1_pos + f1_neg)

        prec = precision_score(y_true, preds, zero_division=0, average="macro")
        rec = recall_score(y_true, preds, zero_division=0, average="macro")

        if f1_macro > best_f1_macro:
            best_f1_macro = f1_macro
            best_th = th

        if logger:
            logger.info(
                f"{th:.2f}       | {f1_pos:.4f}   | {f1_neg:.4f}   | {f1_macro:.4f}   | {prec:.4f}   | {rec:.4f}"
            )

    if logger:
        logger.info("-" * 80)
        logger.info(
            f"Best Threshold (F1_macro): {best_th:.2f} -> Max F1_macro: {best_f1_macro:.4f}"
        )

    return best_th, best_f1_macro

def train_meta_learner(
    oof: np.ndarray,
    y: np.ndarray,
    scores: np.ndarray,
    model_dir: str,
    logger: logging.Logger,
):
    logger.info("\n=== Meta Learner Training (Logistic Regression + Macro F1) ===")
    logger.info("Pruning disabled. Using ALL base models in stacking.")
    
    selected_idx = np.arange(oof.shape[1])
    X_meta = oof[:, selected_idx]

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    meta_oof = np.zeros(len(y))

    for fold, (tr_idx, val_idx) in enumerate(skf.split(X_meta, y), start=1):
        # logger.info(f"Meta fold {fold} / 5 ...")
        meta_model = LogisticRegression(
            penalty="l1",
            solver="liblinear",
            C=0.5,
            class_weight="balanced",
            max_iter=200,
            random_state=SEED + fold, # LR 自身也是确定的
        )
        meta_model.fit(X_meta[tr_idx], y[tr_idx])
        meta_oof[val_idx] = meta_model.predict_proba(X_meta[val_idx])[:, 1]

    best_th, best_f1_macro = scan_thresholds_macro(
        y, meta_oof, logger, prefix="Train Set (Meta OOF)"
    )

    final_meta_model = LogisticRegression(
        penalty="l1",
        solver="liblinear",
        C=0.3,
        class_weight="balanced",
        max_iter=200,
        random_state=SEED,
    )

    calibrated = CalibratedClassifierCV(
        estimator=final_meta_model, method="sigmoid", cv=3, n_jobs=1
    )
    calibrated.fit(X_meta, y)
    
    final_meta_model.fit(X_meta, y)
    final_meta_model = calibrated

    meta_info = {
        "selected_idx": selected_idx,
        "meta_model": final_meta_model,
        "threshold": float(best_th),
        "best_f1_macro_train": float(best_f1_macro),
    }
    joblib.dump(meta_info, os.path.join(model_dir, "meta_info_v7.pkl"))

    with open(os.path.join(model_dir, "best_threshold_macro_v7.txt"), "w") as f:
        f.write(f"{best_th:.6f}\n")

    logger.info(
        f"Meta training done. Best train macro F1 = {best_f1_macro:.4f} @ threshold = {best_th:.2f}"
    )
    return meta_info



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