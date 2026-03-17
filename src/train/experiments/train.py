"""
快捷复现实验结果。
"""
import argparse
import os
import logging
import time
import joblib
import numpy as np
import pandas as pd

from catboost import CatBoostClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import (
    average_precision_score,
    roc_auc_score,
    f1_score,
    precision_score,
    recall_score,
)
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression

# 确保哈希稳定
os.environ['PYTHONHASHSEED'] = '0'
SEED = 2025
np.random.seed(SEED)

def setup_logger(out_dir: str) -> logging.Logger:
    os.makedirs(out_dir, exist_ok=True)
    log_file = os.path.join(out_dir, f"repro_check_{time.strftime('%m%d_%H%M')}.log")
    logger = logging.getLogger("TREPP_REPRO")
    logger.setLevel(logging.INFO)
    if logger.handlers: logger.handlers = []
    fh = logging.FileHandler(log_file); sh = logging.StreamHandler()
    formatter = logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S')
    fh.setFormatter(formatter); sh.setFormatter(formatter)
    logger.addHandler(fh); logger.addHandler(sh)
    return logger

def extract_artifacts(pkl_path: str, logger: logging.Logger):
    """从 pickle 文件中提取特征列表和基模型参数"""
    logger.info(f"Loading artifacts from: {pkl_path}")
    if not os.path.exists(pkl_path):
        raise FileNotFoundError(f"File not found: {pkl_path}")
        
    artifacts = joblib.load(pkl_path)
    base_models = artifacts["base_models"]
    feat_cols = artifacts["feat_cols"]
    
    # 提取每个基模型的参数
    model_params_list = []
    for i, model in enumerate(base_models):
        params = model.estimator.get_params()
        model_params_list.append(params)
        
    logger.info(f"Successfully extracted {len(model_params_list)} model parameter sets.")
    logger.info(f"Feature count locked: {len(feat_cols)}")
    return feat_cols, model_params_list

def load_data_fixed_features(csv_path: str, feat_cols: list, logger: logging.Logger):
    """只加载指定的特征列，顺序严格一致"""
    logger.info(f"Loading {csv_path}...")
    df = pd.read_csv(csv_path)
    y = df["label"].values.astype(int)
    
    # 检查缺失列
    missing = [c for c in feat_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing features in {csv_path}: {missing}")
        
    # 强制排序
    X = df[feat_cols].values.astype(float)
    return X, y

def train_fixed_base_models(X, y, params_list, model_dir, logger):
    """使用提取的参数列表进行训练 (无 Optuna)"""
    logger.info(f"\n=== Training {len(params_list)} Base Models (Fixed Params) ===")

    neg_idx = np.where(y == 0)[0]
    pos_idx = np.where(y == 1)[0]
    pos_n = len(pos_idx)

    models = []
    oof_preds = np.zeros((len(X), len(params_list)))
    scores = []

    for i, params in enumerate(params_list):
        current_seed = SEED + i
        
        np.random.seed(current_seed)
        np.random.shuffle(neg_idx)
        curr_neg = np.random.choice(
            neg_idx, size=min(len(neg_idx), pos_n * 2), replace=False
        )
        curr_idx = np.concatenate([pos_idx, curr_neg])
        
        if params.get('random_state') != current_seed and params.get('random_seed') != current_seed:
            logger.warning(f"Model {i} param seed {params.get('random_state')} != loop seed {current_seed}. Overwriting to ensure reproduction.")
            params['random_state'] = current_seed
            if 'random_seed' in params: del params['random_seed']

        base_est = CatBoostClassifier(**params)

        cv_strategy = StratifiedKFold(n_splits=5, shuffle=True, random_state=current_seed)

        calibrated = CalibratedClassifierCV(
            estimator=base_est, method="sigmoid", cv=cv_strategy, n_jobs=1
        )
        
        calibrated.fit(X[curr_idx], y[curr_idx])

        subset_oof = cross_val_predict(
            calibrated,
            X[curr_idx],
            y[curr_idx],
            cv=cv_strategy,
            method="predict_proba",
            n_jobs=1,
        )[:, 1]

        full_pred = calibrated.predict_proba(X)[:, 1]
        full_pred[curr_idx] = subset_oof

        oof_preds[:, i] = full_pred
        models.append(calibrated)

        score = average_precision_score(y, full_pred)
        scores.append(score)
        logger.info(f"Model {i + 1} (Fixed) OOF AUPRC: {score:.4f}")

    return models, oof_preds, np.array(scores)

def train_meta_learner_repro(oof, y, logger):
    logger.info("\n=== Meta Learner Training (Reproduction) ===")
    
    X_meta = oof
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    meta_oof = np.zeros(len(y))

    for fold, (tr_idx, val_idx) in enumerate(skf.split(X_meta, y), start=1):
        meta_model = LogisticRegression(
            penalty="l1", solver="liblinear", C=0.5, class_weight="balanced", 
            max_iter=200, random_state=SEED + fold
        )
        meta_model.fit(X_meta[tr_idx], y[tr_idx])
        meta_oof[val_idx] = meta_model.predict_proba(X_meta[val_idx])[:, 1]
    
    final_meta_model = LogisticRegression(
        penalty="l1", solver="liblinear", C=0.3, class_weight="balanced",
        max_iter=200, random_state=SEED
    )
    calibrated = CalibratedClassifierCV(estimator=final_meta_model, method="sigmoid", cv=3, n_jobs=1)
    calibrated.fit(X_meta, y)
    
    return calibrated

def evaluate_repro(X_eval, y_eval, models, meta_model, threshold, logger):
    logger.info("\n=== Final Evaluation on Eval Set (Reproduction) ===")
    
    # Base Models
    base_preds = np.zeros((len(X_eval), len(models)))
    for i, m in enumerate(models):
        base_preds[:, i] = m.predict_proba(X_eval)[:, 1]
        
    # Meta Model
    final_prob = meta_model.predict_proba(base_preds)[:, 1]
    
    # Metrics
    auprc = average_precision_score(y_eval, final_prob)
    auroc = roc_auc_score(y_eval, final_prob)
    
    pred_label = (final_prob >= threshold).astype(int)
    f1_macro = f1_score(y_eval, pred_label, average="macro")
    
    logger.info("-" * 40)
    logger.info(f"REPRODUCED METRICS:")
    logger.info(f"AUPRC (eval)    : {auprc:.4f}")
    logger.info(f"AUROC (eval)    : {auroc:.4f}")
    logger.info(f"Macro F1 (th={threshold:.2f}): {f1_macro:.4f}")
    logger.info("-" * 40)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pkl_path", default="trepp_experiments.pkl")
    parser.add_argument("--train", default="../../data/processed/train_features.csv")
    parser.add_argument("--eval", default="../../data/processed/eval_features.csv")
    parser.add_argument("--out_dir", default="results/repro_check")
    args = parser.parse_args()

    logger = setup_logger(args.out_dir)

    # 1. 提取参数
    feat_cols, model_params_list = extract_artifacts(args.pkl_path, logger)

    # 2. 加载数据 (使用提取出来的 feature 列，确保特征对齐)
    X_train, y_train = load_data_fixed_features(args.train, feat_cols, logger)
    X_eval, y_eval = load_data_fixed_features(args.eval, feat_cols, logger)
    
    logger.info(f"Train Shape: {X_train.shape}")
    logger.info(f"Eval Shape : {X_eval.shape}")

    # 3. 复现训练 (Train Base Models)
    models, oof, scores = train_fixed_base_models(
        X_train, y_train, model_params_list, args.out_dir, logger
    )

    # 4. 复现 Meta Learner
    meta_model = train_meta_learner_repro(oof, y_train, logger)
    
    # 5. 验证结果
    threshold = 0.5
    evaluate_repro(X_eval, y_eval, models, meta_model, threshold, logger)

if __name__ == "__main__":
    main()