import joblib
import numpy as np
from catboost import CatBoostClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, RepeatedStratifiedKFold
import optuna
import pandas as pd
import os
import logging
from sklearn.metrics import average_precision_score, f1_score, roc_auc_score, recall_score, precision_score, log_loss, roc_curve, auc, precision_recall_curve, average_precision_score
from sklearn.metrics import log_loss
import xgboost as xgb
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.model_selection import cross_val_score, StratifiedKFold
from imblearn.over_sampling import SMOTE
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import average_precision_score
from catboost import CatBoostClassifier, Pool
from sklearn.model_selection import train_test_split
from imblearn.under_sampling import RandomUnderSampler
import sys
from mylib import calc_score, LoggerMaker, get_xy_drop
from pathlib import Path

def stacking_x(X, model_dir, global_logger):
    model_paths = Path(model_dir).glob("*.bin")
    stacked_x = np.zeros((len(X),len(model_paths)))
    for model_path in model_paths:
        model = joblib.load(model_path)
        preds = np.argmax(probs, axis=1)
        stacked_x[:,i] = preds
    return stacked_x


def tune_logisitic( X, y, model_outdir, n_trials):
   
    def objective(trial):
        params = {
            'C': trial.suggest_float('C', 1e-8, 1000, log=True),
            'penalty': trial.suggest_categorical('penalty', ['l1', 'l2']),
            'solver': trial.suggest_categorical('solver', ['liblinear', 'saga']),
            'max_iter': trial.suggest_int('max_iter', 100, 3000),
            'class_weight': trial.suggest_categorical('class_weight', ['balanced', None])
        }
        # make sure that l1 penalty is only used with liblinear and saga solvers
        if params['penalty'] == 'l1' and params['solver'] not in ['liblinear', 'saga']:
            params['solver'] = 'liblinear'
        
        model = LogisticRegression(**params, verbose=0, random_state=2024)

        cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=2024)
        scores = cross_val_score(model, X, y, cv=cv, scoring='neg_log_loss')
        return scores.mean()

    sampler = optuna.samplers.TPESampler(seed=2024)
    study = optuna.create_study(direction='maximize', sampler=sampler)
    study.optimize(objective, n_trials=n_trials)
    
    best_params = study.best_params.copy()
    best_params.update({
        'random_state': 2024,
        'verbose': 0
    })
    logisitic_model = LogisticRegression(**best_params)
    logisitic_model.fit(X, y)
    joblib.dump(logisitic_model, 'logistic_model.bin')
    
    return logisitic_model


if __name__ == '__main__':
    global_logger = LoggerMaker('./', 'stack_model', 'tune_stack_no_orf.log', logging.INFO).logger

    model_dir = "models"

    for i in range(15):
        model = joblib.load(f"{model_dir}/catboost_model_{i}.bin")
        eval(f'catboost_model_{i}', model, global_logger)
    
    drop_tag = 'no'

    feature_dir = "/home/xuminghua/STRP3/code/strp3_2_hg19/feature0321/train"
    X_train, y_train, ids = get_xy_drop(os.path.join(feature_dir, 'train.csv'), 'train', global_logger, drop_tag)

    X_train = X_train.values    
    y = y_train.values.flatten()
    X = get_stack_x(X_train, model_dir, global_logger)

    from sklearn.ensemble import IsolationForest
    model = IsolationForest(contamination=0.05, random_state=2024)  # contamination 参数控制异常值的比例
    y_pred = model.fit_predict(X) # y_pred 返回1表示正常，-1表示异常
    X = X[y_pred == 1]
    y = y[y_pred == 1]

    # 进行聚类欠采样
    from imblearn.under_sampling import ClusterCentroids
    # cluster_centroids = ClusterCentroids(random_state=2024, sampling_strategy=0.5, voting='soft') # 创建ClusterCentroids对象
    cluster_centroids = ClusterCentroids(random_state=2024)
    X_resampled, y_resampled = cluster_centroids.fit_resample(X, y)

    logisitic_model = train_logisitic(X_resampled, y_resampled, model_dir, 200)
    
    drop_tag = 'no'
    feature_dir = "/home/xuminghua/STRP3/code/strp3_2_hg19/feature0321/test"
    X_test, y_test, ids = get_xy_drop(os.path.join(feature_dir, 'test.csv'), 'train', global_logger, drop_tag)
    
    X_test = X_test.values

    X_test = get_stack_x(X_test, global_logger)

    X_test_copy = np.round(X_test, 2)

    y_test_copy = y_test.astype(int)
    y_test_copy = y_test_copy.values
    print(type(y_test_copy))
    print(y_test_copy.shape)
    # y_test_copy现在是(104,)，要转成(104,1)
    y_test_copy = y_test_copy.reshape(-1, 1)
    print(y_test_copy.shape)
    X_y_test = np.concatenate((y_test_copy, X_test_copy), axis=1)   
    df_test = pd.DataFrame(X_y_test, columns=['label'] + [f'feature_{i}' for i in range(X_test_copy.shape[1])])
    df_test['label'] = df_test['label'].astype(int)
    

    # predictions = rf_model.predict(X_test)
    predictions = ensemble_model.predict_proba(X_test)
    probs = predictions[:, 1]

    df_result = pd.read_csv('/home/xuminghua/STRP3/code/strp3_2_hg19/data/rexprt_input.sorted.bed', sep='\t', header=None)
    df_result['strp3_probs'] = probs
    df_result['strp3_id'] = df_result[0] + '_' + df_result[1].astype(str) + '_' + df_result[2].astype(str) + '_' + df_result[3]
    df_result = df_result[['strp3_id', 'strp3_probs']]
    # df_result['strp3_prob'] = df_result['strp3_prob'].apply(lambda x: round(x, 4))
    df_result.to_csv('strp3_prob.csv', index=False, sep='\t')

    probs = np.round(probs, 2)
    df_test['prob'] = probs

    
    df_test = df_test.sort_values(by='label')
    df_test = df_test[['label', 'prob'] + [f'feature_{i}' for i in range(X_test_copy.shape[1])]]

    roc_auc, pr_auc = calc_score(predictions[:, 1], y_test, './', global_logger, 0.5)

    # joblib.dump(rf_model, 'stack_model.bin')


