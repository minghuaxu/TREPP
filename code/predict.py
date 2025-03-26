import joblib
import numpy as np
from catboost import CatBoostClassifier
import pandas as pd
import os
import logging
from mylib import LoggerMaker


def read_feature(input_file, model_dir, logger):
    def get_stack_x(X, model_dir):
        model_num = 15
        catboost_predictions = np.zeros((len(X), model_num))
        for i in range(model_num):
            model = joblib.load(f"{model_dir}/catboost_model_{i}.bin")
            probs = model.predict_proba(X)
            preds = np.argmax(probs, axis=1)
            catboost_predictions[:, i] = preds
        return catboost_predictions

    random_seed = 2024
    df = pd.read_csv(input_file, sep='\t')
    ids = df['id']
    df['gc_value'] = df['g_value'] + df['c_value']
    df['at_value'] = df['a_value'] + df['t_value']
    # 按照下面的顺序去组织特征的顺序
    feature_columns = ['gc_value', 'at_value', 'motif_len', 'trf_len', 'intron', 'intron_flank500', 'exon', 'exon_flank500', 'strc_1000', 'strc_10000', 'at_1000', 'gc_1000']
    X = df[feature_columns]
    X = X.astype(float)
    X = X.values
    X = get_stack_x(X, model_dir)

    return X, ids
    

if __name__ == '__main__':
    global_logger = LoggerMaker('./', 'predict', 'predict.log', logging.INFO).logger

    model_dir = "/home/xuminghua/STRP3/pictures/ablation_experiment/no_orf/models"
    feature_path = "feature.csv"
    X, ids = read_feature(feature_path, model_dir, global_logger)

    # model = joblib.load('/home/xuminghua/STRP3/pictures/ablation_experiment/no_orf/models/ensemble_model_0321.bin')
    model = joblib.load('/home/xuminghua/STRP3/pictures/ablation_experiment/no_orf/ensemble_model_0324.bin')

    probs = model.predict_proba(X)[: ,1]
    preds = (probs > 0.4).astype(int)
    # preds = np.argmax(probs, axis=1)

    # probs = model.predict_proba(X)[: ,1]

    result_df = pd.DataFrame({'id': ids, 'prob': probs, 'pred': preds})
    result_df.to_csv('result0324.csv', index=False, sep='\t')