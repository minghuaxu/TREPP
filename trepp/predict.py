import joblib
import numpy as np
from catboost import CatBoostClassifier
import pandas as pd
import os
import logging
import argparse

class LoggerMaker:
    def __init__(self, log_name, out_dir='./',  log_file=None, level=logging.INFO):
        self.out_dir = out_dir
        self.logger = None
        self.log_file = log_file
        self.log_name = log_name
        self.level = level
        self._make_logger()
    
    def _make_logger(self):
        logger = logging.getLogger(self.log_name)
        logger.setLevel(self.level)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

        ch = logging.StreamHandler()
        ch.setLevel(self.level)
        ch.setFormatter(formatter)
        logger.addHandler(ch)

        if self.log_file:
            fh = logging.FileHandler(os.path.join(self.out_dir, self.log_file))
            fh.setLevel(self.level)
            fh.setFormatter(formatter)
            logger.addHandler(fh)

        self.logger = logger


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
    feature_columns = ['label', 'gc_value', 'at_value', 'motif_len', 'trf_len',  'strc_1000', 'strc_10000', 'at_1000', 'gc_1000',  'intron', 'intron_flank500', 'exon', 'exon_flank500', 'UTR3', 'UTR3_flank500', 'UTR5', 'UTR5_flank500', 'start_codon', 'start_codon_flank100', 'stop_codon', 'stop_codon_flank100']
    X = df[feature_columns]
    X = X.astype(float)
    X = X.values
    X = get_stack_x(X, model_dir)

    return X, ids
    


def main(args):
    global_logger = LoggerMaker('trepp_test', level=logging.DEBUG).logger

    from feature_tools.feature_main import get_feature
    
    # get_feature('feature_tools', 
    #     args.input_file,
    #     args.db_file, 
    #     args.reference, 
    #     args.outdir, 
    #     args.suffix, 
    #     global_logger, 
    #     resume=False,
    #     reference_name='hg19',
    #     threads=args.threads)

    model_dir = args.model_dir
    feature_path = os.path.join(args.outdir, args.suffix + '.csv')
    X, ids = read_feature(feature_path, model_dir, global_logger)

    model = joblib.load(os.path.join(model_dir, 'ensemble_model.bin'))

    probs = model.predict_proba(X)[: ,1]
    preds = (probs > 0.5).astype(int)

    result_df = pd.DataFrame({'id': ids, 'prob': probs, 'pred': preds})
    result_df.to_csv(os.path.join(args.outdir, args.suffix + '_pred.csv'), index=False, sep='\t')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="TREPP: Tandem Repeat Expansion Pathogenicity Prediction).",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('-i', '--input_file', required=True, help='STR loci site in standard BED format with cloumns: [chrom, start, end, motif]')
    parser.add_argument('-o', '--outdir', required=True, help='Output directory')
    parser.add_argument('-m', '--model_dir', required=True, help='Trained model directory')
    parser.add_argument('-r', '--reference', required=True, help='Reference genome file')
    parser.add_argument('--db_file', required=True, help='Database file for genome annotation')

    parser.add_argument('--suffix', default='trepp', help='Suffix for output files (default: trepp)')
    parser.add_argument("--log_level", type=str, default="INFO", help="Set the logging level (e.g., DEBUG, INFO, WARNING)")
    parser.add_argument("--log_filename", type=str, default="trepp", help="Set the log output file name (default: trepp)")
    parser.add_argument("--refeature", type=str, default="True", help="Whether to re-extract features (default: True)")
    parser.add_argument('-t', '--threads', type=int, default=10, help='Number of threads (default: 10)')

    args = parser.parse_args()

    main(args)

    