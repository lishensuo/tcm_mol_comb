"""
Autogluon CV Related Issues 

- https://github.com/autogluon/autogluon/discussions/3675
- https://github.com/autogluon/autogluon/blob/cb5098137b3bf604b4a800326533e597614a1507/core/src/autogluon/core/utils/utils.py#L185
"""

import pandas as pd
import numpy as np
import logging
from autogluon.tabular import TabularDataset, TabularPredictor

# Configuration
MODE = "medium_quality"
METRICS_LIST = ["ZIP", "HSA", "Loewe", "Bliss"]
FP_TYPES = ["MACCS", "PubChem", "Substructure"]
SAVE_PATH_BASE = "./ML/model/"
TRAIN_DATA_PATH = "./ML/data/train_dat_{}.csv"
OUTPUT_FILE = "./ML/model/Model_test_result.csv"

# Logging setup
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_data(metrics):
    """Load training data for the given metrics."""
    file_path = TRAIN_DATA_PATH.format(metrics)
    logger.info(f"Loading data from {file_path}")
    return pd.read_csv(file_path)

def preprocess_data(dat_train, fp_type):
    """Preprocess data by selecting relevant columns and splitting into CV and test sets."""
    logger.info(f"Preprocessing data for fingerprint type: {fp_type}")
    selected_columns = [col for col in dat_train.columns if col.startswith(fp_type)]
    selected_columns = ['set', 'label'] + selected_columns
    dat_train_fp = dat_train[selected_columns]
    
    dat_cv = dat_train_fp.loc[dat_train_fp.set == 'train', :].drop(["set"], axis=1)
    dat_test = dat_train_fp.loc[dat_train_fp.set != 'train', :].drop(["set"], axis=1)
    
    return TabularDataset(dat_cv), TabularDataset(dat_test)

def train_and_evaluate(dat_cv, dat_test, metrics, fp_type):
    """Train a model and evaluate it on the test set."""
    save_path = f"{SAVE_PATH_BASE}{metrics}_{MODE}/{fp_type}"
    eval_metric = ['root_mean_squared_error', 'mean_squared_error', 'mean_absolute_error', 
                   'median_absolute_error', 'mean_absolute_percentage_error', 'r2']
    
    logger.info(f"Training model for {metrics} with {fp_type} fingerprints")
    predictor = TabularPredictor(label="label", path=save_path, verbosity=2).fit(
        dat_cv, presets=MODE, ag_args_fit={'num_gpus': 2}, num_bag_folds=5
    )
    
    logger.info("Generating leaderboard and evaluating on test set")
    res_test = predictor.leaderboard(dat_test, extra_metrics=eval_metric, silent=True)
    eval_columns = ['model'] + ['score_val', 'score_test'] + eval_metric
    res_test = res_test.loc[:, eval_columns]
    res_test['Metrics'] = metrics
    res_test['FP'] = fp_type
    
    return res_test

def main():
    res_test_list = []
    
    for metrics in METRICS_LIST:
        dat_train = load_data(metrics)
        
        for fp_type in FP_TYPES:
            logger.info(f"Starting evaluation for {metrics} — {fp_type}")
            dat_cv, dat_test = preprocess_data(dat_train, fp_type)
            res_test = train_and_evaluate(dat_cv, dat_test, metrics, fp_type)
            res_test_list.append(res_test)
    
    res_test_df = pd.concat(res_test_list, axis=0)
    res_test_df.to_csv(OUTPUT_FILE, index=False)
    logger.info(f"Results saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()