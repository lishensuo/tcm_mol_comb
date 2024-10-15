import os
import time
import sys
import pandas as pd
import numpy as np
from autogluon.tabular import TabularDataset, TabularPredictor

mode = "medium_quality"
metricses = ["ZIP","HSA","Loewe","Bliss"]

fp_types = ["MACCS","PubChem","Substructure"]
# gene_rank = [958]
# metrics = metricses[0]
# fp_type = fp_types[0]

res_test_list = []

for metrics in metricses:
    dat_train = pd.read_csv("./ML/data/test_bc/train_dat_"+metrics+"_Pan.csv")

    for fp_type in fp_types:
        print("################## Start ", metrics, "————", fp_type)
        selected_columns = [col for col in dat_train.columns if col.startswith(('Exp','Desc',fp_type))]
        # selected_columns = [col for col in dat_train.columns if col.startswith(fp_type)]
        selected_columns = ['set','label'] + selected_columns

        dat_train_fp = dat_train[selected_columns]

        dat_cv = dat_train_fp.loc[dat_train_fp.set=='train',:].drop(["set"], axis=1)
        dat_test = dat_train_fp.loc[dat_train_fp.set!='train',:].drop(["set"], axis=1)

        dat_cv_Tab = TabularDataset(dat_cv)
        dat_test_Tab = TabularDataset(dat_test)

        save_path = "./ML/data/test_bc/model08/" + metrics + "_"+ mode + "/" + fp_type
 
        metric = "root_mean_squared_error"
        eval_metric = ["r2"]	

        predictor = TabularPredictor(label="label", path=save_path, 
            verbosity=1).fit(dat_cv_Tab,presets=mode, 
                                ag_args_fit={'num_gpus': 2})
        
        predictor.leaderboard()
        res_test = predictor.leaderboard(dat_test_Tab, extra_metrics = eval_metric,silent=True)
        eval_columns = ['model'] + ['score_val','score_test'] + eval_metric
        res_test = res_test.loc[:,eval_columns]
        res_test['Metrics'] = metrics
        res_test['FP'] = fp_type

        res_test_list.append(res_test)

res_test_df = pd.concat(res_test_list, axis=0)

file_name_2 = "./ML/data/test_bc/model08/Model_test_result.csv"

res_test_df.to_csv(file_name_2, index=False)
