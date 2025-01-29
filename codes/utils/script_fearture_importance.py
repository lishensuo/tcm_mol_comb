# Calculate feature importance

import pandas as pd
import numpy as np
from autogluon.tabular import TabularDataset, TabularPredictor

mode = "medium_quality"
metricses = ["ZIP"]
fp_types = [
    "MACCS","PubChem","Substructure"
]

res_fi_list = []

for metrics in metricses:
    dat_train = pd.read_csv("./ML/data/train_dat_"+metrics+".csv")

    for fp_type in fp_types:
        print("################## Start ", metrics, "————", fp_type)
        selected_columns = [col for col in dat_train.columns if col.startswith(fp_type)]
        selected_columns = ['set','label'] + selected_columns

        dat_train_fp = dat_train[selected_columns]

        dat_cv = dat_train_fp.loc[dat_train_fp.set=='train',:].drop(["set"], axis=1)
        dat_test = dat_train_fp.loc[dat_train_fp.set!='train',:].drop(["set"], axis=1)

        dat_cv_Tab = TabularDataset(dat_cv)
        dat_test_Tab = TabularDataset(dat_test)

        save_path = "./ML/model/" + metrics + "_"+ mode + "/" + fp_type
 
        predictor = TabularPredictor.load(save_path)

        res_fi = predictor.feature_importance(dat_cv_Tab, "NeuralNetFastAI_BAG_L1")

        res_fi_list.append(res_fi)

res_test_df = pd.concat(res_fi_list, axis=0)

file_name = "./ML/model/Model_fi.csv"

res_test_df.to_csv(file_name, index=False)
