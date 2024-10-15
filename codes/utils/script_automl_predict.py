import os
import time
import sys
import pandas as pd
import numpy as np
from autogluon.tabular import TabularDataset, TabularPredictor
from functools import reduce

mode = "medium_quality"

# cv_best = pd.read_csv("ML/data/test_bc/model04/Model_cv_best.csv")
# cv_best['FP'] = 'PubChem'
# cv_best['model'] = 'XGBoost'


metricses = ['HSA','ZIP','Loewe','Bliss']
fp_types = ["MACCS","PubChem","Substructure"]
models = ["WeightedEnsemble_L2","NeuralNetFastAI","NeuralNetTorch","LightGBM",
          "XGBoost","LightGBMLarge","RandomForestEntr",
          "LightGBMXT","CatBoost","RandomForestGini","ExtraTreesGini","ExtraTreesEntr",
        #   "KNeighborsDist","KNeighborsUnif"
          ]

# models = ["WeightedEnsemble_L2","NeuralNetFastAI","NeuralNetTorch","LightGBM","XGBoost","LightGBMLarge","RandomForestMSE",
#           "LightGBMXT","CatBoost","ExtraTreesMSE"
#         #   ,"KNeighborsDist","KNeighborsUnif"
#           ]


tcm_feat = pd.read_csv("ML/data/test_bc/tcm_dat.csv")
# data_Tab = TabularDataset(tcm_feat)

# i = 0
# metrics = metricses[i]
# fp_type = fp_types[i]
# model = models[i]

def model_predict(metrics, fp_type, model):
    
    selected_columns = [col for col in tcm_feat.columns if col.startswith(('Desc',fp_type))]
    # selected_columns = [col for col in tcm_feat.columns if col.startswith(fp_type)]

    data_Tab = TabularDataset(tcm_feat[selected_columns])

    save_path = "./ML/data/test_bc/model07/" + metrics + "_"+ mode + "/" + fp_type
    predictor = TabularPredictor.load(save_path)
    preds_df = predictor.predict_proba(data_Tab,model=model).iloc[:,[1]]
    # preds_df = pd.DataFrame(predictor.predict(data_Tab,model=model))
    preds_df.columns = ['Prob']
    preds_df['metrics'] = metrics
    preds_df['fp_type'] = fp_type
    preds_df['model'] = model

    preds_df = pd.concat([tcm_feat['Id'], preds_df],axis=1)
    
    return preds_df


# for metrics in metricses:
#     for fp_type in fp_types:

import itertools

metrics_df_merge = []

for metrics, fp_type, model in list(itertools.product(metricses, fp_types, models)):
    print(metrics, fp_type, model)
    metrics_df = model_predict(metrics, fp_type, model)
    metrics_df_merge.append(metrics_df)

merged_df = pd.concat(metrics_df_merge, axis=0)




merged_df.to_csv("./ML/data/test_bc/model07/tcm_pred.csv",index=False)


