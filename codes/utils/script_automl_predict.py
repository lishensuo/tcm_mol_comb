import os
import sys
import pandas as pd
import numpy as np
from autogluon.tabular import TabularDataset, TabularPredictor

tcm_zuhe = pd.read_csv("ML/data/raw/tcm/tcm_zuhe.csv")
model_stat_optimal = pd.read_csv("model_fold/model_stat_optimal.csv")
mode = 'medium_quality'

def model_predict(data_Tab, metrics, fp_type, gene_rank):
	preds_kfold = []
	for k in range(10):  #kFold预测
		# k = 0
		save_path = "model_predict/" + metrics + "_"+ mode + "/" + fp_type + "_" + str(gene_rank) + "_f" + str(k)
		predictor = TabularPredictor.load(save_path)
		preds_df = predictor.predict_proba(data_Tab,model=model).iloc[:,[1]]
		preds_df.columns = ["prob"]
		preds_df['metrics'] = metrics
		preds_df['kFold'] = k
		preds_df = pd.concat([tcm_zuhe,preds_df],axis=1)
		preds_kfold.append(preds_df)
	preds_kfold_df = pd.concat(preds_kfold, axis=0)
	return preds_kfold_df

metrics_df_merge = []

for j in range(model_stat_optimal.shape[0]):
	print(j)
	metrics = model_stat_optimal.metrics[j]
	model = model_stat_optimal.model[j]
	fp_type = model_stat_optimal.fp_type[j]
	gene_rank = model_stat_optimal.gene_rank[j]

	tcm_feat = pd.read_csv("tcm_feature_"+metrics+".csv")
	data_Tab = TabularDataset(tcm_feat)

	metrics_df = model_predict(data_Tab, metrics, fp_type, gene_rank)

	metrics_df_merge.append(metrics_df)


metrics_df_merge = pd.concat(metrics_df_merge, axis=0)
metrics_df_merge.to_csv("tcm_predict.csv",index=False)


