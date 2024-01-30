import os
import time
import sys
import pandas as pd
import numpy as np
from autogluon.tabular import TabularDataset, TabularPredictor


mode = "medium_quality"
metricses = ['ZIP','HSA','Loewe','Bliss']
fp_types = ["AtomPairs2D","CDK","CDKextended","CDKgraphonly","EState","KlekotaRoth","MACCS","PubChem","Substructure"]
gene_ranks = [100,300,500]


for metrics in metricses:
	for fp_type in fp_types:
		for gene_rank in gene_ranks:
			data_kfold = pd.read_csv("data_fold/"+metrics+".csv")
			fl = "data_fold/"+metrics+"_"+fp_type+"_"+str(gene_rank)+"_raw.csv"
			zuhe_raw = pd.read_csv(fl)
			zuhe_raw = zuhe_raw.drop(["CID.A","CID.B","DepMap_ID"],axis=1)
			test_res_kfold = []
			for k in range(10):
				print(k)
				# k = 0
				sle_kfold = pd.DataFrame(data_kfold.iloc[:,k])
				sle_kfold.columns=["Set"]
				zuhe_dat = pd.concat([sle_kfold,zuhe_raw],axis=1)

				zuhe_train = zuhe_dat.loc[zuhe_dat.Set==0,:].drop(["Set"], axis=1).reset_index(drop=True)
				zuhe_test = zuhe_dat.loc[zuhe_dat.Set!=0,:].drop(["Set"], axis=1).reset_index(drop=True)
				zuhe_test.shape

				zuhe_train_Tab = TabularDataset(zuhe_train)
				zuhe_test_Tab = TabularDataset(zuhe_test)

				save_path = "model_fold/" + metrics + "_"+ mode + "/" + fp_type + "_" + str(gene_rank) + "_f" + str(k)

				metric = "accuracy"
				eval_metric = ["roc_auc","average_precision","precision","recall","f1"]
				predictor = TabularPredictor(label="label", path=save_path, 
					verbosity=1).fit(zuhe_train_Tab,presets=mode, excluded_model_types = [""])

				test_res = predictor.leaderboard(zuhe_test_Tab, extra_metrics = eval_metric,silent=True)
				test_res = test_res.iloc[:,0:8]
				test_res["fold"] = k
				test_res_kfold.append(test_res)
				# time.sleep(5)

			test_res_df = pd.concat(test_res_kfold, axis=0)
			file_name = "model_fold/eval/"+metrics+"_" + fp_type + "_" + str(gene_rank) + "_kfold.csv"
			test_res_df.to_csv(file_name, index=False)


