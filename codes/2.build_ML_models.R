library(tidyverse)
library(parallel)


## experiment label
assays = read.csv("synergxdb_combos_breast.csv")
assays$ZIP = as.numeric(assays$ZIP)
assays$Bliss = as.numeric(assays$Bliss)
assays$Loewe = as.numeric(assays$Loewe)
assays$HSA = as.numeric(assays$HSA)
assays$Comboscore = as.numeric(assays$Comboscore)
assays = na.omit(assays)
assays = assays[,c("Cell.Line","Compound.A","Compound.B","ZIP", "Bliss","Loewe","HSA","Comboscore")]

assays2 = assays
tmp = t(apply(assays[,c("Compound.A","Compound.B")],1,sort))
assays2$Compound.A = tmp[,1]
assays2$Compound.B = tmp[,2]
assays2$Cell.Line[assays2$Cell.Line=="MCF-7"]="MCF7"

ccle_meta = read.csv("CCLE_meta.csv")
ccle_meta_br = ccle_meta %>%
	dplyr::filter(primary_disease=="Breast Cancer")
drug_info = read.csv("synergydb_compound.csv")
assays2 = assays2 %>%
	dplyr::left_join(drug_info[,c(1,3)],by=c("Compound.A"="Name")) %>% 
	dplyr::rename(CID.A=PubChem.CID) %>%
	dplyr::left_join(drug_info[,c(1,3)],by=c("Compound.B"="Name")) %>% 
	dplyr::rename(CID.B=PubChem.CID) %>% 
	dplyr::mutate(zuhe=paste0(CID.A,"_",CID.B))  %>%
	dplyr::mutate(DepMap_ID = ccle_meta_br$DepMap_ID[match(Cell.Line, ccle_meta_br$cell_line_name)]) %>%
	dplyr::select(DepMap_ID,Compound.A,Compound.B,CID.A,CID.B,zuhe,
				  ZIP, Bliss, Loewe, HSA, Comboscore)

assays2 = assays2 %>%
	dplyr::mutate(ZIP_c = ifelse(ZIP >= quantile(ZIP, 0.9), "Syner",
		ifelse(ZIP <= quantile(ZIP, 0.1), "Antag","Non"))) %>%
	dplyr::mutate(Bliss_c = ifelse(Bliss >= quantile(Bliss, 0.9), "Syner",
		ifelse(Bliss <= quantile(Bliss, 0.1), "Antag","Non"))) %>%
	dplyr::mutate(Loewe_c = ifelse(Loewe >= quantile(Loewe, 0.9), "Syner",
		ifelse(Loewe <= quantile(Loewe, 0.1), "Antag","Non"))) %>%
	dplyr::mutate(HSA_c = ifelse(HSA >= quantile(HSA, 0.9), "Syner",
		ifelse(HSA <= quantile(HSA, 0.1), "Antag","Non"))) %>%
	dplyr::select(DepMap_ID,Compound.A,Compound.B,CID.A,CID.B,zuhe,
				  ZIP_c, Bliss_c, Loewe_c, HSA_c)




## feature
# utils/script_fingerprint_padelpy.py
fp_types = c("AtomPairs2D", "CDK", "CDKextended", "CDKgraphonly", "EState",
			 "KlekotaRoth", "MACCS", "PubChem", "Substructure")



ccle_meta = read.csv("CCLE_meta.csv")
ccle_meta_br = ccle_meta %>%
	dplyr::filter(primary_disease=="Breast Cancer")
ccle_exp = data.table::fread("CCLE_expression.csv") %>%
	as.data.frame() %>%
	tibble::column_to_rownames("V1")
colnames(ccle_exp) = str_split(colnames(ccle_exp)," ", simplify=T)[,1]
ccle_exp_br = ccle_exp[rownames(ccle_exp) %in% ccle_meta_br$DepMap_ID,]
sd_stat = apply(ccle_exp_br, 2, sd)
ccle_exp_br = ccle_exp_br[colnames(ccle_exp_br)[order(sd_stat, decreasing=T)]]
write.csv(ccle_exp_br, file="feature_ccle/CCLE_expression_breast.csv")



## 10 fold--based on drug
merge_feat_func = function(DepMap_ID, CID.A, CID.B, fp_type="MACCS", gene_rank=100){
	fp_tmpA = read.csv(paste0("feature_compound/",CID.A,"/",fp_type,".csv"),
					  row.names=1) %>% t()
	colnames(fp_tmpA) = "CID.A"
	fp_tmpB = read.csv(paste0("feature_compound/",CID.B,"/",fp_type,".csv"),
					  row.names=1) %>% t()
	colnames(fp_tmpB) = "CID.B"
	fp_tmp = cbind(fp_tmpA, fp_tmpB) %>%
		as.data.frame() %>%
		dplyr::mutate(merge=CID.A+CID.B) %>% as.matrix()
	ccle_tmp = t(ccle_exp_br[DepMap_ID,1:gene_rank])
	feat_merge = c(fp_tmp[,3],ccle_tmp[,1])
	return(feat_merge)
}


for (metrics_sle in c("ZIP","HSA","Loewe","Bliss")){
	assays2_sub = assays2[,c('DepMap_ID','CID.A','CID.B',paste0(metrics_sle,'_c'),'zuhe')]
	assays2_sub = assays2_sub[assays2_sub[,paste0(metrics_sle,'_c')]!="Non",]
	colnames(assays2_sub)[4] = "label" 

	data_kfold = lapply(1:10, function(fold){
		# fold = 1
		external_drugs = drug_fold$drug[drug_fold$fold==fold]
		data_split = apply(assays2_sub[,c("CID.A","CID.B")], 1, function(x){
			sum(x %in% external_drugs)
		})
		return(data_split)
	}) %>% do.call(cbind, .)
	colnames(data_kfold) = paste0("f",1:10)
	write.csv(data_kfold, row.names=F,
		file=paste0("data_fold/",metrics_sle,".csv"))

	for (fp_type in fp_types){
		for (gene_rank in c(100,300,500)){
			print(c(fp_type, gene_rank))
			feat_merge_df = mclapply(seq(nrow(assays2_sub)),function(i){
				# i = 1
				# print(i)
				DepMap_ID = assays2_sub$DepMap_ID[i]
				CID.A = assays2_sub$CID.A[i]
				CID.B = assays2_sub$CID.B[i]
				fp_type=fp_type
				gene_rank=gene_rank

				feat_merge = merge_feat_func(DepMap_ID, CID.A, CID.B, fp_type, gene_rank)
				return(feat_merge)
			}, mc.cores=40) %>% do.call(rbind,.) %>% as.data.frame()
			feat_merge_df = cbind(assays2_sub[,c("CID.A","CID.B","label","DepMap_ID"),drop=F],feat_merge_df)
			feat_merge_df$label = ifelse(feat_merge_df$label=="Syner",1,0)
			table(feat_merge_df$label)
			feat_merge_df[1:6,1:6]
			print(dim(feat_merge_df))
			write.csv(feat_merge_df, row.names=F,
				file=paste0("data_fold/",metrics_sle,"_",fp_type,"_",gene_rank,"_raw.csv"))
		}
	}
}




## AutoML modeling
# utils/script_automl_train.py


## model compare
fls = list.files("model_fold/eval", full.names=TRUE)
model_stat = lapply(seq(fls), function(i){
	metrics = str_split(fls[i], "_")[[1]][2]
	fp_type = str_split(fls[i], "_")[[1]][3]
	gene_rank = str_split(fls[i], "_")[[1]][4]
	read.csv(fls[i]) %>%
		dplyr::group_by(model) %>%
		dplyr::summarise(train_acc = mean(score_val),
						 external_acc = mean(score_test),
						 external_auc = mean(roc_auc),
						 external_aupr = mean(average_precision),
						 external_prec = mean(precision),
						 external_recall = mean(recall),
						 external_f1 = mean(f1)) %>%
		dplyr::mutate(metrics=metrics,fp_type=fp_type,gene_rank=gene_rank, ) %>%
		dplyr::mutate(feat = paste0(fp_type,"_",gene_rank))
}) %>% do.call(rbind, .) %>% as.data.frame() 

model_stat_optimal = model_stat %>% 
	dplyr::slice_max(order_by=external_acc, n=1, by=metrics)

