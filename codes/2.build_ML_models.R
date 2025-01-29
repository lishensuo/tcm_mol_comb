library(tidyverse)
library(parallel)

## 特征
fp_types=gsub(".csv","",list.files("ML/data/raw/fp_comb/10114/"))
fp_types=grep("Count",fp_types, value=T,invert=T)
fp_types
# [1] "AtomPairs2D"  "CDK"          "CDKextended"  "CDKgraphonly" "EState"
# [6] "KlekotaRoth"  "MACCS"        "PubChem"      "Substructure"




# BC cell line exp
assays2 = read.csv("ML/data/middle/comb_df_pan.csv")
assays2 = subset(assays2, Cellosaurus.ID=='CVCL_0031') #MCF7


# Regression
assays2 = assays2 %>% 
    dplyr::mutate(ZIP_c = ZIP, Bliss_c = Bliss,
                  Loewe_c = Loewe, HSA_c = HSA)
summary(assays2$ZIP_c)


discard_outliers = function(data){
    Q1 <- quantile(data, 0.25)
    Q3 <- quantile(data, 0.75)
    IQR <- Q3 - Q1
    # 计算离群点的下界和上界
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR 
    which(!(data < lower_bound | data > upper_bound))
}


# metric = "ZIP"
# test_ratio = 0.1
# seed = 42
# data = assays2
func_data_split = function(data, metric='Loewe', test_ratio = 0.1, seed = 42){
    data_sub = data[,c('Cellosaurus.ID','CID.A','CID.B',paste0(metric,'_c'),'zuhe')]
    # data_sub = data_sub[data_sub[,paste0(metric,'_c')]!="Non",]

    data_sub = data_sub[discard_outliers(data_sub[,paste0(metric,'_c')]),]
    # Min-Max Norm
    data_sub[,paste0(metric,'_c')] = scales::rescale(data_sub[,paste0(metric,'_c')], to = c(0, 1))
    
    colnames(data_sub)[4]='label'
    
    # data_sub$label = ifelse(data_sub$label=="Antag",0,1)
    
    set.seed(seed)
    data_sub$set = sample(c('train','test'),nrow(data_sub),
                          replace=T,prob=c(1-test_ratio, test_ratio))
    return(data_sub)
}


# data_sub = func_data_split(assays2, metric='Loewe', test_ratio = 0.1, seed = 42)
# CL_ID = data_sub$Cellosaurus.ID[1]
# CL_Exp = ccle_exp_pan_bc
# CID.A = data_sub$CID.A[1]
# CID.B = data_sub$CID.B[1]
# fp_type="MACCS"
# gene_rank=ncol(ccle_exp_pan_bc)

cmap_gene = data.table::fread("PW/data/raw/CMAP_geneinfo_beta.txt", data.table=F)
cmap_gene = cmap_gene %>%
	dplyr::filter(feature_space=="landmark") %>%
	dplyr::pull(gene_symbol)
ccle_exp = read.csv("ML/data/raw/RNA__RNA_seq_composite_expression.csv") %>% 
	tibble::column_to_rownames("X")
ccle_exp = ccle_exp[rownames(ccle_exp) %in% cmap_gene,]


get_feat_data = function(data_sub, fp_sle = NULL, CL_ID = NULL){
    rdkit_desc = read.csv("ML/data/test_bc/RDKit_Descript_101_drugs.csv") %>% 
        dplyr::select(!mol_smile)

    feat_desc_1 = rdkit_desc[match(data_sub$CID.A, rdkit_desc$CID), -1]
    colnames(feat_desc_1) = paste0("Desc_",colnames(feat_desc_1),"_1")
    rownames(feat_desc_1) = seq(nrow(feat_desc_1))

    feat_desc_2 = rdkit_desc[match(data_sub$CID.B, rdkit_desc$CID), -1]
    colnames(feat_desc_2) = paste0("Desc_",colnames(feat_desc_2),"_2")
    rownames(feat_desc_2) = seq(nrow(feat_desc_2))

    # test = lapply(rdkit_desc$CID, function(x){
    #     fp_tmp2 = lapply(fp_types, function(fp_type){
    #         fp_tmp = read.csv(paste0("ML/data/raw/fp_comb/",x,"/",fp_type,".csv"),
    #                             row.names=1)
    #         colnames(fp_tmp) = paste0(fp_type,"_",colnames(fp_tmp))
    #         fp_tmp
    #     }) %>% do.call(cbind, .)
    #     fp_tmp2
    # }) %>% do.call(rbind, .) %>% 
    #     tibble::rownames_to_column("CID")
    # write.csv(test, file = "ML/data/test_bc/FP_9Merge_101_drugs.csv", row.names=FALSE)

    fp_merge = read.csv("ML/data/test_bc/FP_9Merge_101_drugs.csv")
	if(!is.null(fp_sle)){
		fp_merge = fp_merge %>% 
			dplyr::select("CID",starts_with(fp_sle))
	}
    fp_merge_1 = fp_merge[match(data_sub$CID.A, fp_merge$CID), -1]
    colnames(fp_merge_1) = paste0(colnames(fp_merge_1),"_1")
    rownames(fp_merge_1) = seq(nrow(fp_merge_1))

    fp_merge_2 = fp_merge[match(data_sub$CID.B, fp_merge$CID), -1]
    colnames(fp_merge_2) = paste0(colnames(fp_merge_2),"_2")
    rownames(fp_merge_2) = seq(nrow(fp_merge_2))

    dat_feat_desc = rbind(as.matrix(cbind(feat_desc_1, feat_desc_2)),
                          as.matrix(cbind(feat_desc_2, feat_desc_1)))
    dat_fp_merge = rbind(as.matrix(cbind(fp_merge_1, fp_merge_2)),
                         as.matrix(cbind(fp_merge_2, fp_merge_1)))                      
    dat_feat = as.data.frame(cbind(dat_feat_desc, dat_fp_merge))

	if(!is.null(CL_ID)){
		exp_dat = as.data.frame(t(ccle_exp[, rep(CL_ID,nrow(dat_feat))]))
		colnames(exp_dat) = paste0("Exp_", colnames(exp_dat))
		dat_feat = cbind(dat_feat, exp_dat)
	}
	return(dat_feat)
}




# set.seed(123)
# CL_IDs = sample(setdiff(colnames(ccle_exp),'CVCL_0031'),9)
# CL_IDs = c(CL_IDs, "CVCL_0031")


CL_IDs = "CVCL_0031"
CL_ID = "CVCL_0031"




for (metrics_sle in c("ZIP","HSA","Loewe","Bliss")){

    print(paste0("### ", metrics_sle))

	data_final = lapply(CL_IDs, function(CL_ID){
		# CL_ID = CL_IDs[1]
		print(CL_ID)
		assays2_sub = subset(assays2, Cellosaurus.ID==CL_ID) #MCF7

		data_sub = func_data_split(assays2_sub, metric=metrics_sle, test_ratio = 0.2, seed = 123)

		dat_feat = get_feat_data(data_sub, fp_sle = "PubChem", CL_ID = CL_ID)

		data_final = rbind(data_sub,data_sub)[,c("label","set")] %>% tibble::remove_rownames()
		data_final = cbind(data_final, dat_feat)
		data_final
	}) %>% do.call(rbind, .)
	
	print(dim(data_final))

    write.csv(data_final, row.names = FALSE,
                file = paste0("./ML/data/test_bc/train_dat_",metrics_sle,"_Pan_MCF7.csv"))
}


data_final_2 = read.csv(paste0("./ML/data/test_bc/train_dat_",metrics_sle,"_Pan_MCF7.csv"))
dim(data_final_2)


data_final[1:4, 1:4]

data_final_2[1:4, 1:4]




#### TCM data
library(tidyverse)
library(parallel)

tcm_df = read.csv("ML/data/raw/tcm/tcm_df_all.csv")
tcm_df = tcm_df %>% 
    dplyr::select(!DepMap_ID) %>% 
    dplyr::distinct()
    
tcm_df$zuhe = apply(tcm_df[,c('tcmA','tcmB')],1,function(x){ paste(sort(x),collapse = "_")})


tcm_df_sle <- data.frame(
  tcmA = c("S14S25", "S14S25", "S14S25", "S14S25", "S14S25", 
           "S14S25", "S14S25", "S14S25", "S14S25", "S10S12", "S10S12"),
  tcmB = c("S4S21", "S4S14", "S4S15", "S2S15", "S5S27", 
           "S4S18", "S4S2", "S4S20", "S5S24", "S2S15", "S2S3")
)
tcm_df_sle$zuhe = apply(tcm_df_sle[,c('tcmA','tcmB')],1,function(x){ paste(sort(x),collapse = "_")})



tcm_df2 = tcm_df %>%
    dplyr::filter(zuhe %in% tcm_df_sle$zuhe) %>% 
    dplyr::rename('CID.A'='CID_A','CID.B'='CID_B')



rdkit_desc_tcm = read.csv("ML/data/test_bc/RDKit_Descript_12_tcms.csv") %>% 
        dplyr::select(!mol_smile)

feat_desc_tcm_1 = rdkit_desc_tcm[match(tcm_df2$CID.A, rdkit_desc_tcm$CID), -1]
colnames(feat_desc_tcm_1) = paste0("Desc_",colnames(feat_desc_tcm_1),"_1")
rownames(feat_desc_tcm_1) = seq(nrow(feat_desc_tcm_1))

feat_desc_tcm_2 = rdkit_desc_tcm[match(tcm_df2$CID.B, rdkit_desc_tcm$CID), -1]
colnames(feat_desc_tcm_2) = paste0("Desc_",colnames(feat_desc_tcm_2),"_2")
rownames(feat_desc_tcm_2) = seq(nrow(feat_desc_tcm_2))



fp_merge_tcm = read.csv("ML/data/test_bc/FP_9Merge_12_tcms.csv")

fp_merge_tcm_1 = fp_merge_tcm[match(tcm_df2$CID.A, fp_merge_tcm$CID), -1]
colnames(fp_merge_tcm_1) = paste0(colnames(fp_merge_tcm_1),"_1")
rownames(fp_merge_tcm_1) = seq(nrow(fp_merge_tcm_1))

fp_merge_tcm_2 = fp_merge_tcm[match(tcm_df2$CID.B, fp_merge_tcm$CID), -1]
colnames(fp_merge_tcm_2) = paste0(colnames(fp_merge_tcm_2),"_2")
rownames(fp_merge_tcm_2) = seq(nrow(fp_merge_tcm_2))

dat_feat_desc_tcm = rbind(as.matrix(cbind(feat_desc_tcm_1, feat_desc_tcm_2)),
                          as.matrix(cbind(feat_desc_tcm_2, feat_desc_tcm_1)))

dat_fp_merge_tcm = rbind(as.matrix(cbind(fp_merge_tcm_1, fp_merge_tcm_2)),
                         as.matrix(cbind(fp_merge_tcm_2, fp_merge_tcm_1)))                      
dat_feat_tcm = as.data.frame(cbind(dat_feat_desc_tcm, dat_fp_merge_tcm))
dat_feat_tcm = dat_feat_tcm %>% 
    dplyr::mutate(Id = c(tcm_df2$zuhe,tcm_df2$zuhe), .before = 1)
write.csv(dat_feat_tcm, row.names = FALSE,
            file = paste0("./ML/data/test_bc/tcm_dat.csv"))






library(tidyverse)

test_res = read.csv("./ML/data/test_bc/model05/Model_test_result.csv")


cv_res_best = test_res %>%
	dplyr::filter(FP=='PubChem') %>% 
    dplyr::filter(model!='WeightedEnsemble_L2') %>%
    # dplyr::filter(!model %in% c('NeuralNetFastAI','NeuralNetTorch')) %>%
	dplyr::slice_max(order_by=score_val, n=1, by=Metrics) %>% 
    dplyr::distinct(Metrics, .keep_all = T) 
	