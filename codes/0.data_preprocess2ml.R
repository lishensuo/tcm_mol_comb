################# For machine learning #########################
library(tidyverse)
library(parallel)

## Chemical fingerprint
fp_types = c("MACCS", "PubChem", "Substructure")

# SynergyDB data for breast cancer
synergxdb_BC = read.csv("synergydb/comb_df_pan_bc.csv")
# SynergyDB data for MCF7
synergxdb_mcf7 = subset(synergxdb_BC, Cellosaurus.ID=='CVCL_0031') 


# Copy the raw synergy score
synergxdb_mcf7 = synergxdb_mcf7 %>% 
    dplyr::mutate(ZIP_c = ZIP, 
                  Bliss_c = Bliss,
                  Loewe_c = Loewe, 
                  HSA_c = HSA)
summary(synergxdb_mcf7$ZIP_c)
summary(synergxdb_mcf7$Bliss_c)
summary(synergxdb_mcf7$Loewe_c)
summary(synergxdb_mcf7$HSA_c)


# Discard the outlier according to ± 1.5*IQR
discard_outliers = function(data){
    Q1 <- quantile(data, 0.25)
    Q3 <- quantile(data, 0.75)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR 
    which(!(data < lower_bound | data > upper_bound))
}


# metric = "Loewe"
# test_ratio = 0.1
# seed = 42
# data = synergxdb_mcf7
func_data_split = function(data, metric='Loewe', test_ratio = 0.1, seed = 42){
    data_sub = data[,c('Cellosaurus.ID','CID.A','CID.B',paste0(metric,'_c'),'zuhe')]

    # Discard outliers
    data_sub = data_sub[discard_outliers(data_sub[,paste0(metric,'_c')]),]
    # Min-Max Norm
    data_sub[,paste0(metric,'_c')] = scales::rescale(data_sub[,paste0(metric,'_c')], to = c(0, 1))
    
    colnames(data_sub)[4]='label'

    set.seed(seed)
    data_sub$set = sample(c('train','test'),nrow(data_sub),
                          replace=T,prob=c(1-test_ratio, test_ratio))
    return(data_sub)
}


# data_sub = func_data_split(synergxdb_mcf7, metric='Loewe', test_ratio = 0.1, seed = 42)
# CID.A = data_sub$CID.A[1]
# CID.B = data_sub$CID.B[1]
# fp_type="MACCS"



get_feat_data = function(data_sub){

    fp_merge = read.csv("synergydb/fingerprint_compounds.csv")

    fp_merge_1 = fp_merge[match(data_sub$CID.A, fp_merge$CID), -1]
    colnames(fp_merge_1) = paste0(colnames(fp_merge_1),"_1")
    rownames(fp_merge_1) = seq(nrow(fp_merge_1))

    fp_merge_2 = fp_merge[match(data_sub$CID.B, fp_merge$CID), -1]
    colnames(fp_merge_2) = paste0(colnames(fp_merge_2),"_2")
    rownames(fp_merge_2) = seq(nrow(fp_merge_2))

    dat_fp_merge = rbind(as.matrix(cbind(fp_merge_1, fp_merge_2)),
                         as.matrix(cbind(fp_merge_2, fp_merge_1)))                      
    dat_feat = dat_fp_merge
    as.data.frame(dat_feat)
}





for (metrics_sle in c("ZIP","HSA","Loewe","Bliss")){
    # metrics_sle = "ZIP"
    print(paste0("### ", metrics_sle))
    # Split data
    data_sub = func_data_split(synergxdb_mcf7, metric=metrics_sle, test_ratio = 0.2, seed = 123)
    # data feature
    dat_feat = get_feat_data(data_sub)
    dat_feat = dat_feat %>% dplyr::select(starts_with(c("MACCS", "PubChem", "Substructure")))
    # Merge
    data_final = rbind(data_sub,data_sub)[,c("label","set")] %>% tibble::remove_rownames()
    data_final = cbind(data_final, dat_feat)
    print(dim(data_final))
    write.csv(data_final, row.names = FALSE,
                file = paste0("./ML/data/train_dat_",metrics_sle,".csv"))
}






####### Prepare the TCM combination feature
library(tidyverse)
library(parallel)

tcm_df = read.csv("ML/tcm/tcm_df_all.csv")
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
#     tcmA  tcmB         zuhe
# 1 S14S25 S4S21 S14S25_S4S21
# 2 S14S25 S4S14 S14S25_S4S14


tcm_df2 = tcm_df %>%
    dplyr::filter(zuhe %in% tcm_df_sle$zuhe) %>% 
    dplyr::rename('CID.A'='CID_A','CID.B'='CID_B')
#    tcmA   tcmB   CID.A    CID.B         zuhe
# 1  S2S3 S10S12 5280633    72303  S10S12_S2S3
# 2 S2S15 S10S12  185617    72303 S10S12_S2S15



# ##### TCM descriptor feature
# rdkit_desc_tcm = read.csv("ML/data/test_bc/RDKit_Descript_12_tcms.csv") %>% 
#         dplyr::select(!mol_smile)

# feat_desc_tcm_1 = rdkit_desc_tcm[match(tcm_df2$CID.A, rdkit_desc_tcm$CID), -1]
# colnames(feat_desc_tcm_1) = paste0("Desc_",colnames(feat_desc_tcm_1),"_1")
# rownames(feat_desc_tcm_1) = seq(nrow(feat_desc_tcm_1))

# feat_desc_tcm_2 = rdkit_desc_tcm[match(tcm_df2$CID.B, rdkit_desc_tcm$CID), -1]
# colnames(feat_desc_tcm_2) = paste0("Desc_",colnames(feat_desc_tcm_2),"_2")
# rownames(feat_desc_tcm_2) = seq(nrow(feat_desc_tcm_2))


# #### TCM fingerprint feature
# fp_merge_tcm = read.csv("ML/tcm/fingerprint_tcms.csv")

# fp_merge_tcm_1 = fp_merge_tcm[match(tcm_df2$CID.A, fp_merge_tcm$CID), -1]
# colnames(fp_merge_tcm_1) = paste0(colnames(fp_merge_tcm_1),"_1")
# rownames(fp_merge_tcm_1) = seq(nrow(fp_merge_tcm_1))

# fp_merge_tcm_2 = fp_merge_tcm[match(tcm_df2$CID.B, fp_merge_tcm$CID), -1]
# colnames(fp_merge_tcm_2) = paste0(colnames(fp_merge_tcm_2),"_2")
# rownames(fp_merge_tcm_2) = seq(nrow(fp_merge_tcm_2))




## Merge above two types of feature
# dat_feat_desc_tcm = rbind(as.matrix(cbind(feat_desc_tcm_1, feat_desc_tcm_2)),
#                           as.matrix(cbind(feat_desc_tcm_2, feat_desc_tcm_1)))

# dat_fp_merge_tcm = rbind(as.matrix(cbind(fp_merge_tcm_1, fp_merge_tcm_2)),
#                          as.matrix(cbind(fp_merge_tcm_2, fp_merge_tcm_1)))                      
dat_feat_tcm = dat_fp_merge_tcm
dat_feat_tcm = dat_feat_tcm %>% 
    dplyr::mutate(Id = c(tcm_df2$zuhe,tcm_df2$zuhe), .before = 1)

write.csv(dat_feat_tcm, row.names = FALSE,
            file = paste0("./ML/data/tcm_dat.csv"))

