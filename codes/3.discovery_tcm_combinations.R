

## TCM calculation
pw_meta = pathway_state
pw_final = pw_meta[,"Term"]
pw_ID = pw_meta[,"ID"]

tcm_fd = read.csv("./PW/data/raw/TCM500_deg_P0.01.csv", row.names=1)
tcm_fd = na.omit(tcm_fd)
tcm_fd_FC1 = tcm_fd
tcm_fd_FC1[abs(tcm_fd_FC1)<0.5] = 0
tcm_fd_FC1_sign = sign(tcm_fd_FC1)
tcm_fd_FC1_sign_neg = -1*tcm_fd_FC1_sign
tcm_fd_FC1_sign_neg$gene = rownames(tcm_fd_FC1_sign_neg)

AC_tcm.list = lapply(seq(pw_final), function(k){
	print(paste0("Pathway-",k))
	gset_sub = subset(gset, term %in% pw_final[k])
	genes = gset_sub$gene
	ppi_sub = subset(ppi, SYMBOL_1 %in% genes & SYMBOL_2 %in% genes)
	dim(ppi_sub)
	g1 = graph_from_data_frame(d = ppi_sub, directed = FALSE)
	g1_BW = betweenness(g1) %>% as.data.frame() %>% 
		tibble::rownames_to_column("gene") %>% 
		dplyr::rename("BW"=".") %>% 
		dplyr::mutate(BW_log = log2(BW+1)) %>% 
		dplyr::arrange(desc(BW_log))
	bp_genes_sign = gset %>% 
		dplyr::filter(term %in% pw_final[k],
					  gene %in% rownames(exp_cmap),
					  gene %in% rownames(tcm_fd)) %>% 
		dplyr::inner_join(g1_BW[,c("gene","BW_log")]) %>%
		dplyr::inner_join(subset(deg, padj < 0.01 & abs(log2FoldChange)>0.5)[,c("X","log2FoldChange")],
						 by=c("gene"="X")) %>% 
		dplyr::rename(BRCA_log2FC=log2FoldChange) %>% 
		dplyr::mutate(BRCA_sign = sign(BRCA_log2FC)) %>% 
		dplyr::left_join(tcm_fd_FC1_sign_neg) %>% 
		dplyr::mutate(across(!c(term,gene,BW_log,BRCA_log2FC), ~ factor(.x, levels=c(-1,0,1)))) 
	AC_res = mclapply(seq(ncol(tcm_fd_FC1)), function(j){
				cp_sle = colnames(tcm_fd_FC1)[j]
				ct_raw = bp_genes_sign[,c("BW_log", "BRCA_sign", cp_sle)]
				colnames(ct_raw)[3] = "CP_sign"
				ct_stat = ct_raw %>% 
					dplyr::group_by(BRCA_sign, CP_sign) %>% 
					dplyr::summarise(sums = sum(BW_log)) %>% 
					dplyr::mutate(BRCA_sign2 = case_when(BRCA_sign == "-1" ~ "d_neg",
														BRCA_sign == "1" ~ "d_posi",
														BRCA_sign == "0" ~ "d_zero")) %>% 
					dplyr::mutate(CP_sign2 = case_when(CP_sign == "-1" ~ "cp_neg",
													  CP_sign == "1" ~ "cp_posi",
													  CP_sign == "0" ~ "cp_zero"))
				ct = matrix(0, nrow=3, ncol=3,
							dimnames=list(c("d_neg","d_zero","d_posi"),
										  c("cp_neg", "cp_zero", "cp_posi")))
				for(i in seq(nrow(ct_stat))){
					ct[ct_stat$BRCA_sign2[i], ct_stat$CP_sign2[i]] = ct_stat$sums[i]
				}
				ac2 = gwet.ac1.table(ct, weights)$coeff.val
				ac2
				return(ac2)
				}, mc.cores=15) %>% unlist
	return(AC_res)
})

AC_tcm = do.call(cbind, AC_tcm.list)
colnames(AC_tcm) = pw_final
rownames(AC_tcm) = colnames(tcm_fd_FC1)


AC_ttd_brca = AC_ttd[meta_cmap$sig_id[meta_cmap$Anno=="BRCA"],pw_final]
cut_1 = apply(AC_ttd_brca, 2, function(x){quantile(x, prob=0.5)})
cut_2 = apply(AC_ttd_brca, 2, function(x){quantile(x, prob=0.8)})


AC_tcm_reshaped = AC_tcm %>%
	reshape2::melt() %>% 
	dplyr::rename(TCM=Var1, PW=Var2)

AC_tcm_reshaped = AC_tcm_reshaped %>%
	dplyr::left_join(as.data.frame(cut_1) %>% 
					 tibble::rownames_to_column("PW")) %>% 
	dplyr::left_join(as.data.frame(cut_2) %>% 
					 tibble::rownames_to_column("PW"))

AC_tcm_reshaped = AC_tcm_reshaped %>%
	dplyr::mutate(score = case_when(
		value>cut_2 ~ 1,
		value>cut_1 ~ 0.5,
		TRUE ~ 0
		)
	)



rank_single = AC_tcm_reshaped %>%
	dplyr::group_by(TCM) %>%
	dplyr::summarise(Score=sum(score), PW = sum(score>0)) %>% 
	dplyr::rename(score=Score) %>% 
	dplyr::arrange(desc(score), desc(PW)) %>%
	as.data.frame() %>%
	dplyr::mutate(TCM=as.character(TCM))


candi_tcm = rank_single$TCM[rank_single$score>0]
candi_zuhe = combn(candi_tcm,2) %>% 
	t() %>% as.data.frame() %>% 
	dplyr::rename(tcmA=V1, tcmB=V2)
candi_zuhe$score = mclapply(seq(nrow(candi_zuhe)), function(i) {
	print(i)
	score_merge = dplyr::inner_join(
		AC_tcm_reshaped %>%
			dplyr::filter(TCM==candi_zuhe$tcmA[i]) %>%
			dplyr::select(PW, score) %>%
			dplyr::rename(score_1=score),
		AC_tcm_reshaped %>%
			dplyr::filter(TCM==candi_zuhe$tcmB[i]) %>%
			dplyr::select(PW, score) %>%
			dplyr::rename(score_2=score)
	) %>% dplyr::select(!PW) %>%
		dplyr::mutate(merge = score_1+score_2) %>%
		dplyr::mutate(merge2 = ifelse(merge>1, 1, merge))
	sum(score_merge$merge2)
}, mc.cores=40) %>% unlist()

candi_zuhe = candi_zuhe %>%
	dplyr::arrange(desc(score)) %>% 
	dplyr::mutate(testA = ifelse(tcmA>tcmB,tcmB,tcmA)) %>% 
	dplyr::mutate(testB = ifelse(tcmA<tcmB,tcmB,tcmA)) %>% 
	dplyr::mutate(tcmAB = paste0(testA,"_",testB)) %>% 
	dplyr::select(-testA, -testB)



## machine learning prediction

# utils/script_fingerprint_padelpy.py
tcm_meta = readxl::read_xlsx("tcm_data/metadata.xlsx")
tcm_meta_sig = tcm_meta %>% 
	dplyr::filter(id %in% candi_tcm)
tcm_zuhe = t(combn(tcm_meta_sig$id,2)) %>% as.data.frame()
colnames(tcm_zuhe) = c("tcmA","tcmB")
tcm_zuhe = tcm_zuhe %>%
	dplyr::left_join(tcm_meta_sig[,2:3],by=c("tcmA"="id")) %>%
	dplyr::rename(CID_A=pubchem_CID) %>%
	dplyr::left_join(tcm_meta_sig[,2:3],by=c("tcmB"="id")) %>%
	dplyr::rename(CID_B=pubchem_CID) %>% 
	dplyr::mutate(DepMap_ID="ACH-000019")

for (i in seq(nrow(model_stat_optimal))){
	metrics = model_stat_optimal$metrics[i]
	fp_type = model_stat_optimal$fp_type[i]
	gene_rank = model_stat_optimal$gene_rank[i]
	DepMap_ID = "ACH-000019"

	feat_merge_df = mclapply(seq(nrow(tcm_zuhe)), function(j){
		print(j)
		CID_A = tcm_zuhe$CID_A[j]
		CID_B = tcm_zuhe$CID_B[j]
		merge_feat_func(DepMap_ID, CID_A, CID_B, fp_type, gene_rank)
	}, mc.cores=20) %>% do.call(rbind, .)
	dim(feat_merge_df)
	write.csv(feat_merge_df, row.names=F,
		 file=paste0("tcm_feature/",metrics,".csv"))
}

write.csv(tcm_zuhe, file = "tcm_data/tcm_zuhe.csv")


# utils/script_automl_predict.py

tcm_predict = data.table::fread("tcm_predict.csv", data.table=FALSE)


prob_mean = tcm_predict %>%
	dplyr::mutate(testA = ifelse(tcmA>tcmB,tcmB,tcmA)) %>% 
	dplyr::mutate(testB = ifelse(tcmA<tcmB,tcmB,tcmA)) %>% 
	dplyr::mutate(tcmAB = paste0(testA,"_",testB)) %>% 
	dplyr::group_by(tcmAB, metrics) %>%
	dplyr::summarise(prob_mean = mean(prob)) %>%
	tidyr::pivot_wider(names_from =metrics, values_from =prob_mean) %>% as.data.frame()

ml_predict$Overall = apply(ml_predict[,-1],1,mean)

ml_predict_sub = ml_predict[match(candi_zuhe$tcmAB[1:11],ml_predict$tcmAB),] %>%
	dplyr::arrange(desc(Overall))
#           tcmAB     Bliss       HSA     Loewe       ZIP   Overall
# 1   S10S12_S2S3 0.8262045 0.8082842 0.4470000 0.9921178 0.7684016
# 2  S14S25_S4S15 0.7447488 0.7480226 0.4441389 0.9121419 0.7122631
# 3  S10S12_S2S15 0.7179576 0.3941357 0.4436111 0.9810124 0.6341792
# 4  S14S25_S4S18 0.8621383 0.6955093 0.4187778 0.5354554 0.6279702
# 5  S14S25_S4S21 0.3337803 0.6430628 0.4453611 0.8985180 0.5801806
# 6  S14S25_S5S27 0.3680591 0.5835467 0.3543040 0.8965764 0.5506215
# 7  S14S25_S5S24 0.3680591 0.5835467 0.3543040 0.8965764 0.5506215
# 8  S14S25_S4S20 0.3237390 0.3908321 0.4386111 0.9962140 0.5373490
# 9   S14S25_S4S2 0.3365237 0.5615735 0.2714206 0.8671755 0.5091733
# 10 S14S25_S2S15 0.3564575 0.3070088 0.3360556 0.9846147 0.4960341
# 11 S14S25_S4S14 0.2241391 0.4080410 0.2919206 0.8484404 0.4431353