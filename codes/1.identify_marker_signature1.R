library(clusterProfiler)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(tidyverse)
library(RandomWalkRestartMH)
library(igraph)
library(parallel)
library(irrCAC)
library(cmapR)
library(ggpubr)


## signature collection from MSigDB
gset1 = read.gmt("MsigDB_c2.cp.reactome.v7.5.1.symbols.gmt") 
gset2 = read.gmt("MsigDB_c5.go.bp.v7.5.1.symbols.gmt") 
gset3 = read.gmt("MsigDB_c4.cm.v7.5.1.symbols.gmt")
gset4= read.gmt("MsigDB_c2.cp.kegg.v7.5.1.symbols.gmt") #
gset5 = read.gmt("MsigDB_h.all.v7.5.1.symbols.gmt") #
gset = rbind(gset1, gset2, gset3, gset4, gset5)

gset_meta = gset %>% 
	dplyr::count(term, name="size") %>% 
	dplyr::mutate(ID=paste0("PW",seq(nrow(.)))) %>% 
	dplyr::mutate(type = str_split(term,"_",simplify=T)[,1]) %>%
	dplyr::mutate(tmp = seq(nrow(.))) %>% 
	dplyr::group_by(type) %>%
	dplyr::mutate(tmp2 = rank(tmp)) %>% as.data.frame() %>%
	dplyr::mutate(type_id = paste0(type,"-",tmp2)) %>% 
	dplyr::select(-tmp, -tmp2)


## DEG identification and OAR enrichment of breast cancer
query <- GDCquery(project = "TCGA-BRCA",
                  legacy = FALSE,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification")
GDCdownload(query, files.per.chunk = 10)
data <- GDCprepare(query, save = T, save.filename = "CHOL_RNAseq_hg38.rda")
rownames(data) = rowData(data)$gene_name
data = data[!duplicated(rownames(data)),]
data = data[rowData(data)$gene_type=="protein_coding",]
data = data[,data$sample_type!="Metastatic"]
data = data[,!data$barcode  %in% c("TCGA-A7-A0DC-01B-04R-A22O-07",
								   "TCGA-A7-A0DB-01C-02R-A277-07",
								   "TCGA-A7-A0DB-01A-11R-A277-07",
								   "TCGA-A7-A13E-01A-11R-A277-07",
								   "TCGA-A7-A13E-01B-06R-A277-07")]
count = assay(data, "unstranded")
meta = colData(data)[,c("patient", "shortLetterCode", "barcode", "sample")] 
meta$Group = factor(meta$shortLetterCode, levels = c("NT", "TP"))
dds = DESeqDataSetFromMatrix(countData = count,
                              colData = meta,
                              design= ~ Group)
dds = DESeq(dds)  
res = results(dds, name="Group_TP_vs_NT")
deg = as.data.frame(res)
deg = deg[order(deg$padj), ]
deg_sig = dplyr::filter(deg, padj<0.01 & abs(log2FoldChange)>1) 

enrich_obj = enricher(deg_sig$X, TERM2GENE = gset, pvalueCutoff = 0.01,
	minGSSize = 5, maxGSSize = 500)
enrich_brca_df = enrich_obj@result %>% 
	dplyr::filter(p.adjust<0.01) %>% 
	dplyr::select(ID, p.adjust, Count) %>% 
	tibble::remove_rownames() %>% 
	dplyr::mutate(logP= -log10(p.adjust)) %>%
	dplyr::mutate(Type=str_split(ID,"_",simplify=T)[,1]) %>%
	dplyr::left_join(gset_meta[,c("term","size")], by=c("ID"="term")) 
dim(enrich_brca_df)
pw_filt_1 = enrich_brca_df$ID




## RWR Proximity to BC targets
t2d = read.csv("TTD_target2disease_v8101.csv")
tid =  read.csv("TTD_target_info_v8101.csv",row.names=1)
t2d_brca = t2d %>%
	dplyr::filter(DiseaseName=="Breast cancer") %>%
	dplyr::left_join(tid) %>%
	dplyr::filter(!is.na(GENENAME)) %>%
	dplyr::filter(TARGTYPE=="Successful target")

brca_targets = t2d_brca$GENENAME
brca_targets = str_split(brca_targets, "; ") %>% unlist %>% unique
brca_targets = str_split(brca_targets, "-") %>% unlist %>% unique
brca_targets = setdiff(brca_targets, c("Candi TMP1","HIV RT"))
brca_targets = sort(brca_targets)

ppiDat= data.table::fread("PPI_datasets_symbol.csv")[,c("SYMBOL_1","SYMBOL_2")] %>% na.omit()

g1 = graph_from_data_frame(d = ppiDat, directed = FALSE)
PPI_MultiplexObject <- create.multiplex(list(PPI=g1))
AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)

brca_targets = brca_targets[brca_targets %in% as_ids(V(g1))]
SeedGene <- brca_targets
RWR_PPI_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI,
                                                 PPI_MultiplexObject,SeedGene)
rwr_score=RWR_PPI_Results$RWRM_Results %>% 
			rbind(data.frame(NodeNames=SeedGene, Score=0.7/length(SeedGene)),.)
rwr_score$Scale = scale(rwr_score$Score)


pw2targets_p = mclapply(1:length(pw_filt_1), function(i){
	print(i)
	set.seed(123)
	pw2targets = gset %>% 
			dplyr::filter(term==pw_filt_1[i]) %>% 
			dplyr::inner_join(rwr_score,.,by=c("NodeNames"="gene"))

	perms = lapply(1:1000, function(x){
		sp_idx = sample(1:nrow(rwr_score), nrow(pw2targets))
		temp = sum(rwr_score[sp_idx,"Score"])
		return(temp)
	}) %>% unlist

	temp2 = dnorm((sum(pw2targets[,"Score"]) - mean(perms))/sd(perms))
	return(temp2)
}, mc.cores=10) %>% unlist

names(pw2targets_p) = pw_filt_1
pw2targets_p_adj = p.adjust(pw2targets_p, method="BH")
pw_filt_2 = pw2targets_p_adj[pw2targets_p_adj<0.01] %>% names()


## signatue de-redundancy
jccard_min = matrix(NA, nrow=length(pw_filt_2), ncol=length(pw_filt_2))
colnames(jccard_min)=pw_filt_2
rownames(jccard_min)=pw_filt_2
gset.list = split(gset$gene, gset$term)
for(i in seq(pw_filt_2)){
	print(i)
	for(j in seq(pw_filt_2)){
		id1 = gset.list[[pw_filt_2[i]]]
		id2 = gset.list[[pw_filt_2[j]]]
		jccard_min[i,j] = length(intersect(id1, id2))/min(length(id1),length(id2))
	}
}
dup_terms = lapply(seq(nrow(jccard_min)-1), function(i){
	# i = 1
	print(i)
	test = jccard_min[i, (i+1):nrow(jccard_min)]
	return(names(test[test>=0.5]))
}) %>% unlist %>% unique()

pw_filt_3 = setdiff(colnames(jccard_min), dup_terms)




## wAC index analysis
d2d = read.csv("TTD_drug2disease_v8101.csv", row.names = 1)
d2BRCA = d2d[grep("breast", d2d$DiseaseName, ignore.case = T),]
dim(d2BRCA) # [1] 419   4
druginfo = read.csv("TTD_drug_info_v8101.csv",row.names = 1)
drugBRCA_info = druginfo %>% 
  dplyr::filter(DrugID %in% d2BRCA$TTDDRUID)
cmap_cp = data.table::fread("CMAP_compoundinfo_beta.txt")
cmap_cp_BRCA = cmap_cp %>%
	dplyr::filter(inchi_key %in% drugBRCA_info$DRUGINKE) %>%
	dplyr::distinct(pert_id, inchi_key, compound_aliases)
sig_TTD = sig_info_sub[,c("pert_id","sig_id")] %>%
	dplyr::filter(pert_id %in% (cmap_cp %>%
		dplyr::filter(inchi_key %in% druginfo$DRUGINKE) %>% 
		dplyr::pull(pert_id)))
meta_cmap = data.frame(sig_id=sig_info_sub$sig_id,
					   BRCA = ifelse(sig_info_sub$pert_id %in% cmap_cp_BRCA$pert_id,
									1, 0))				
meta_cmap = meta_cmap[match(colnames(exp_cmap), meta_cmap$sig_id),] %>% 
				dplyr::left_join(sig_TTD[,c(1,2)]) %>% 
				dplyr::mutate(Anno=case_when(
					BRCA == 0 & is.na(pert_id) ~ "Lincs-remain",
					BRCA == 0 & !is.na(pert_id) ~ "TTD-remain",
					BRCA == 1 ~ "BRCA")) %>% dplyr::select(!c(pert_id, BRCA))
rownames(meta_cmap) = meta_cmap$sig_id


col_meta = read_gctx_meta("CMAP_level5_beta_trt_cp_n720216x12328.gctx", 
						   dim="col")
sig_info = data.table::fread("CMAP_siginfo_beta.txt",data.table=FALSE)
table(sig_info$cell_mfc_name %in% "MCF7")
sig_info_sub = sig_info %>% 
				dplyr::filter(pert_type == "trt_cp") %>%
				dplyr::filter(cell_mfc_name %in% "MCF7") %>%
				dplyr::filter(pert_itime == "24 h") %>%
				dplyr::arrange(pert_id, desc(tas)) %>% 
				dplyr::distinct(pert_id, .keep_all=T)
gene_meta = data.table::fread("CMAP_geneinfo_beta.txt")
gene_meta = subset(gene_meta, ensembl_id != "ENSG00000150526")


row_meta <- read_gctx_meta("CMAP_level5_beta_trt_cp_n720216x12328.gctx", 
						   dim="row")
gctx = parse_gctx("CMAP_level5_beta_trt_cp_n720216x12328.gctx",
				  cid=which(col_meta$id %in% sig_info_sub$sig_id))
exp_cmap = gctx@mat %>% as.data.frame() %>% 
 				dplyr::filter(rownames(.) %in% gene_meta$gene_id)
rownames(exp_cmap) = gene_meta$gene_symbol[match(rownames(exp_cmap),
													 gene_meta$gene_id)]

exp_cmap_FC1 = exp_cmap[,meta_cmap$sig_id[meta_cmap$Anno != "Lincs-remain"]]
exp_cmap_FC1[abs(exp_cmap_FC1)<0.5] = 0
exp_cmap_FC1_sign = sign(exp_cmap_FC1)
exp_cmap_FC1_sign_neg = -1*exp_cmap_FC1_sign
exp_cmap_FC1_sign_neg$gene = rownames(exp_cmap_FC1_sign_neg)
exp_cmap_FC1_sign_neg[1:4,1:4]
dim(exp_cmap_FC1_sign_neg)


ppi= data.table::fread("PPI_experiment_based.csv")[,c("SYMBOL_1","SYMBOL_2")] %>% na.omit()
weights = matrix(c(1, 0.5, 0,
				   0, 0, 0,
				   0, 0.5, 1),
				nrow=3, byrow = TRUE)
AC_ttd.list = mclapply(seq(pw_filt_3), function(k){
	gset_sub = subset(gset, term %in% pw_filt_3[k])
	genes = gset_sub$gene
	ppi_sub = subset(ppi, SYMBOL_1 %in% genes & SYMBOL_2 %in% genes)
	g1 = graph_from_data_frame(d = ppi_sub, directed = FALSE)
	g1_BW = igraph::betweenness(g1) %>% as.data.frame() %>% 
		tibble::rownames_to_column("gene") %>% 
		dplyr::rename("BW"=".") %>% 
		dplyr::mutate(BW_log = log2(BW+2)) %>% 
		dplyr::arrange(desc(BW_log))
	bp_genes_sign = gset %>% 
		dplyr::filter(term %in% pw_filt_3[k],
					  gene %in% rownames(exp_cmap)) %>% 
		dplyr::inner_join(g1_BW[,c("gene","BW_log")],by=c("gene"="gene")) %>%
		dplyr::inner_join(subset(deg, padj < 0.01 & abs(log2FoldChange)>0.5)[,c("X","log2FoldChange")],
						 by=c("gene"="X")) %>% 
		dplyr::rename(BRCA_log2FC=log2FoldChange) %>% 
		dplyr::mutate(BRCA_sign = sign(BRCA_log2FC)) %>% 
		dplyr::left_join(exp_cmap_FC1_sign_neg,by=c("gene"="gene")) %>% 
		dplyr::mutate(across(!c(term,gene,BW_log,BRCA_log2FC), ~ factor(.x, levels=c(-1,0,1)))) 
	AC_res = mclapply(seq(ncol(exp_cmap_FC1)), function(j){
				cp_sle = colnames(exp_cmap_FC1)[j]
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
				}, mc.cores=30) %>% unlist
	return(AC_res)
},mc.cores=5)


AC_ttd = do.call(cbind, AC_ttd.list)
colnames(AC_ttd) = pw_filt_3
rownames(AC_ttd) = colnames(exp_cmap_FC1)

AC_ttd_reshaped = reshape2::melt(AC_ttd) %>%
	dplyr::left_join(meta_cmap, by=c("Var1"="sig_id")) %>% 
	dplyr::filter(Anno!="Lincs-remain")


wilcox_stat = compare_means(value ~ Anno, data = AC_ttd_reshaped,
              group.by = "Var2", alternative = "greater")
wilcox_stat$Var2 = as.character(wilcox_stat$Var2) 
pw_filt_4 = wilcox_stat$Var2[wilcox_stat$p<0.01]


sig_info_submeta = meta_cmap %>%
				dplyr::select(sig_id, pert_mfc_id) %>%
				dplyr::left_join(cmap_cp[,c("inchi_key","pert_id")],
					 by=c("pert_mfc_id"="pert_id")) %>% na.omit() %>% 
				dplyr::left_join(meta_cmap) %>% distinct()

GDSC1 = read.csv("GDSC1_MCF7_drug_IC50.csv") %>%
	dplyr::group_by(inchi_key) %>%
	dplyr::summarise(IC50=mean(IC50))
GDSC1_ttd = sig_info_submeta %>% 
	dplyr::filter(inchi_key %in% unique(GDSC1$inchi_key)) %>% 
	dplyr::filter(Anno != "Lincs-remain") %>% 
	dplyr::left_join(GDSC1)


AC_ttd_GDSC1 = AC_ttd[GDSC1_ttd$sig_id, pw_filt_4]
cor_res = lapply(1:ncol(AC_ttd_GDSC1), function(i){
	cor_p = cor.test(AC_ttd_GDSC1[,i], GDSC1_ttd$IC50)$p.value
	cor_r = cor.test(AC_ttd_GDSC1[,i], GDSC1_ttd$IC50)$estimate
	return(c(cor_p, cor_r))
}) %>% do.call(rbind, .) %>% as.data.frame() %>%
	dplyr::rename(pval=V1) %>%
	dplyr::mutate(PW=pw_filt_4)
gdsc1_pw = cor_res %>%
	dplyr::filter(cor<0, pval<0.01) %>%
	dplyr::pull(PW)

GDSC2 = read.csv("GDSC2_MCF7_drug_IC50.csv") %>%
	dplyr::group_by(inchi_key) %>%
	dplyr::summarise(IC50=mean(IC50))
GDSC2_ttd = sig_info_submeta %>% 
	dplyr::filter(inchi_key %in% unique(GDSC2$inchi_key)) %>% 
	dplyr::filter(Anno != "Lincs-remain") %>% 
	dplyr::left_join(GDSC2)

AC_ttd_GDSC2 = AC_ttd[GDSC2_ttd$sig_id, pw_filt_4]
cor_res2 = lapply(1:ncol(AC_ttd_GDSC2), function(i){
	cor_p = cor.test(AC_ttd_GDSC2[,i], GDSC2_ttd$IC50)$p.value
	cor_r = cor.test(AC_ttd_GDSC2[,i], GDSC2_ttd$IC50)$estimate
	return(c(cor_p, cor_r))
}) %>% do.call(rbind, .) %>% as.data.frame() %>%
	dplyr::rename(pval=V1) %>%
	dplyr::mutate(PW=pw_filt_4)
gdsc2_pw = cor_res2 %>%
	dplyr::filter(cor<0, pval<0.01) %>%
	dplyr::pull(PW)
pw_filt_5 = intersect(gdsc1_pw, gdsc2_pw)


pw_filt_ok = pw_filt_5
gset_sub = subset(gset, term %in% pw_filt_ok)
gset_sub$term = as.character(gset_sub$term)
gset_sub.list = split(gset_sub$gene, gset_sub$term)
pathway_state = lapply(gset_sub.list, length) %>% unlist() %>% 
	reshape2::melt() %>% 
	dplyr::rename(size=value) %>% 
	tibble::rownames_to_column("Term")
pathway_state$DEGs = NA
pathway_state$DEG_UPs = NA
pathway_state$DEG_DOWNs = NA
pathway_state$Jacc_max = NA
pathway_state$Jacc_PW=NA

for(i in seq(nrow(pathway_state))){
	pw_sle = pathway_state$Term[i]
	jaccs = lapply(setdiff(pathway_state$Term, pw_sle), function(x){
		r1 = intersect(gset_sub.list[[pw_sle]],
				  gset_sub.list[[x]])
		r2 = length(r1)/length(gset_sub.list[[pw_sle]])
		return(r2)
	}) %>% unlist
	names(jaccs) = setdiff(pathway_state$Term, pw_sle)
	pw_deg = deg %>% 
		dplyr::filter(padj<0.01, abs(log2FoldChange)>0.5) %>% 
		dplyr::filter(X %in% gset_sub.list[[pw_sle]])
	pw_deg_up = subset(pw_deg, Direct=="Up")
	pw_deg_down = subset(pw_deg, Direct=="Down")
		
	pathway_state$DEGs[i]=nrow(pw_deg)
	pathway_state$DEG_UPs[i]=nrow(pw_deg_up)
	pathway_state$DEG_DOWNs[i]=nrow(pw_deg_down)
	pathway_state$Jacc_max[i]=max(jaccs)
	pathway_state$Jacc_PW[i]=names(jaccs)[jaccs==max(jaccs)]
}
pathway_state = pathway_state[match(pw_filt_ok, pathway_state$Term), ]
rownames(pathway_state) = 1:nrow(pathway_state)