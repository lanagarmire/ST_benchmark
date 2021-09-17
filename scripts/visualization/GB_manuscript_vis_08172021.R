rm(list=ls())

setwd("/home/liyijun/ST_benchmark_01082020")
library(data.table)

library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(gridExtra)
library(viridis)
library(gtools)
library(grDevices)
library(ggpubr)
library(gtable)
library(cowplot)

##### Figure 1: example of ST_MOB1 ################# 
if(!dir.exists("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig1")){
  dir.create("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig1")
}

load("data_analysis/real_data/data/ST_MOB1.RData")
spatial_annot = cbind(spatial, annotation)
ggplot(data = spatial_annot, aes(x=sdimx, y=sdimy, color = factor(group))) +
       geom_point(size=7) + scale_color_manual(values=manual_col_pallete[1:5]) +
       theme_bw() + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), legend.position="none") + 
       labs(x="",y="", title = "")
ggsave("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig1/ST_example_spat.png", width = 8, height = 8)

seed=12345
set.seed(seed)
rand_cells = sample(1:256, 25)
rand_genes = sample(1:16573, 50)
norm_expr_gather = norm_expr[rand_genes, rand_cells] %>% t() %>% as.matrix() %>% as_tibble() %>% rowid_to_column(var="cell") %>% gather(key="gene",value="exp",-1) %>% mutate(cell=paste("cell",cell,sep="_"))
head(norm_expr_gather)
ggplot(data = norm_expr_gather, aes(cell, gene, fill = exp)) + 
       geom_tile() + scale_fill_viridis(discrete=FALSE) +
       theme_bw() + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), legend.position="none") + 
       labs(x="", y="", title = "")
ggsave("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig1/ST_example_exp.png", width = 5, height = 8)



######## Figure 2: real dataset results #################
if(!dir.exists("data_analysis/real_data/vis_manuscript/Fig2_main")){
  dir.create("data_analysis/real_data/vis_manuscript/Fig2_main")
}

params_df = read.csv("data_analysis/real_data/real_data_names.csv", row.names = 1)
data_groups = list(seqfish = c(params_df$name[14]),
                   visium = c(params_df$name[c(1, 15:19)]),
                   merfish = c(params_df$name[2:13]))

metrics_sc = NULL

### load merFISH results 
load("data_analysis/real_data/merFISH/ndims_10/results/cluster/metrics_dim_10.RData")
metrics = metrics %>% filter(model != "HVG")
metrics_avg = metrics %>% group_by(data_name, model) %>% dplyr::summarise(ARI = mean(ARI), AMI = mean(AMI), DBI = mean(DBI), DBI_truth = mean(DBI_truth), 
									  AIC = mean(AIC), AIC0 = mean(AIC0), VN = mean(VN), TPR = mean(TPR), PPV = mean(PPV), FDR = mean(FDR), 
									  FM = mean(FM), CSI = mean(CSI), ACC = mean(ACC), F1 = mean(F1))
metrics_avg = metrics_avg %>% dplyr::group_by(data_name) %>%
                                      dplyr::mutate(ARI_rank = rank(-1.0*ARI, ties.method="first"),
                                                    AMI_rank = rank(-1.0*AMI, ties.method="first"),
                                                    DBI_rank = rank(DBI, ties.method="first"),
                                                    DBI_truth_rank = rank(DBI_truth, ties.method="first"),
                                                    AIC_rank = rank(-1.0*AIC, ties.method="first"),
                                                    AIC0_rank = rank(-1.0*AIC0, ties.method="first"),			
                                                    VN_rank = rank(-1.0*VN, ties.method="first"),
						    TPR_rank = rank(-1.0*TPR, ties.method = "first"),
						    PPV_rank = rank(-1.0*PPV, ties.method = "first"),
						    FDR_rank = rank(FDR, ties.method = "first"),
 						    FM_rank = rank(-1.0*FM, ties.method = "first"),
						    CSI_rank = rank(-1.0*CSI, ties.method = "first"),
						    ACC_rank = rank(-1.0*ACC, ties.method = "first"),
						    F1_rank = rank(-1.0*F1, ties.method = "first"))
metrics_avg$data_type = "merFISH"
metrics_sc = rbind(metrics_sc, metrics_avg)

### load SeqFISH+ results
load("data_analysis/real_data/ndims_10/results/cluster/metrics_dim_10.RData")
metrics = metrics %>% filter(data_name == "sfp_SScortex") %>% filter(model != "HVG")
metrics_avg = metrics %>% group_by(data_name, model) %>% dplyr::summarise(ARI = mean(ARI), AMI = mean(AMI), DBI = mean(DBI), DBI_truth = mean(DBI_truth), 
									  AIC = mean(AIC), AIC0 = mean(AIC0), VN = mean(VN), TPR = mean(TPR), PPV = mean(PPV), FDR = mean(FDR), 
									  FM = mean(FM), CSI = mean(CSI), ACC = mean(ACC), F1 = mean(F1))
metrics_avg = metrics_avg %>% dplyr::group_by(data_name) %>%
                                      dplyr::mutate(ARI_rank = rank(-1.0*ARI, ties.method="first"),
                                                    AMI_rank = rank(-1.0*AMI, ties.method="first"),
                                                    DBI_rank = rank(DBI, ties.method="first"),
                                                    DBI_truth_rank = rank(DBI_truth, ties.method="first"),
                                                    AIC_rank = rank(-1.0*AIC, ties.method="first"),
                                                    AIC0_rank = rank(-1.0*AIC0, ties.method="first"),			
                                                    VN_rank = rank(-1.0*VN, ties.method="first"),
						    TPR_rank = rank(-1.0*TPR, ties.method = "first"),
						    PPV_rank = rank(-1.0*PPV, ties.method = "first"),
						    FDR_rank = rank(FDR, ties.method = "first"),
 						    FM_rank = rank(-1.0*FM, ties.method = "first"),
						    CSI_rank = rank(-1.0*CSI, ties.method = "first"),
						    ACC_rank = rank(-1.0*ACC, ties.method = "first"),
						    F1_rank = rank(-1.0*F1, ties.method = "first"))
metrics_avg$data_type = "SeqFISH+"
metrics_sc = rbind(metrics_sc, metrics_avg)
metrics_sc$data_class = "single-cell"
save(metrics_sc, file = fs::path("data_analysis/real_data/vis_manuscript", "metrics_sc", ext="RData"))

metrics_non_sc = NULL

### load ST_MOB1 results
load("data_analysis/real_data/ndims_5/results/cluster/metrics_dim_5.RData")
metrics = metrics %>% filter(model != "HVG") %>% filter(data_name == "ST_MOB1")
metrics_avg = metrics %>% group_by(data_name, model) %>% dplyr::summarise(ARI = mean(ARI), AMI = mean(AMI), DBI = mean(DBI), DBI_truth = mean(DBI_truth),
									  AIC = mean(AIC), AIC0 = mean(AIC0), VN = mean(VN), TPR = mean(TPR), PPV = mean(PPV), FDR = mean(FDR), 
									  FM = mean(FM), CSI = mean(CSI), ACC = mean(ACC), F1 = mean(F1))

metrics_avg = metrics_avg %>% dplyr::group_by(data_name) %>%
                                      dplyr::mutate(ARI_rank = rank(-1.0*ARI, ties.method="first"),
                                                    AMI_rank = rank(-1.0*AMI, ties.method="first"),
                                                    DBI_rank = rank(DBI, ties.method="first"),
                                                    DBI_truth_rank = rank(DBI_truth, ties.method="first"),
                                                    AIC_rank = rank(-1.0*AIC, ties.method="first"),
                                                    AIC0_rank = rank(-1.0*AIC0, ties.method="first"),			
                                                    VN_rank = rank(-1.0*VN, ties.method="first"),
						    TPR_rank = rank(-1.0*TPR, ties.method = "first"),
						    PPV_rank = rank(-1.0*PPV, ties.method = "first"),
						    FDR_rank = rank(FDR, ties.method = "first"),
						    FM_rank = rank(-1.0*FM, ties.method = "first"),
						    CSI_rank = rank(-1.0*CSI, ties.method = "first"),
						    ACC_rank = rank(-1.0*ACC, ties.method = "first"),
						    F1_rank = rank(-1.0*F1, ties.method = "first"))
#metrics_avg$data_type = "ST"
metrics_avg$data_type = "Visium"
metrics_non_sc = rbind(metrics_non_sc, metrics_avg)

### load Visium results
load("data_analysis/real_data/ndims_15/results/cluster/metrics_dim_15.RData")
metrics = metrics %>% filter(data_name %in% data_groups$visium) %>% filter(model != "HVG") %>% filter(data_name != "ST_MOB1")
metrics_avg = metrics %>% group_by(data_name, model) %>% dplyr::summarise(ARI = mean(ARI), AMI = mean(AMI), DBI = mean(DBI), DBI_truth = mean(DBI_truth), 
									  AIC = mean(AIC), AIC0 = mean(AIC0), VN = mean(VN), TPR = mean(TPR), PPV = mean(PPV), FDR = mean(FDR), 
									  FM = mean(FM), CSI = mean(CSI), ACC = mean(ACC), F1 = mean(F1))

metrics_avg = metrics_avg %>% dplyr::group_by(data_name) %>%
                                      dplyr::mutate(ARI_rank = rank(-1.0*ARI, ties.method="first"),
                                                    AMI_rank = rank(-1.0*AMI, ties.method="first"),
                                                    DBI_rank = rank(DBI, ties.method="first"),
                                                    DBI_truth_rank = rank(DBI_truth, ties.method="first"),
                                                    AIC_rank = rank(-1.0*AIC, ties.method="first"),
                                                    AIC0_rank = rank(-1.0*AIC0, ties.method="first"),			
                                                    VN_rank = rank(-1.0*VN, ties.method="first"),
						    TPR_rank = rank(-1.0*TPR, ties.method = "first"),
						    PPV_rank = rank(-1.0*PPV, ties.method = "first"),
						    FDR_rank = rank(FDR, ties.method = "first"),
 						    FM_rank = rank(-1.0*FM, ties.method = "first"),
						    CSI_rank = rank(-1.0*CSI, ties.method = "first"),
						    ACC_rank = rank(-1.0*ACC, ties.method = "first"),
						    F1_rank = rank(-1.0*F1, ties.method = "first"))
metrics_avg$data_type = "Visium"
metrics_non_sc = rbind(metrics_non_sc, metrics_avg)
metrics_non_sc$data_class = "non-single-cell"
save(metrics_non_sc, file = fs::path("data_analysis/real_data/vis_manuscript", "metrics_non_sc", ext="RData"))

metrics_all = rbind(metrics_sc, metrics_non_sc)

#### get KS p values (SV genes < method)
p_value_label = function(x){
  if(x>=0 & x<0.001){
    res = "***"
  }else if(x>=0.001 & x<0.01){
    res = "**"
  }else if(x>=0.01 & x<0.05){
    res = "*"
  }else if(x>=0.05 & x<0.1){
    res = "•" #Alt+0149
  }else if(x>=0.1 & x<=1){
    res=""
  }
 
  return(res)
}

methods_order = c("SG", "CIMLR", "concatenation", "MOFA+", "scVI", "SNF", "WNN")
#AMI_KS_pval = VN_KS_pval = F1_KS_pval = AMI_KS_pval_label = VN_KS_pval_label = F1_KS_pval_label = rep(NA, length(methods_order))
AMI_W_pval = VN_W_pval = F1_W_pval = AMI_W_pval_label = VN_W_pval_label = F1_W_pval_label = rep(NA, length(methods_order))

SG_metrics = metrics_all %>% filter(model == "SG") %>% select(data_name, model, AMI, VN, F1)
for(i in 2:length(methods_order)){
  method_metrics = metrics_all %>% filter(model == methods_order[i]) %>% select(data_name, model, AMI, VN, F1)
  #AMI_KS = ks.test(x = SG_metrics$AMI, y=method_metrics$AMI, alternative = "greater")
  #AMI_KS_pval[i] = AMI_KS$p.value
  #AMI_KS_pval_label[i] = p_value_label(AMI_KS_pval[i])

  AMI_W = wilcox.test(x = SG_metrics$AMI, y=method_metrics$AMI, alternative = "two.sided", paired=T)
  AMI_W_pval[i] = AMI_W$p.value
  AMI_W_pval_label[i] = p_value_label(AMI_W_pval[i])

  #VN_KS = ks.test(x = SG_metrics$VN, y=method_metrics$VN, alternative = "greater")
  #VN_KS_pval[i] = VN_KS$p.value
  #VN_KS_pval_label[i] = p_value_label(VN_KS_pval[i])

  VN_W = wilcox.test(x = SG_metrics$VN, y=method_metrics$VN, alternative = "two.sided", paired=T)
  VN_W_pval[i] = VN_W$p.value
  VN_W_pval_label[i] = p_value_label(VN_W_pval[i])  

  #F1_KS = ks.test(x = SG_metrics$F1, y=method_metrics$F1, alternative = "greater")
  #F1_KS_pval[i] = F1_KS$p.value
  #F1_KS_pval_label[i] = p_value_label(F1_KS_pval[i])

  F1_W = wilcox.test(x = SG_metrics$F1, y=method_metrics$F1, alternative = "two.sided", paired=T)
  F1_W_pval[i] = F1_W$p.value
  F1_W_pval_label[i] = p_value_label(F1_W_pval[i])
}

#model_pval = data.frame(model = methods_order, AMI_KS_pval, AMI_KS_pval_label, VN_KS_pval, VN_KS_pval_label, F1_KS_pval, F1_KS_pval_label, AMI_W_pval, AMI_W_pval_label, VN_W_pval, VN_W_pval_label, F1_W_pval, F1_W_pval_label)
model_pval = data.frame(model = methods_order, AMI_W_pval, AMI_W_pval_label, VN_W_pval, VN_W_pval_label, F1_W_pval, F1_W_pval_label)
#model_pval[1,c("AMI_KS_pval_label", "VN_KS_pval_label", "F1_KS_pval_label", "AMI_W_pval_label", "VN_W_pval_label", "F1_W_pval_label")] = ""
model_pval[1,c("AMI_W_pval_label", "VN_W_pval_label", "F1_W_pval_label")] = ""
metrics_all = full_join(metrics_all, model_pval, by="model")

save(metrics_all, file = fs::path("data_analysis/real_data/vis_manuscript", "metrics_all", ext="RData"))

main_data = c("mf_hypo_60", "mf_hypo_160", "mf_hypo_-240", data_groups$seqfish, data_groups$visium)
supp_data = c("mf_hypo_-140", "mf_hypo_-190", "mf_hypo_-290", "mf_hypo_-40", "mf_hypo_-90", "mf_hypo_10", "mf_hypo_110", "mf_hypo_210", "mf_hypo_260")

metrics_all_main = metrics_all %>% filter(data_name %in% main_data)
metrics_all_supp = metrics_all %>% filter(data_name %in% supp_data)
save(metrics_all_main, file = fs::path("data_analysis/real_data/vis_manuscript", "metrics_all_main", ext="RData"))
save(metrics_all_supp, file = fs::path("data_analysis/real_data/vis_manuscript", "metrics_all_supp", ext="RData"))

###### Figure 2: barplots
manual_col_pallete0 = c("WNN" = "aquamarine3", "concatenation" = "coral4", "SG"= "black", "MOFA+" = "blue", "CIMLR" = "blueviolet", "scVI" = "darkgoldenrod1", "SNF" = "darkgreen")

############## main 
if(!dir.exists("data_analysis/real_data/vis_manuscript/Fig2_main")){
  dir.create("data_analysis/real_data/vis_manuscript/Fig2_main")
}

methods_labels = c("SV genes", "CIMLR", "concatenation", "MOFA+", "scVI", "SNF", "Seurat v4")

#AMI main
plot_main_list = list()
for(i in 1:length(main_data)){
  model_levels_main_AMI = metrics_all_main  %>% filter(model != "BREM-SC") %>% filter(data_name %in% main_data[i]) %>% arrange(desc(AMI)) %>% pull(model)
  SG_idx = which(model_levels_main_AMI == "SG")
  model_levels_main_AMI = c("SG", model_levels_main_AMI[-SG_idx])
  plot_main_list[[i]] = metrics_all_main %>% filter(model != "BREM-SC") %>% filter(data_name %in% main_data[i]) %>% mutate(model = factor(model, levels = model_levels_main_AMI)) %>%
                        ggplot(aes(x=data_name, y=AMI, fill=model)) +
                        geom_bar(position="dodge", stat="identity") +
                        scale_fill_manual(values = manual_col_pallete0, breaks = methods_order, labels = methods_labels) +
                        labs(x="", y="", title="")+ ylim(0,1)+
                        theme_bw() + theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.position="none", 
			axis.text.x = element_text(size = 5.5), plot.margin = unit(c(0, 0, 0, 0),"cm"))
#                        geom_text(aes(label = AMI_W_pval_label), position = position_dodge(0.9), vjust = -0.5, size = 1.5)
}
plot_main_list[[1]] = plot_main_list[[1]]+theme(axis.ticks.y=element_line(), axis.text.y=element_text(size=7))

legend_plot = plot_main_list[[10]]+theme(legend.position="bottom")+scale_fill_manual(values = manual_col_pallete0, breaks = methods_order, labels = methods_labels) +
                                                                           theme(legend.text=element_text(size=7), legend.title = element_text(size = 8)) + guides(fill = guide_legend(nrow=1))
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig2_main", "fig2_legend", ext="png"), legend_plot, width = 210, height = 58, units = "mm")

AMI_real_data_main = grid.arrange(plot_main_list[[1]], plot_main_list[[2]], plot_main_list[[3]], plot_main_list[[4]], plot_main_list[[5]], plot_main_list[[6]], plot_main_list[[7]], plot_main_list[[8]], plot_main_list[[9]], plot_main_list[[10]], widths = c(1.2, rep(1,9)), ncol=10)
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig2_main", "AMI_real_data_main", ext="png"), AMI_real_data_main, width = 212, height = 61, units = "mm")

#VN main
plot_main_list = list()
for(i in 1:length(main_data)){
  model_levels_main_VN = metrics_all_main %>% filter(model != "BREM-SC") %>% filter(data_name %in% main_data[i]) %>% arrange(desc(VN)) %>% pull(model)
  SG_idx = which(model_levels_main_VN == "SG")
  model_levels_main_VN = c("SG", model_levels_main_VN[-SG_idx])
  plot_main_list[[i]] = metrics_all_main %>% filter(model != "BREM-SC") %>% filter(data_name %in% main_data[i]) %>% mutate(model = factor(model, levels = model_levels_main_VN)) %>%
                        ggplot(aes(x=data_name, y=VN, fill=model)) +
                        geom_bar(position="dodge", stat="identity") +
                        scale_fill_manual(values = manual_col_pallete0, breaks = methods_order, labels = methods_labels) +
                        labs(x="", y="", title="")+ ylim(0,1)+
                        theme_bw() + theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.position="none", 
			axis.text.x = element_text(size = 5.5), plot.margin = unit(c(0, 0, 0, 0),"cm"))
#			geom_text(aes(label = VN_W_pval_label), position = position_dodge(0.9), vjust = -0.5, size = 1.5)
}
plot_main_list[[1]] = plot_main_list[[1]]+theme(axis.ticks.y=element_line(), axis.text.y=element_text(size=7))
VN_real_data_main = grid.arrange(plot_main_list[[1]], plot_main_list[[2]], plot_main_list[[3]], plot_main_list[[4]], plot_main_list[[5]], plot_main_list[[6]], plot_main_list[[7]], plot_main_list[[8]], plot_main_list[[9]], plot_main_list[[10]], widths = c(1.2, rep(1,9)), ncol=10)
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig2_main", "VN_real_data_main", ext="png"), VN_real_data_main, width = 212, height = 61, units = "mm")

#F1 main
plot_main_list = list()
for(i in 1:length(main_data)){
  model_levels_main_F1 = metrics_all_main %>% filter(model != "BREM-SC") %>% filter(data_name %in% main_data[i]) %>% arrange(desc(F1)) %>% pull(model)
  SG_idx = which(model_levels_main_F1 == "SG")
  model_levels_main_F1 = c("SG", model_levels_main_F1[-SG_idx])
  plot_main_list[[i]] = metrics_all_main %>% filter(model != "BREM-SC") %>% filter(data_name %in% main_data[i]) %>% mutate(model = factor(model, levels = model_levels_main_F1)) %>%
                        ggplot(aes(x=data_name, y=F1, fill=model)) +
                        geom_bar(position="dodge", stat="identity") +
                        scale_fill_manual(values = manual_col_pallete0, breaks = methods_order, labels = methods_labels) +
                        labs(x="", y="", title="")+ ylim(0,1)+
                        theme_bw() + theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.position="none", 
			axis.text.x = element_text(size = 5.5), plot.margin = unit(c(0, 0, 0, 0),"cm"))
#                        geom_text(aes(label = F1_W_pval_label), position = position_dodge(0.9), vjust = -0.5, size = 1.5)+
}
plot_main_list[[1]] = plot_main_list[[1]]+theme(axis.ticks.y=element_line(), axis.text.y=element_text(size=7))
F1_real_data_main = grid.arrange(plot_main_list[[1]], plot_main_list[[2]], plot_main_list[[3]], plot_main_list[[4]], plot_main_list[[5]], plot_main_list[[6]], plot_main_list[[7]], plot_main_list[[8]], plot_main_list[[9]], plot_main_list[[10]], widths = c(1.2, rep(1,9)), ncol=10)
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig2_main", "F1_real_data_main", ext="png"), F1_real_data_main, width = 212, height = 61, units = "mm")



############# Figure 2: detailed heatmaps
# AMI
metrics_all_main %>% filter(model != "BREM-SC") %>%
ggplot(aes(x=factor(model, level = methods_order), y=data_name, fill=AMI))+
       geom_tile()+
       geom_text(aes(label = AMI_rank), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = -1, 
            breaks = c(min(metrics_all_main$AMI),max(metrics_all_main$AMI)), labels=c("7","1"))+
       scale_x_discrete(position="top", breaks = methods_order, labels = methods_labels)+
       guides(colour=F)+
       labs(x="", y="", fill="rank")+
       theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 8, hjust=0.05, vjust=0.2), axis.text.y = element_text(size=8), 
        legend.title = element_text(size=8), legend.text = element_text(size=8))
ggsave(fs::path("data_analysis/real_data/vis_manuscript/Fig2_main", "AMI_heatmap", ext="png"), width = 85, height = 76, units = "mm") #change this

# VN
metrics_all_main %>% filter(model != "BREM-SC") %>%
ggplot(aes(x=factor(model, level = methods_order), y=data_name, fill=VN))+
       geom_tile()+
       geom_text(aes(label = VN_rank), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = -1, 
            breaks = c(min(metrics_all_main$VN),max(metrics_all_main$VN)), labels=c("7","1"))+
       scale_x_discrete(position="top", breaks = methods_order, labels = methods_labels)+
       guides(colour=F)+
       labs(x="", y="", fill="rank")+
       theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 8, hjust=0.05, vjust=0.2), axis.text.y = element_text(size=8), 
        legend.title = element_text(size=8), legend.text = element_text(size=8))
ggsave(fs::path("data_analysis/real_data/vis_manuscript/Fig2_main", "VN_heatmap", ext="png"), width = 85, height = 76, units = "mm") #change this

# F1
metrics_all_main %>% filter(model != "BREM-SC") %>%
ggplot(aes(x=factor(model, level = methods_order), y=data_name, fill=F1))+
       geom_tile()+
       geom_text(aes(label = F1_rank), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = -1, 
            breaks = c(min(metrics_all_main$F1),max(metrics_all_main$F1)), labels=c("7","1"))+
       scale_x_discrete(position="top", breaks = methods_order, labels = methods_labels)+
       guides(colour=F)+
       labs(x="", y="", fill="rank")+
       theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 8, hjust=0.05, vjust=0.2), axis.text.y = element_text(size=8), 
        legend.title = element_text(size=8), legend.text = element_text(size=8))
ggsave(fs::path("data_analysis/real_data/vis_manuscript/Fig2_main", "F1_heatmap", ext="png"), width = 85, height = 76, units = "mm") #change this



########### Figure2: boxplots
### all (main + supp)
#AMI 
metrics_all %>% filter(model != "BREM-SC") %>% mutate(model = factor(model, levels = methods_order, ordered = T)) %>%
ggplot(aes(x=model, y=AMI, fill=model)) + 
       geom_boxplot(outlier.shape = NA)+
       scale_fill_manual(values = manual_col_pallete0)+
       scale_x_discrete(breaks = methods_order, labels = methods_labels)+
       theme_bw()+
       labs(y="", x="", title="") + ylim(0,1)+
       theme(legend.position="none",  axis.text.x = element_text(size=8, angle=30, vjust = 0.4), axis.text.y = element_text(size=8))+
       stat_compare_means(method = "wilcox", paired = T, method.args = list(alternative = "two.sided"), ref.group = "SG", hide.ns = T, label = "p.signif", vjust = 1,
	symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig2_main", "AMI_real_data_all_boxplot", ext="png"), width = 85, height = 85, units= "mm")

#VN
metrics_all %>% filter(model != "BREM-SC") %>% mutate(model = factor(model, levels = methods_order, ordered = T)) %>%
ggplot(aes(x=model, y=VN, fill=model)) + 
       geom_boxplot(outlier.shape = NA)+
       scale_fill_manual(values = manual_col_pallete0)+
       scale_x_discrete(breaks = methods_order, labels = methods_labels)+
       theme_bw()+
       labs(y="", x="", title="") + ylim(0,1)+
       theme(legend.position="none",  axis.text.x = element_text(size=8, angle=30, vjust = 0.4), axis.text.y = element_text(size=8))+
       stat_compare_means(method = "wilcox", paired = T, method.args = list(alternative = "two.sided"), ref.group = "SG", hide.ns = T, label = "p.signif", vjust = 1,
	symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig2_main", "VN_real_data_all_boxplot", ext="png"), width = 85, height = 85, units= "mm")

#F1
metrics_all %>% filter(model != "BREM-SC") %>% mutate(model = factor(model, levels = methods_order, ordered = T)) %>%
ggplot(aes(x=model, y=F1, fill=model)) + 
       geom_boxplot(outlier.shape = NA)+
       scale_fill_manual(values = manual_col_pallete0)+
       scale_x_discrete(breaks = methods_order, labels = methods_labels)+
       theme_bw()+
       labs(y="", x="", title="") + ylim(0,1)+
       theme(legend.position="none",  axis.text.x = element_text(size=8, angle=30, vjust = 0.4), axis.text.y = element_text(size=8))+
       stat_compare_means(method = "wilcox", paired = T, method.args = list(alternative = "two.sided"), ref.group = "SG", hide.ns = T, label = "p.signif", vjust = 0.5,
	symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig2_main", "F1_real_data_all_boxplot", ext="png"), width = 85, height = 85, units= "mm")



############### Figure 4: simulation results #########################
if(!dir.exists("simulation/simulation_07202021/ST_MOB1/vis_manuscript/Fig4")){
  dir.create("simulation/simulation_07202021/ST_MOB1/vis_manuscript/Fig4")
}

R_models = c("MOFA+","scVI")
ndims = c(5, 15, 30)
#ndims = c(5)
sim_data_path = "simulation/simulation_07202021"
sim_data_ref = "ST_MOB1"

sim_metrics_avg_all = sim_metrics_reps_all = NULL

for(i in 1:length(ndims)){
  sim_metrics_path = fs::path(sim_data_path, sim_data_ref, paste("ndims", ndims[i], sep="_")) #change this
  load(fs::path(sim_metrics_path, "results", "cluster", paste("metrics_dim", ndims[i], "k", 30, sep="_"), ext = "RData")) #change here sometimes
  
  metrics = metrics %>% filter(model != "HVG") %>% filter(spat_prob %in% c(0.6, 0.7, 0.8, 0.9))

  #stochasticity
  sim_metrics_avg = metrics %>% group_by(spat_prob, rep_id, model) %>% dplyr::summarise(ARI = mean(ARI), AMI = mean(AMI), DBI = mean(DBI), 
                                                                          DBI_truth = mean(DBI_truth), AIC = mean(AIC), AIC0 = mean(AIC0), VN = mean(VN), 
                                                                          TPR = mean(TPR), PPV = mean(PPV), FDR = mean(FDR), FM = mean(FM), 
									  CSI = mean(CSI), ACC = mean(ACC), F1 = mean(F1))
  
  #account for repeated experiments
  sim_metrics_avg_rep = sim_metrics_avg %>% group_by(spat_prob, model) %>% dplyr::summarise(ARI_mean = mean(ARI), ARI_sd = sd(ARI), 
                                                                                    AMI_mean = mean(AMI), AMI_sd = sd(AMI),
                                                                                    DBI_mean = mean(DBI), DBI_sd = sd(DBI),
                                                                                    DBI_truth_mean = mean(DBI_truth), DBI_truth_sd = sd(DBI_truth),
                                                                                    AIC_mean = mean(AIC), AIC_sd = sd(AIC),
                                                                                    AIC0_mean = mean(AIC0), AIC0_sd = sd(AIC0),
                                                                                    VN_mean = mean(VN), VN_sd = sd(VN),
                                                                                    TPR_mean = mean(TPR), TPR_sd = sd(TPR),
                                                                                    PPV_mean = mean(PPV), PPV_sd = sd(PPV),
                                                                                    FDR_mean = mean(FDR), FDR_sd = sd(FDR),
                                                                                    FM_mean = mean(FM), FM_sd = sd(FM),
                                                                                    CSI_mean = mean(CSI), CSI_sd = sd(CSI),
                                                                                    ACC_mean = mean(ACC), ACC_sd = sd(ACC),
                                                                                    F1_mean = mean(F1), F1_sd = sd(F1))
  sim_metrics_avg_rep = sim_metrics_avg_rep %>% dplyr::group_by(spat_prob) %>%
                                        dplyr::mutate(ARI_rank = rank(-1.0*ARI_mean, ties.method="first"),
                                                    AMI_rank = rank(-1.0*AMI_mean, ties.method="first"),
                                                    DBI_rank = rank(DBI_mean, ties.method="first"),
                                                    DBI_truth_rank = rank(DBI_truth_mean, ties.method="first"),
                                                    AIC_rank = rank(-1.0*AIC_mean, ties.method="first"),
                                                    AIC0_rank = rank(-1.0*AIC0_mean, ties.method="first"),			
                                                    VN_rank = rank(-1.0*VN_mean, ties.method="first"),
						    TPR_rank = rank(-1.0*TPR_mean, ties.method = "first"),
						    PPV_rank = rank(-1.0*PPV_mean, ties.method = "first"),
						    FDR_rank = rank(-1.0*FDR_mean, ties.method = "first"),
 						    FM_rank = rank(-1.0*FM_mean, ties.method = "first"),
						    CSI_rank = rank(-1.0*CSI_mean, ties.method = "first"),
						    ACC_rank = rank(-1.0*ACC_mean, ties.method = "first"),
						    F1_rank = rank(-1.0*F1_mean, ties.method = "first"))
  sim_metrics_avg_rep$dims = ndims[i]

  #concatenate metrics
  sim_metrics_avg_all = rbind(sim_metrics_avg_all, sim_metrics_avg) #adjusted for stochastic runs
  sim_metrics_reps_all = rbind(sim_metrics_reps_all, sim_metrics_avg_rep) #adjusted for stochastic and repeated runs 
}

sim_AMI_W_pval = sim_VN_W_pval = sim_F1_W_pval = sim_AMI_W_pval_label = sim_VN_W_pval_label = sim_F1_W_pval_label = rep(NA, length(methods_order))

sim_SG_metrics = sim_metrics_reps_all %>% filter(model == "SG") %>% select(spat_prob, model, AMI_mean, VN_mean, F1_mean)
for(i in 2:length(methods_order)){
  sim_method_metrics = sim_metrics_reps_all %>% filter(model == methods_order[i]) %>% select(spat_prob, model, AMI_mean, VN_mean, F1_mean)
  sim_AMI_W = wilcox.test(x = sim_SG_metrics$AMI_mean, y=sim_method_metrics$AMI_mean, alternative = "two.sided", paired=T)
  sim_AMI_W_pval[i] = sim_AMI_W$p.value
  sim_AMI_W_pval_label[i] = p_value_label(sim_AMI_W_pval[i])

  sim_VN_W = wilcox.test(x = sim_SG_metrics$VN_mean, y=sim_method_metrics$VN_mean, alternative = "two.sided", paired=T)
  sim_VN_W_pval[i] = sim_VN_W$p.value
  sim_VN_W_pval_label[i] = p_value_label(sim_VN_W_pval[i])  

  sim_F1_W = wilcox.test(x = sim_SG_metrics$F1_mean, y=sim_method_metrics$F1_mean, alternative = "two.sided", paired=T)
  sim_F1_W_pval[i] = sim_F1_W$p.value
  sim_F1_W_pval_label[i] = p_value_label(sim_F1_W_pval[i])
}

sim_model_pval = data.frame(model = methods_order, sim_AMI_W_pval, sim_AMI_W_pval_label, sim_VN_W_pval, sim_VN_W_pval_label, sim_F1_W_pval, sim_F1_W_pval_label)
sim_model_pval[1,c("sim_AMI_W_pval_label", "sim_VN_W_pval_label", "sim_F1_W_pval_label")] = ""

sim_metrics_reps_all = full_join(sim_metrics_reps_all, sim_model_pval, by="model")

save(sim_metrics_avg_all, sim_metrics_reps_all, file = "simulation/simulation_07202021/ST_MOB1/vis_manuscript/sim_metrics.RData")




################### Figure 4: line plots ######################
# AMI for ISCB poster
sim_metrics_reps_all %>% filter(dims == 5) %>%
ggplot(aes(x = spat_prob, y = AMI_mean, color = factor(model, levels=methods_order)))+
       geom_line() + geom_point()+
       geom_errorbar(aes(ymin=AMI_mean-AMI_sd, ymax=AMI_mean+AMI_sd), width = 0.02)+
       scale_color_manual(values = manual_col_pallete0, labels = methods_labels)+
       labs(x = "p", y="", color = "model", title = "")+
       theme_bw() +
       theme(axis.text.x=element_text(size=5), axis.text.y=element_text(size=5), legend.position="bottom", axis.title.x = element_text(size = 6),
             legend.text=element_text(size=6), legend.title = element_text(size = 6))
#       guides(colour = guide_legend(nrow = 1))
ggsave(fs::path("simulation/simulation_07202021/ST_MOB1/vis_manuscript/Fig4", "AMI_sim_ISCB", ext="png"), width = 105, height = 90, units = "mm") #change this

# VN
sim_metrics_reps_all %>% filter(dims == 5) %>%
  ggplot(aes(x = spat_prob, y = VN_mean, color = factor(model, levels=methods_order)))+
  geom_line() + geom_point()+
  geom_errorbar(aes(ymin=VN_mean-VN_sd, ymax=VN_mean+VN_sd), width = 0.02)+
  scale_color_manual(values = manual_col_pallete0, labels = methods_labels)+
  labs(x = "p", y="", color = "model", title = "")+
  theme_bw() +
  theme(axis.text.x=element_text(size=5), axis.text.y=element_text(size=5), legend.position="bottom", axis.title.x = element_text(size = 6),
             legend.text=element_text(size=6), legend.title = element_text(size = 6))
ggsave(fs::path("simulation/simulation_07202021/ST_MOB1/vis_manuscript/Fig4", "VN_sim", ext="png"), width = 105, height = 90, units = "mm") #change this

# F1
sim_metrics_reps_all %>% filter(dims == 5) %>%
  ggplot(aes(x = spat_prob, y = F1_mean, color = factor(model, levels=methods_order)))+
  geom_line() + geom_point()+
  geom_errorbar(aes(ymin=F1_mean-F1_sd, ymax=F1_mean+F1_sd), width = 0.02)+
  scale_color_manual(values = manual_col_pallete0, labels = methods_labels)+
  labs(x = "p", y="", color = "model", title = "")+
  theme_bw() +
  theme(axis.text.x=element_text(size=5), axis.text.y=element_text(size=5), legend.position="right", axis.title.x = element_text(size = 8),
             legend.text=element_text(size=8), legend.title = element_text(size = 8))
ggsave(fs::path("simulation/simulation_07202021/ST_MOB1/vis_manuscript/Fig4", "F1_sim", ext="png"), width = 140, height = 74.5, units = "mm") #change this



################# Fig_6: runtime and memory line plots ##############
## load cell runtime and memory
load("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/runtime_cells_df.RData")
load("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/memory_cells_df.RData")
load("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/runtime_features_df.RData")
load("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/memory_features_df.RData")

runtime_cells_df %>% mutate(method = factor(method, levels=methods_order)) %>%
ggplot(aes(x=num_cells, y=runtime, group = method, color = method))+
       geom_line() + geom_point()+
       scale_color_manual(values = manual_col_pallete0, breaks = methods_order, labels = methods_labels)+
       ylim(0, 350) +
       labs(x = "number of cells", y="", color = "model", title = "")+
       theme_bw()+
       theme(axis.text.x=element_text(size=5), legend.text=element_text(size=6), legend.title = element_text(size = 6), legend.position = "bottom", axis.title.x = element_text(size = 6))
ggsave(fs::path("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig5", "runtime_cells", ext="png"), width = 105, height = 100, units = "mm")

runtime_features_df %>% mutate(method = factor(method, levels=methods_order)) %>%
ggplot(aes(x=num_features, y=runtime, group = method, color = method))+
       geom_line() + geom_point()+
       scale_color_manual(values = manual_col_pallete0, breaks = methods_order, labels = methods_labels)+
       ylim(0, 350) +
       labs(x = "number of features", y="", color = "model", title = "")+
       theme_bw()+
       theme(axis.text.x=element_text(size=5), legend.text=element_text(size=6), legend.title = element_text(size = 6), legend.position = "bottom", axis.title.x = element_text(size = 6))
ggsave(fs::path("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig5", "runtime_features", ext="png"), width = 105, height = 100, units = "mm")

memory_cells_df %>% mutate(method = factor(method, levels=methods_order)) %>%
ggplot(aes(x=num_cells, y=memory, group = method, color = method))+
       geom_line() + geom_point()+
       scale_color_manual(values = manual_col_pallete0, breaks = methods_order, labels = methods_labels)+
       ylim(0, 35000) +
       labs(x = "number of cells", y="", color = "model", title = "")+
       theme_bw()+
       theme(axis.text.x=element_text(size=5), legend.text=element_text(size=6), legend.title = element_text(size = 6), legend.position = "bottom", axis.title.x = element_text(size = 6))
ggsave(fs::path("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig5", "memory_cells", ext="png"), width = 105, height = 100, units = "mm")

memory_features_df %>% mutate(method = factor(method, levels=methods_order)) %>%
ggplot(aes(x=num_features, y=memory, group = method, color = method))+
       geom_line() + geom_point()+
       scale_color_manual(values = manual_col_pallete0, breaks = methods_order, labels = methods_labels)+
       ylim(0, 35000) +
       labs(x = "number of features", y="", color = "model", title = "")+
       theme_bw()+
       theme(axis.text.x=element_text(size=5), legend.text=element_text(size=6), legend.title = element_text(size = 6), legend.position = "bottom", axis.title.x = element_text(size = 6))
ggsave(fs::path("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig5", "memory_features", ext="png"), width = 105, height = 100, units = "mm")



############# Figure 7: heatmaps #######################
if(!dir.exists("data_analysis/real_data/vis_manuscript/Fig7")){
  dir.create("data_analysis/real_data/vis_manuscript/Fig7")
}

##################### real data
metrics_all_tech_sum = metrics_all %>% filter(model != "BREM-SC") %>%
group_by(data_type, model) %>% dplyr::summarise(ARI = mean(ARI), AMI = mean(AMI), DBI = mean(DBI), DBI_truth = mean(DBI_truth), 
									  AIC = mean(AIC), AIC0 = mean(AIC0), VN = mean(VN), TPR = mean(TPR), PPV = mean(PPV), FDR = mean(FDR), 
									  FM = mean(FM), CSI = mean(CSI), ACC = mean(ACC), F1 = mean(F1))
metrics_all_tech_sum = metrics_all_tech_sum %>% dplyr::group_by(data_type) %>%
                                      dplyr::mutate(ARI_rank = rank(-1.0*ARI, ties.method="first"),
                                                    AMI_rank = rank(-1.0*AMI, ties.method="first"),
                                                    DBI_rank = rank(DBI, ties.method="first"),
                                                    DBI_truth_rank = rank(DBI_truth, ties.method="first"),
                                                    AIC_rank = rank(-1.0*AIC, ties.method="first"),
                                                    AIC0_rank = rank(-1.0*AIC0, ties.method="first"),			
                                                    VN_rank = rank(-1.0*VN, ties.method="first"),
						    TPR_rank = rank(-1.0*TPR, ties.method = "first"),
						    PPV_rank = rank(-1.0*PPV, ties.method = "first"),
						    FDR_rank = rank(FDR, ties.method = "first"),
 						    FM_rank = rank(-1.0*FM, ties.method = "first"),
						    CSI_rank = rank(-1.0*CSI, ties.method = "first"),
						    ACC_rank = rank(-1.0*ACC, ties.method = "first"),
						    F1_rank = rank(-1.0*F1, ties.method = "first"))

# AMI: real data
metrics_all_tech_sum %>% 
ggplot(aes(x=factor(model, level = methods_order), y=data_type, fill=AMI))+
       geom_tile()+
       geom_text(aes(label = AMI_rank), size = 3)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = -1, 
            breaks = c(min(metrics_all_tech_sum$AMI),max(metrics_all_tech_sum$AMI)), labels=c("7","1"))+
       scale_x_discrete(position="top", breaks = methods_order, labels = methods_labels)+
       scale_y_discrete(position="left")+
       guides(colour=F)+
       labs(x="", y="", fill="rank")+ 
       theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 8, hjust=0.05, vjust=0.2), axis.text.y = element_text(size=8), legend.position = "right",
       legend.title = element_text(size = 8), legend.text = element_text(size = 8))
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig7", "AMI_avg_heatmap_real_complete", ext="png"), width = 50, height = 80, units = "mm")

# AMI real data version 2
AMI_avg_heatmap_real = metrics_all_tech_sum %>% 
ggplot(aes(x=factor(model, level = methods_order0), y=data_type, fill=AMI))+
       geom_tile()+
       geom_text(aes(label = AMI_rank), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = -1, 
            breaks = c(min(metrics_all_tech_sum$AMI),max(metrics_all_tech_sum$AMI)), labels=c("7","1"))+
       scale_x_discrete(position="top", breaks = methods_order0, labels = methods_labels0)+
       scale_y_discrete(position="left")+
       guides(colour=F)+
       labs(x="", y="", fill="rank")+ 
       theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 9, hjust=0.05, vjust=0.2), axis.ticks.y=element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size=9),
       legend.position = "none")

# VN: real data
VN_avg_heatmap_real = metrics_all_tech_sum %>% 
ggplot(aes(x=factor(model, level = methods_order), y=data_type, fill=VN))+
       geom_tile()+
       geom_text(aes(label = VN_rank), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = -1, 
            breaks = c(min(metrics_all_tech_sum$VN),max(metrics_all_tech_sum$VN)), labels=c("7","1"))+
       scale_x_discrete(position="top", breaks = methods_order, labels = methods_labels)+
       scale_y_discrete(position="left")+
       guides(colour=F)+
       labs(x="", y="", fill="rank")+
       theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 9, hjust=0.05, vjust=0.2), axis.ticks.y=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=9),
       legend.position = "none")
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig7", "VN_avg_heatmap_real", ext="png"), VN_avg_heatmap_real, width = 105, height = 50, units = "mm")

# VN real data version 2
VN_avg_heatmap_real = metrics_all_tech_sum %>% 
ggplot(aes(x=factor(model, level = methods_order0), y=data_type, fill=VN))+
       geom_tile()+
       geom_text(aes(label = VN_rank), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = -1, 
            breaks = c(min(metrics_all_tech_sum$VN),max(metrics_all_tech_sum$VN)), labels=c("7","1"))+
       scale_x_discrete(position="top", breaks = methods_order0, labels = methods_labels0)+
       scale_y_discrete(position="left")+
       guides(colour=F)+
       labs(x="", y="", fill="rank")+ 
       theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.y=element_blank(), 
	axis.ticks.x = element_blank(), axis.text.y = element_text(size=9), legend.position = "none")

# F1: real data
F1_avg_heatmap_real = metrics_all_tech_sum %>% 
ggplot(aes(x=factor(model, level = methods_order), y=data_type, fill=F1))+
       geom_tile()+
       geom_text(aes(label = F1_rank), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = -1, 
            breaks = c(min(metrics_all_tech_sum$F1),max(metrics_all_tech_sum$F1)), labels=c("7","1"))+
       scale_x_discrete(position="top", breaks = methods_order, labels = methods_labels)+
       scale_y_discrete(position="left")+
       guides(colour=F)+
       labs(x="", y="", fill="rank")+
       theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 9, hjust=0.05, vjust=0.2), axis.ticks.y=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=9),
       legend.position = "none")
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig7", "F1_avg_heatmap_real", ext="png"), F1_avg_heatmap_real, width = 105, height = 50, units = "mm")

# F1 real data version 2
F1_avg_heatmap_real = metrics_all_tech_sum %>% 
ggplot(aes(x=factor(model, level = methods_order0), y=data_type, fill=F1))+
       geom_tile()+
       geom_text(aes(label = F1_rank), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = -1, 
            breaks = c(min(metrics_all_tech_sum$F1),max(metrics_all_tech_sum$F1)), labels=c("7","1"))+
       scale_x_discrete(position="top", breaks = methods_order0, labels = methods_labels0)+
       scale_y_discrete(position="left")+
       guides(colour=F)+
       labs(x="", y="", fill="rank")+ 
       theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.y=element_blank(), 
	axis.ticks.x = element_blank(), axis.text.y = element_text(size=9), legend.position = "none")

############### simulation 
spat_prob_order = c(0.9, 0.8, 0.7, 0.6)

# AMI: simulation
AMI_heatmap_sim = sim_metrics_reps_all %>% filter(dims == 5) %>%
ggplot(aes(x=factor(model, level = methods_order), y=factor(spat_prob,spat_prob_order), fill=AMI_mean))+
       geom_tile()+
       geom_text(aes(label = AMI_rank), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = -1)+
       scale_x_discrete(position="top", breaks = methods_order, labels = methods_labels)+
       scale_y_discrete(labels = c("0.6" = "p = 0.6", "0.7" = "p = 0.7", "0.8" = "p = 0.8", "0.9" = "p = 0.9"), position = "left")+
       guides(colour=F)+
       labs(x="", y="", fill="rank")+
       theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 9, hjust=0.05, vjust=0.2), axis.ticks.y=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size = 9),
       legend.position = "none")
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig7", "AMI_heatmap_sim", ext="png"), AMI_heatmap_sim, width = 100, height = 56, units = "mm")

# VN: simulation
VN_heatmap_sim = sim_metrics_reps_all %>% filter(dims == 5) %>%
ggplot(aes(x=factor(model, level = methods_order), y=factor(spat_prob,spat_prob_order), fill=VN_mean))+
       geom_tile()+
       geom_text(aes(label = VN_rank), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = -1)+
       scale_x_discrete(position="top", breaks = methods_order, labels = methods_labels)+
       scale_y_discrete(labels = c("0.6" = "p = 0.6", "0.7" = "p = 0.7", "0.8" = "p = 0.8", "0.9" = "p = 0.9"), position = "left")+
       guides(colour=F)+
       labs(x="", y="", fill="VN rank")+
       theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 9, hjust=0.05, vjust=0.2), axis.ticks.y=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=9),
       legend.position = "none")
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig7", "VN_heatmap_sim", ext="png"), VN_heatmap_sim, width = 100, height = 56, units = "mm")

# F1: simulation
F1_heatmap_sim = sim_metrics_reps_all %>% filter(dims == 5) %>%
ggplot(aes(x=factor(model, level = methods_order), y=factor(spat_prob,spat_prob_order), fill=F1_mean))+
       geom_tile()+
       geom_text(aes(label = F1_rank), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = -1)+
       scale_x_discrete(position="top", breaks = methods_order, labels = methods_labels)+
       scale_y_discrete(labels = c("0.6" = "p = 0.6", "0.7" = "p = 0.7", "0.8" = "p = 0.8", "0.9" = "p = 0.9"), position="left")+
       guides(colour=F)+
       labs(x="", y="", fill="F1 rank")+
       theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 9, hjust=0.05, vjust=0.2), axis.ticks.y=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=9),
       legend.position = "none")
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig7", "F1_heatmap_sim", ext="png"), F1_heatmap_sim, width = 100, height = 56, units = "mm")

################# runtime and memory
# runtime_cells
runtime_cells_heatmap = runtime_cells_df %>% mutate(method = factor(method, levels = methods_order)) %>%
ggplot(aes(x = method, y = factor(num_cells, levels = c(4500, 3000, 1500, 1000, 500)), fill = runtime)) +
	geom_tile() +
	geom_text(aes(label = runtime_rank), size = 4)+
	scale_fill_viridis(discrete = F, option = "C", begin=0.25, direction=1)+
	scale_x_discrete(position = "top")+
	scale_y_discrete(labels = c("500"="500 cells", "1000"="1000 cells", "1500"="1500 cells", "3000"="3000 cells", "4500"="4500 cells"), position="left")+
	guides(colour=F)+
	labs(x="", y="", fill = "")+
	theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 9, hjust=0.05, vjust=0.2), axis.ticks.y=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=9),
       legend.position = "none")
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig7", "runtime_cells_heatmap", ext="png"), runtime_cells_heatmap, width = 105, height = 62, units = "mm")

# memory_cells
memory_cells_heatmap = memory_cells_df %>% mutate(method = factor(method, levels = methods_order)) %>%
ggplot(aes(x = method, y = factor(num_cells, levels = c(4500, 3000, 1500, 1000, 500)), fill = memory)) +
	geom_tile() +
	geom_text(aes(label = memory_rank), size = 4)+
	scale_fill_viridis(discrete = F, option = "C", begin=0.25, direction=1)+
	scale_x_discrete(position = "top")+
	scale_y_discrete(labels = c("500"="500 cells", "1000"="1000 cells", "1500"="1500 cells", "3000"="3000 cells", "4500"="4500 cells"), position="left")+
	guides(colour=F)+
	labs(x="", y="", fill = "")+
	theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 9, hjust=0.05, vjust=0.2), axis.ticks.y=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=9),
       legend.position = "none")
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig7", "memory_cells_heatmap", ext="png"), memory_cells_heatmap, width = 105, height = 62, units = "mm")

# runtime features
runtime_features_heatmap = runtime_features_df %>% mutate(method = factor(method, levels = methods_order)) %>%
ggplot(aes(x=method, y = factor(num_features, levels = c(10000, 5000, 3000, 1000, 500)), fill=runtime))+
       geom_tile()+
       geom_text(aes(label = runtime_rank), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = 1)+
       scale_x_discrete(position="top")+
       scale_y_discrete(labels = c("500"="500 features", "1000"="1000 features", "3000"="3000 features", "5000"="5000 features", "10000"="10000 features"), position="left")+
       guides(colour=F)+
       labs(x="", y="", fill="")+
       theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 9, hjust=0.05, vjust=0.2), axis.ticks.y=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=9),
       legend.position = "none")
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig7", "runtime_features_heatmap", ext="png"), runtime_features_heatmap, width = 110, height = 62, units = "mm")

# memory_features
memory_features_heatmap = memory_features_df %>% mutate(method = factor(method, levels = methods_order)) %>%
  ggplot(aes(x = method, y = factor(num_features, levels = c(10000, 5000, 3000,	1000, 500)), fill = memory)) +
  geom_tile() +
  geom_text(aes(label = memory_rank), size = 4)+
  scale_fill_viridis(discrete = F, option = "C", begin=0.25, direction=1)+
  scale_x_discrete(position = "top")+
  scale_y_discrete(labels = c("500"="500 features", "1000"="1000 features", "3000"="3000 features", "5000"="5000 features", "10000"="10000 features"), position="left")+
  guides(colour=F)+
  labs(x="", y="", fill = "")+
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 9, hjust=0.05, vjust=0.2), axis.ticks.y=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=9),
       legend.position = "none")
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig7", "memory_features_heatmap", ext="png"), memory_features_heatmap, width = 110, height = 62, units = "mm")

############ Figure 7; 08/24/2021
### real data
metrics_all_gather_values = metrics_all %>% select(data_name, model, AMI, VN, F1) %>% gather("metric","value", -data_name, -model)
metrics_all_gather_ranking = metrics_all %>% select(data_name, model, AMI_rank, VN_rank, F1_rank) %>% gather("metric","ranking",-data_name, -model) 
metrics_all_gather = data.frame(data_name = metrics_all_gather_values$data_name, 
			 model = metrics_all_gather_values$model,
			 metric = paste(metrics_all_gather_values$metric, "real", sep="_"),
			 value = metrics_all_gather_values$value,
			 ranking = metrics_all_gather_ranking$ranking)
metrics_all_sum = metrics_all %>% select()

### simulation datasets
sim_metrics_reps_all_gather_values = sim_metrics_reps_all %>% filter(dims == 5) %>% select(spat_prob, model, AMI_mean, VN_mean, F1_mean) %>% gather("metric", "value", -spat_prob, -model)
sim_metrics_reps_all_gather_ranking = sim_metrics_reps_all %>% filter(dims == 5) %>% select(spat_prob, model, AMI_rank, VN_rank, F1_rank) %>% gather("metric", "ranking", -spat_prob, -model)
sim_metrics_reps_all_gather = data.frame(data_name = paste("p =", sim_metrics_reps_all_gather_values$spat_prob, sep=" "),
				model = sim_metrics_reps_all_gather_values$model,
				metric = paste(c(rep("AMI", 28), rep("VN",28), rep("F1", 28)), "sim", sep="_"),
				value = sim_metrics_reps_all_gather_values$value,
				ranking = sim_metrics_reps_all_gather_ranking$ranking)

sim_metrics_reps_all_sum = sim_metrics_reps_all %>% filter(dims == 5) %>% dplyr::group_by(model) %>% 
	dplyr::summarise(AMI = mean(AMI_mean), VN = mean(VN_mean), F1 = mean(F1_mean)) %>% 
	dplyr::mutate(AMI_rank = rank(-1.0*AMI, ties.method = "first"), VN_rank = rank(-1.0*VN, ties.method = "first"),
		F1_rank = rank(-1.0*F1, ties.method = "first"))
sim_metrics_reps_all_sum$data_name = "simulation"
sim_metrics_reps_all_sum_gather_values = sim_metrics_reps_all_sum %>% select(data_name, model, AMI, VN, F1) %>% gather("metric", "value", -data_name, -model)
sim_metrics_reps_all_sum_gather_ranking = sim_metrics_reps_all_sum %>% select(data_name, model, AMI_rank, VN_rank, F1_rank) %>% gather("metric", "ranking", -data_name, -model)
sim_metrics_reps_all_sum_gather = data.frame(data_name = sim_metrics_reps_all_sum_gather_values$data_name,
	model = sim_metrics_reps_all_sum_gather_values$model,
	metric = paste(c(rep("AMI",7), rep("VN",7), rep("F1",7)), "sim", sep="_"),
	value = sim_metrics_reps_all_sum_gather_values$value,
	ranking = sim_metrics_reps_all_sum_gather_ranking$ranking)

runtime_cells_df_copy = runtime_cells_df
memory_cells_df_copy = memory_cells_df
runtime_features_df_copy = runtime_features_df
memory_features_df_copy = memory_features_df

metrics_level_2_attach_cells = data.frame(data_name = c(500, 1000, 1500, 3000, 4500),
			          model = rep("SG", 5),
			          value = NA,
			          ranking = NA)

colnames(runtime_cells_df_copy) = c("data_name", "model", "value", "ranking")
runtime_cells_df_copy = rbind(runtime_cells_df_copy, metrics_level_2_attach_cells)
runtime_cells_df_copy$data_name = paste(runtime_cells_df_copy$data_name, "cells", sep=" ")
runtime_cells_df_copy$metric = rep("runtime_cells", nrow(runtime_cells_df_copy))

colnames(memory_cells_df_copy) = c("data_name", "model", "value", "ranking")
memory_cells_df_copy = rbind(memory_cells_df_copy, metrics_level_2_attach_cells)
memory_cells_df_copy$data_name = paste(memory_cells_df_copy$data_name, "cells", sep=" ")
memory_cells_df_copy$metric = rep("memory_cells", nrow(memory_cells_df_copy))

metrics_level_2_attach_features = data.frame(data_name = c(500, 1000, 3000, 5000, 10000),
			          model = rep("SG", 5),
			          value = NA,
			          ranking = NA)

colnames(runtime_features_df_copy) = c("data_name", "model", "value", "ranking")
runtime_features_df_copy = rbind(runtime_features_df_copy, metrics_level_2_attach_features)
runtime_features_df_copy$data_name = paste(runtime_features_df_copy$data_name, "features", sep=" ")
runtime_features_df_copy$metric = rep("runtime_features", nrow(runtime_features_df_copy))

colnames(memory_features_df_copy) = c("data_name", "model", "value", "ranking")
memory_features_df_copy = rbind(memory_features_df_copy, metrics_level_2_attach_features)
memory_features_df_copy$data_name = paste(memory_features_df_copy$data_name, "features", sep=" ")
memory_features_df_copy$metric = rep("memory_features", nrow(memory_features_df_copy))

metrics_level_2_attach_features_sum = data.frame(data_name = c(500, 1000, 3000, 5000, 10000),
			          model = rep("SG", 5),
			          value = NA,
			          ranking = NA)

runtime_cells_df_sum = runtime_cells_df %>% select(-runtime_rank) %>% 
	dplyr::rename(data_name = num_cells, model = method, value = runtime) %>% 
	dplyr::group_by(model) %>% dplyr::summarise(value = mean(value)) %>%
	dplyr::mutate(ranking = rank(value, ties.method = "first"), 
	data_name = "runtime (cells)") %>%
	add_row(data_name = "runtime (cells)", model = "SG", value = NA, ranking = NA) %>%
	dplyr::mutate(metric = "runtime_cells") %>%
	select(data_name, model, metric, value, ranking) 

memory_cells_df_sum = memory_cells_df %>% select(-memory_rank) %>% 
	dplyr::rename(data_name = num_cells, model = method, value = memory) %>% 
	dplyr::group_by(model) %>% dplyr::summarise(value = mean(value)) %>%
	dplyr::mutate(ranking = rank(value, ties.method = "first"),
	data_name = "memory (cells)") %>%
	add_row(data_name = "memory (cells)", model = "SG", value = NA, ranking = NA) %>%
	dplyr::mutate(metric = "memory_cells") %>%
	select(data_name, model, metric, value, ranking)

runtime_features_df_sum = runtime_features_df %>% select(-runtime_rank) %>% 
	dplyr::rename(data_name = num_features, model = method, value = runtime) %>% 
	dplyr::group_by(model) %>% dplyr::summarise(value = mean(value)) %>%
	dplyr::mutate(ranking = rank(value, ties.method = "first"),
	data_name = "runtime (features)") %>%
	add_row(data_name = "runtime (features)", model = "SG", value = NA, ranking = NA) %>%
	dplyr::mutate(metric = "runtime_features") %>%
	select(data_name, model, metric, value, ranking)

memory_features_df_sum = memory_features_df %>% select(-memory_rank) %>% 
	dplyr::rename(data_name = num_features, model = method, value = memory) %>% 
	dplyr::group_by(model) %>% dplyr::summarise(value = mean(value)) %>%
	dplyr::mutate(ranking = rank(value, ties.method = "first"),
	data_name = "memory (features)") %>%
	add_row(data_name = "memory (features)", model = "SG", value = NA, ranking = NA) %>%
	dplyr::mutate(metric = "memory_features") %>%
	select(data_name, model, metric, value, ranking)

metrics_fig7 = full_join(full_join(full_join(full_join(full_join(metrics_all_gather, sim_metrics_reps_all_gather),
	          runtime_cells_df_copy), 
	          memory_cells_df_copy),
    	          runtime_features_df_copy), 
	          memory_features_df_copy)
save(metrics_fig7, file = fs::path("data_analysis/real_data/vis_manuscript", "metrics_fig7", ext="RData"))

metrics_fig7_v2 = full_join(full_join(full_join(full_join(full_join(metrics_all_gather, sim_metrics_reps_all_sum_gather),
	          runtime_cells_df_sum), 
	          memory_cells_df_sum),
    	          runtime_features_df_sum), 
	          memory_features_df_sum)
save(metrics_fig7_v2, file = fs::path("data_analysis/real_data/vis_manuscript", "metrics_fig7_v2", ext="RData"))

fig7_metrics = c("AMI_real", "VN_real", "F1_real", "AMI_sim", "VN_sim", "F1_sim", 
	           "runtime_cells", "memory_cells", "runtime_features", "memory_features")
fig7_y_order = list(c("Visium_MK_cor", "Visium_MBSP1", "Visium_MBSA1", "Visium_MB_cor", "Visium_HCere", "ST_MOB1", 
	"sfp_SScortex", "mf_hypo_-290", "mf_hypo_-240", "mf_hypo_-190", "mf_hypo_-140", "mf_hypo_-90", 
	"mf_hypo_-40", "mf_hypo_10", "mf_hypo_60", "mf_hypo_110", "mf_hypo_160", "mf_hypo_210", 
	"mf_hypo_260"), 
	c("Visium_MK_cor", "Visium_MBSP1", "Visium_MBSA1", "Visium_MB_cor", "Visium_HCere", "ST_MOB1", 
	"sfp_SScortex", "mf_hypo_-290", "mf_hypo_-240", "mf_hypo_-190", "mf_hypo_-140", "mf_hypo_-90", 
	"mf_hypo_-40", "mf_hypo_10", "mf_hypo_60", "mf_hypo_110", "mf_hypo_160", "mf_hypo_210", 
	"mf_hypo_260"),
	c("Visium_MK_cor", "Visium_MBSP1", "Visium_MBSA1", "Visium_MB_cor", "Visium_HCere", "ST_MOB1", 
	"sfp_SScortex", "mf_hypo_-290", "mf_hypo_-240", "mf_hypo_-190", "mf_hypo_-140", "mf_hypo_-90", 
	"mf_hypo_-40", "mf_hypo_10", "mf_hypo_60", "mf_hypo_110", "mf_hypo_160", "mf_hypo_210", 
	"mf_hypo_260"),
	c("p = 0.6", "p = 0.7", "p = 0.8", "p = 0.9"), 
	c("p = 0.6", "p = 0.7", "p = 0.8", "p = 0.9"),
	c("p = 0.6", "p = 0.7", "p = 0.8", "p = 0.9"),
	c("500 cells", "1000 cells", "1500 cells", "3000 cells", "4500 cells"), 
	c("500 cells", "1000 cells", "1500 cells", "3000 cells", "4500 cells"),
	c("500 features", "1000 features", "3000 features", "5000 features", "10000 features"),
	c("500 features", "1000 features", "3000 features", "5000 features", "10000 features"))

fig7_y_order_sum = list(c("Visium_MK_cor", "Visium_MBSP1", "Visium_MBSA1", "Visium_MB_cor", "Visium_HCere", "ST_MOB1", 
	"sfp_SScortex", "mf_hypo_-290", "mf_hypo_-240", "mf_hypo_-190", "mf_hypo_-140", "mf_hypo_-90", 
	"mf_hypo_-40", "mf_hypo_10", "mf_hypo_60", "mf_hypo_110", "mf_hypo_160", "mf_hypo_210", 
	"mf_hypo_260"), 
	c("Visium_MK_cor", "Visium_MBSP1", "Visium_MBSA1", "Visium_MB_cor", "Visium_HCere", "ST_MOB1", 
	"sfp_SScortex", "mf_hypo_-290", "mf_hypo_-240", "mf_hypo_-190", "mf_hypo_-140", "mf_hypo_-90", 
	"mf_hypo_-40", "mf_hypo_10", "mf_hypo_60", "mf_hypo_110", "mf_hypo_160", "mf_hypo_210", 
	"mf_hypo_260"),
	c("Visium_MK_cor", "Visium_MBSP1", "Visium_MBSA1", "Visium_MB_cor", "Visium_HCere", "ST_MOB1", 
	"sfp_SScortex", "mf_hypo_-290", "mf_hypo_-240", "mf_hypo_-190", "mf_hypo_-140", "mf_hypo_-90", 
	"mf_hypo_-40", "mf_hypo_10", "mf_hypo_60", "mf_hypo_110", "mf_hypo_160", "mf_hypo_210", 
	"mf_hypo_260"),
	c("simulation"), 
	c("simulation"), 
	c("simulation"), 
	c("runtime (cells)"), 
	c("memory (cells)"),
	c("runtime (features)"),
	c("memory (features)"))

fig7_list = list()
metrics_fig7_temp_list = list()

methods_order0 = c("CIMLR", "concatenation", "MOFA+", "scVI", "SNF", "WNN", "SG")
methods_labels0 = c("CIMLR", "concatenation", "MOFA+", "scVI", "SNF", "Seurat v4", "SV genes")

for(i in 4:length(fig7_metrics)){ #length(fig7_metrics)
  metrics_fig7_temp_list[[i]] = metrics_fig7_v2 %>% filter(metric == fig7_metrics[i]) %>% 
	mutate(data_name = factor(data_name, levels = rev(fig7_y_order_sum[[i]]))) %>%
	mutate(model = factor(model, levels = methods_order0))
  #y_order = rev(unique(metrics_fig7_temp$data_name))

  if(i == 1){
    fig7_list[[i]] = metrics_fig7_temp_list[[i]] %>%
       ggplot(aes(x=model, y= data_name, fill=value))+
       geom_tile()+
       geom_text(aes(label = ranking), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = -1, 
            breaks = c(min(metrics_fig7_temp_list[[i]]$value),max(metrics_fig7_temp_list[[i]]$value)), labels=c("7","1"))+
       scale_x_discrete(position="top", breaks = methods_order0, labels = methods_labels0)+
       scale_y_discrete(position="left")+
       guides(colour=F)+
       labs(x="", y="", fill="rank")+
       theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 9, hjust=0.05, vjust=0.2), axis.ticks.y=element_blank(), 
	axis.ticks.x=element_blank(), axis.text.y = element_text(size=9), legend.position = "none",
	plot.margin = unit(c(0, 0, 0, 0),"cm"))
  }else if(i > 6){
    fig7_list[[i]] = metrics_fig7_temp_list[[i]] %>%
       ggplot(aes(x=model, y=data_name, fill=value))+
       geom_tile()+
       geom_text(aes(label = ranking), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = 1, 
            breaks = c(min(metrics_fig7_temp_list[[i]]$value),max(metrics_fig7_temp_list[[i]]$value)), labels=c("7","1"))+
       scale_x_discrete(position="top", breaks = methods_order0, labels = methods_labels0)+
       scale_y_discrete(position="left")+
       guides(colour=F)+
       labs(x="", y="", fill="rank")+
       theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.y=element_blank(), 
	axis.ticks.x=element_blank(), axis.text.y = element_text(size=9), legend.position = "none",
	plot.margin = unit(c(0, 0, 0, 0),"cm"))
  }else{
    fig7_list[[i]] = metrics_fig7_temp_list[[i]] %>%
       ggplot(aes(x=model, y=data_name, fill=value))+
       geom_tile()+
       geom_text(aes(label = ranking), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = -1, 
            breaks = c(min(metrics_fig7_temp_list[[i]]$value),max(metrics_fig7_temp_list[[i]]$value)), labels=c("7","1"))+
       scale_x_discrete(position="top", breaks = methods_order0, labels = methods_labels0)+
       scale_y_discrete(position="left")+
       guides(colour=F)+
       labs(x="", y="", fill="rank")+
       theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.y=element_blank(), 
	axis.ticks.x=element_blank(), axis.text.y = element_text(size=9), legend.position = "none",
	plot.margin = unit(c(0, 0, 0, 0),"cm"))
  }
}
#unit: top, right, bottom, left
#fig7 = plot_grid(fig7_list[[1]], fig7_list[[2]], fig7_list[[3]], fig7_list[[4]],
#	         align="v", ncol=1, rel_heights=c(3, 2.3, 2.3, 0.7))
#fig7 = plot_grid(fig7_list[[1]], fig7_list[[2]], fig7_list[[3]], fig7_list[[4]], fig7_list[[5]],
#	         fig7_list[[6]], fig7_list[[7]], fig7_list[[8]], fig7_list[[9]], fig7_list[[10]],
#	         align="v", ncol=1, rel_heights=c(7, 5, 5, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65))
#ggsave(fig7, file = fs::path("data_analysis/real_data/vis_manuscript/Fig7", "Fig7_0824", ext="png"), width = 110, height = 297, units = "mm")
fig7 = plot_grid(AMI_avg_heatmap_real, VN_avg_heatmap_real, F1_avg_heatmap_real, 
	fig7_list[[4]], fig7_list[[5]], fig7_list[[6]], fig7_list[[7]], fig7_list[[8]], fig7_list[[9]], fig7_list[[10]],
	align = "v", ncol = 1, rel_heights = c(5, 3, 3, 1.3, 1.3, 1.3 ,1.3, 1.3, 1.3, 1.3))
ggsave(fig7, file = fs::path("data_analysis/real_data/vis_manuscript/Fig7", "Fig7_0826", ext="png"), width = 110, height = 180, units = "mm")

################ supplementary fig 4
supp_fig4 = plot_grid(fig7_list[[4]], fig7_list[[5]],fig7_list[[6]], align = "h", nrow = 1)
ggsave(supp_fig4, 
	file = fs::path("data_analysis/real_data/vis_manuscript/Fig2_supp", "fig_4_supp_0825", ext="png"), width = 255, height = 60, units = "mm")

################ Figure 7: spider plots
metrics_all_radar_sum = metrics_all %>% filter(model != "BREM-SC") %>% group_by(model) %>% dplyr::summarise(ARI_real = mean(ARI), AMI_real = mean(AMI), DBI_real = mean(DBI), 
	DBI_truth_real = mean(DBI_truth), AIC_real = mean(AIC), AIC0_real = mean(AIC0), VN_real = mean(VN), TPR_real = mean(TPR), PPV_real = mean(PPV), FDR_real = mean(FDR), 
	FM_real = mean(FM), CSI_real = mean(CSI), ACC_real = mean(ACC), F1_real = mean(F1))

sim_metrics_radar_sum = sim_metrics_reps_all %>% group_by(model) %>% dplyr::summarise(ARI_sim = mean(ARI_mean), AMI_sim = mean(AMI_mean), DBI_sim = mean(DBI_mean), 
	DBI_truth_sim = mean(DBI_truth_mean), AIC_sim = mean(AIC_mean), AIC0_sim = mean(AIC0_mean), VN_sim = mean(VN_mean), TPR_sim = mean(TPR_mean), PPV_sim= mean(PPV_mean), 
	FDR_sim = mean(FDR_mean), FM_sim = mean(FM_mean), CSI_sim = mean(CSI_mean), ACC_sim = mean(ACC_mean), F1_sim = mean(F1_mean))
radar_combine_df = merge(metrics_all_radar_sum, sim_metrics_radar_sum, by="model")

## average real and sim
radar_combine_df$AMI = (radar_combine_df$AMI_real*19 + radar_combine_df$AMI_sim*20)/39
radar_combine_df$VN = (radar_combine_df$VN_real*19 + radar_combine_df$VN_sim*20)/39
radar_combine_df$F1 = (radar_combine_df$F1_real*19 + radar_combine_df$F1_sim*20)/39

## load cell runtime and memory
runtime_cells_df_sum = runtime_cells_df %>% group_by(method) %>% dplyr::summarise(runtime_cells = -1.0*mean(runtime))
radar_combine_df = merge(radar_combine_df, runtime_cells_df_sum, by.x="model", by.y="method")

memory_cells_df_sum = memory_cells_df %>% group_by(method) %>% dplyr::summarise(memory_cells = -1.0*mean(memory))
radar_combine_df = merge(radar_combine_df, memory_cells_df_sum, by.x="model", by.y="method")

#### load features runtime and memory
runtime_features_df_sum = runtime_features_df %>% group_by(method) %>% dplyr::summarise(runtime_features = -1.0*mean(runtime))
radar_combine_df = merge(radar_combine_df, runtime_features_df_sum, by.x="model", by.y="method")

memory_features_df_sum = memory_features_df %>% group_by(method) %>% dplyr::summarise(memory_features = -1.0*mean(memory))
radar_combine_df = merge(radar_combine_df, memory_features_df_sum, by.x="model", by.y="method")

radar_combine_df %>% select(model, AMI, VN, F1, runtime_cells, memory_cells, runtime_features, memory_features)

save(radar_combine_df, file = "/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/radar_combine_df.RData")

library(fmsb)
radar_max = c(rep(1,3), 0, 0, 0, 0)
radar_min = c(rep(0,3), -93, -11383, -27, -1277)

for(i in 1:length(methods_order)){
  radar_tmp_df = rbind(radar_max, radar_min)
  #colnames(radar_tmp_df) = c("AMI_real", "VN_real", "F1_real", "AMI_sim", "VN_sim", "F1_sim", "runtime_cells", "memory_cells", "runtime_features", "memory_features")
  colnames(radar_tmp_df) = c("AMI", "VN", "F1", "runtime_cells", "memory_cells", "runtime_features", "memory_features")
  #radar_tmp_row = radar_combine_df %>% filter(model == methods_order[i]) %>% select(AMI_real, VN_real, F1_real, AMI_sim, VN_sim, F1_sim, runtime_cells, memory_cells, runtime_features, memory_features)
  radar_tmp_row = radar_combine_df %>% filter(model == methods_order[i]) %>% select(AMI, VN, F1, runtime_cells, memory_cells, runtime_features, memory_features)
  radar_tmp_df = rbind(radar_tmp_df, radar_tmp_row[1,])  
  #colnames(radar_tmp_df) = c("AMI (real data)", "VN (real data)", "F1 (real data)", "AMI (simulation)", "VN (simulation)", "F1 (simulation)", "runtime (cells)", "memory (cells)", "runtime (features)", "memory (features)")
  colnames(radar_tmp_df) = c("AMI", "VN", "F1", "runtime (cells)", "memory (cells)", "runtime (features)", "memory (features)")

  png(file=fs::path("data_analysis/real_data/vis_manuscript/Fig7", paste(methods_labels[i], "spiderplot", sep="_"), ext="png"), width=6, height=6, units="in", res=1200)
  col = c(col2rgb(manual_col_pallete0[methods_order[i]]))/255
  radarchart(radar_tmp_df, axistype=0,  pcol=rgb(red = col[1], green = col[2], blue = col[3], alpha = 0.9) , pfcol=rgb(red = col[1], green = col[2], blue = col[3], alpha=0.5) , plwd=4 , cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8, vlcex=0.1)
  dev.off()
}



############################### supplementary figures ################
if(!dir.exists("data_analysis/real_data/vis_manuscript/Fig2_supp")){
  dir.create("data_analysis/real_data/vis_manuscript/Fig2_supp")
}

##################### barplots 
#AMI supp
plot_supp_list = list()
for(i in 1:length(supp_data)){
  model_levels_supp_AMI = metrics_all_supp %>% filter(model != "BREM-SC") %>%filter(data_name %in% supp_data[i]) %>% arrange(desc(AMI)) %>% pull(model)
  SG_idx = which(model_levels_supp_AMI == "SG")
  model_levels_supp_AMI = c("SG", model_levels_supp_AMI[-SG_idx])
  plot_supp_list[[i]] = metrics_all_supp %>% filter(model != "BREM-SC") %>%filter(data_name %in% supp_data[i]) %>% mutate(model = factor(model, levels = model_levels_supp_AMI)) %>%
    ggplot(aes(x=data_name, y=AMI, fill=model)) +
    geom_bar(position="dodge", stat="identity") +
    scale_fill_manual(values = manual_col_pallete0, breaks = methods_order, labels = methods_labels) +
    labs(x="", y="", title="")+ ylim(0,1)+
    theme_bw() + theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.position="none", axis.text.x = element_text(size=5.5),
    plot.margin = unit(c(0, 0, 0, 0),"cm")) 
}
plot_supp_list[[1]] = plot_supp_list[[1]]+theme(axis.ticks.y=element_line(), axis.text.y=element_text(size=7))
AMI_real_data_supp = grid.arrange(plot_supp_list[[1]], plot_supp_list[[2]], plot_supp_list[[3]], plot_supp_list[[4]], plot_supp_list[[5]], plot_supp_list[[6]], plot_supp_list[[7]], plot_supp_list[[8]], plot_supp_list[[9]], widths = c(1.2, rep(1,8)), ncol=9)
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig2_supp", "AMI_real_data_supp", ext="png"), AMI_real_data_supp, width = 212, height = 64, units = "mm")

#VN supp
plot_supp_list = list()
for(i in 1:length(supp_data)){
  model_levels_supp_VN = metrics_all_supp %>% filter(model != "BREM-SC") %>% filter(data_name %in% supp_data[i]) %>% arrange(desc(VN)) %>% pull(model)
  SG_idx = which(model_levels_supp_VN == "SG")
  model_levels_supp_VN = c("SG", model_levels_supp_VN[-SG_idx])
  plot_supp_list[[i]] = metrics_all_supp %>% filter(model != "BREM-SC") %>% filter(data_name %in% supp_data[i]) %>% mutate(model = factor(model, levels = model_levels_supp_VN)) %>%
    ggplot(aes(x=data_name, y=VN, fill=model)) +
    geom_bar(position="dodge", stat="identity") +
    scale_fill_manual(values = manual_col_pallete0, breaks = methods_order, labels = methods_labels) +
    labs(x="", y="", title="")+ ylim(0,1)+
    theme_bw() + theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.position="none", axis.text.x = element_text(size=5.5),
    plot.margin = unit(c(0, 0, 0, 0),"cm")) 
}
plot_supp_list[[1]] = plot_supp_list[[1]]+theme(axis.ticks.y=element_line(), axis.text.y=element_text(size=7))
VN_real_data_supp = grid.arrange(plot_supp_list[[1]], plot_supp_list[[2]], plot_supp_list[[3]], plot_supp_list[[4]], plot_supp_list[[5]], plot_supp_list[[6]], plot_supp_list[[7]], plot_supp_list[[8]], plot_supp_list[[9]], widths = c(1.2, rep(1,8)), ncol=9)
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig2_supp", "VN_real_data_supp", ext="png"), VN_real_data_supp, width = 212, height = 64, units = "mm")

#F1 supp
plot_supp_list = list()
for(i in 1:length(supp_data)){
  model_levels_supp_F1 = metrics_all_supp %>% filter(model != "BREM-SC") %>% filter(data_name %in% supp_data[i]) %>% arrange(desc(F1)) %>% pull(model)
  SG_idx = which(model_levels_supp_F1 == "SG")
  model_levels_supp_F1 = c("SG", model_levels_supp_F1[-SG_idx])
  plot_supp_list[[i]] = metrics_all_supp %>% filter(model != "BREM-SC") %>% filter(data_name %in% supp_data[i]) %>% mutate(model = factor(model, levels = model_levels_supp_F1)) %>%
    ggplot(aes(x=data_name, y=F1, fill=model)) +
    geom_bar(position="dodge", stat="identity") +
    scale_fill_manual(values = manual_col_pallete0, breaks = methods_order, labels = methods_labels) +
    labs(x="", y="", title="")+ ylim(0,1)+
    theme_bw() + theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.position="none", axis.text.x = element_text(size=5.5),
    plot.margin = unit(c(0, 0, 0, 0),"cm")) 
}
plot_supp_list[[1]] = plot_supp_list[[1]]+theme(axis.ticks.y=element_line(), axis.text.y=element_text(size=7))
F1_real_data_supp = grid.arrange(plot_supp_list[[1]], plot_supp_list[[2]], plot_supp_list[[3]], plot_supp_list[[4]], plot_supp_list[[5]], plot_supp_list[[6]], plot_supp_list[[7]], plot_supp_list[[8]], plot_supp_list[[9]], widths = c(1.2, rep(1,8)), ncol=9)
ggsave(file = fs::path("data_analysis/real_data/vis_manuscript/Fig2_supp", "F1_real_data_supp", ext="png"), F1_real_data_supp, width = 212, height = 64, units = "mm")

############################# heatmaps
# AMI supp
metrics_all_supp %>% filter(model != "BREM-SC") %>%
ggplot(aes(x=factor(model, level = methods_order), y=data_name, fill=AMI))+
       geom_tile()+
       geom_text(aes(label = AMI_rank), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = -1, 
            breaks = c(min(metrics_all_supp$AMI),max(metrics_all_supp$AMI)), labels=c("7","1"))+
       scale_x_discrete(position="top", breaks = methods_order, labels = methods_labels)+
       guides(colour=F)+
       labs(x="", y="", fill="rank")+
       theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 8, hjust=0.05, vjust=0.2), axis.text.y = element_text(size=8), 
        legend.title = element_text(size=8), legend.text = element_text(size=8))
ggsave(fs::path("data_analysis/real_data/vis_manuscript/Fig2_supp", "AMI_heatmap_supp", ext="png"), width = 85, height = 76, units = "mm") #change this

# VN supp
metrics_all_supp %>% filter(model != "BREM-SC") %>%
ggplot(aes(x=factor(model, level = methods_order), y=data_name, fill=VN))+
       geom_tile()+
       geom_text(aes(label = VN_rank), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = -1, 
            breaks = c(min(metrics_all_supp$VN),max(metrics_all_supp$VN)), labels=c("7","1"))+
       scale_x_discrete(position="top", breaks = methods_order, labels = methods_labels)+
       guides(colour=F)+
       labs(x="", y="", fill="rank")+
       theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 8, hjust=0.05, vjust=0.2), axis.text.y = element_text(size=8), 
        legend.title = element_text(size=8), legend.text = element_text(size=8))
ggsave(fs::path("data_analysis/real_data/vis_manuscript/Fig2_supp", "VN_heatmap_supp", ext="png"), width = 85, height = 76, units = "mm") #change this

# F1 supp
metrics_all_supp %>% filter(model != "BREM-SC") %>%
ggplot(aes(x=factor(model, level = methods_order), y=data_name, fill=F1))+
       geom_tile()+
       geom_text(aes(label = F1_rank), size = 4)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = -1, 
            breaks = c(min(metrics_all_supp$F1),max(metrics_all_supp$F1)), labels=c("7","1"))+
       scale_x_discrete(position="top", breaks = methods_order, labels = methods_labels)+
       guides(colour=F)+
       labs(x="", y="", fill="rank")+
       theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 8, hjust=0.05, vjust=0.2), axis.text.y = element_text(size=8), 
        legend.title = element_text(size=8), legend.text = element_text(size=8))
ggsave(fs::path("data_analysis/real_data/vis_manuscript/Fig2_supp", "F1_heatmap_supp", ext="png"), width = 85, height = 76, units = "mm") #change this

