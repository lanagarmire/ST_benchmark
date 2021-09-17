#rm(list=ls())
setwd("/home/liyijun/ST_benchmark_01082020")
source("functions/spat_metric_03192021.R")
library(aricode)
library(clusterSim)
library(data.table)
library(dplyr)
library(Giotto)

#### get arguments from bash scripts
args = commandArgs(trailingOnly = TRUE)
data_path = as.character(args[1])
thres = as.numeric(args[2])
data_ref = as.character(args[3])
diff_in_cov = as.numeric(args[4])
svg_pval_thres = as.numeric(args[5])
num_nn = as.numeric(args[6])
ndims = as.numeric(args[7])
scVI_error = as.character(args[8])
SNF_metric = as.character(args[9])
WNN_smooth = as.numeric(args[10])
CIMLR_k = as.numeric(args[11])
save_name = as.character(args[12])
n_random_reps = as.numeric(args[13])

### set parameters
if(save_name == "thres"){
  save_path = fs::path(data_path, data_ref, paste(save_name, format(thres,drop0Trailing=F), sep="_"))
}else if(save_name == "diff_in_cov"){
  save_path = fs::path(data_path, data_ref, paste(save_name, format(diff_in_cov,drop0Trailing=F), sep="_"))
}else if(save_name == "svg_pval_thres"){
  save_path = fs::path(data_path, data_ref, paste(save_name, format(svg_pval_thres,drop0Trailing=F), sep="_"))
}else if(save_name == "ndims"){
  save_path = fs::path(data_path, data_ref, paste(save_name, format(ndims,drop0Trailing=F), sep="_"))
}else if(save_name == "scVI_error"){
  save_path = fs::path(data_path, data_ref, paste(save_name, format(scVI_error,drop0Trailing=F), sep="_"))
}else if(save_name == "SNF_metric"){
  save_path = fs::path(data_path, data_ref, paste(save_name, format(SNF_metric,drop0Trailing=F), sep="_"))
}else if(save_name == "WNN_smooth"){
  save_path = fs::path(data_path, data_ref, paste(save_name, format(WNN_smooth,drop0Trailing=F), sep="_"))
}else if(save_name == "CIMLR_k"){
  save_path = fs::path(data_path, data_ref, paste(save_name, format(CIMLR_k,drop0Trailing=F), sep="_"))
}

if(!dir.exists(save_path)){
  dir.create(save_path, recursive = T)
}

seed_list=c(1:n_random_reps)+12345-1
#model_levels = c("MOFA+", "scVI", "concatenation", "WNN", "SNF", "CIMLR", "SG", "HVG")
R_model_levels = c("MOFA+", "scVI")
NR_model_levels = c("concatenation", "WNN", "SNF", "CIMLR", "SG", "HVG")

#model_type = c("integration","integration","integration","integration","integration","integration","control","control")
R_model_type = c("integration", "integration")
NR_model_type = c("integration", "integration", "integration", "integration", "control", "control")

#### load simulation parameters
params_df = read.csv(fs::path(data_path,"sim_params_df",ext="csv"),row.names = 1)

##### save parameters
spat_prob = rep_id = random_rep_id = model = type = ARI = AMI = DBI = DBI_truth = AIC = AIC0 = VN = TPR = PPV = FDR = FM = CSI = ACC = F1 = c()
print(nrow(params_df))

#### create Giotto Object
load("data_analysis/real_data/data/ST_MOB1.RData")
instrs = createGiottoInstructions(python_path = python_path)	
GO_test = createGiottoObject(raw_exprs=raw_expr, norm_expr=norm_expr, 	
                            spatial_locs=spatial, instructions=instrs)	
GO_test = filterGiotto(gobject = GO_test, expression_threshold = 1, gene_det_in_min_cells = 3, min_det_genes_per_cell = 200, expression_values = c("raw"), verbose = T)
GO_test@cell_metadata$clust_truth = annotation$group	
markers_truth_df = findMarkers_one_vs_all(GO_test, method = "mast", min_genes=0, cluster_column="clust_truth")	
markers_truth = markers_truth_df$genes

### compute metrics
for(i in 1:nrow(params_df)){#nrow(params_df)
  ## load data
  fname_attach = paste("sim", "spat", params_df$spat_prob[i], "rep", params_df$rep_id[i], sep = "_")
  load(fs::path(save_path, "data", fname_attach, ext = "RData"))

  for(t in 1:length(R_model_levels)){
    for(h in 1:length(seed_list)){
      clus_name = paste0("LC_",R_model_levels[t],"_clus")

       if(R_model_levels[t] == "MOFA+"){
        load(fs::path(save_path, "results", "cluster", "MOFA+",
                    paste(fname_attach, R_model_levels[t], "dim", ndims, "k", num_nn, "bin_LC", seed_list[h], sep = "_"), ext = "RData"))
        GO = MOFAp_bin_LC
        label_df = GO$cluster_label
      }else if(R_model_levels[t] == "scVI"){
        load(fs::path(save_path, "results", "cluster", "scVI", 
                    paste(fname_attach, R_model_levels[t], "sil_dim", ndims, "k", num_nn, "bin_LC", seed_list[h], sep = "_"),ext = "RData"))
        GO = scVI_bin_LC
        label_df = GO$cluster_label
      }
    
      #### aggregate values
      spat_prob = c(spat_prob, params_df$spat_prob[i])
      rep_id = c(rep_id, params_df$rep_id[i])
      model = c(model, R_model_levels[t])
      type = c(type, R_model_type[t])
      random_rep_id = c(random_rep_id, seed_list[h])
    
      DBI = c(DBI, index.DB(x=spatial, cl = as.numeric(label_df[,clus_name]))$DB)
      ARI = c(ARI, ARI(label_df[,clus_name], annotation$group))
      AMI = c(AMI, AMI(label_df[,clus_name], annotation$group))
      DBI_truth = c(DBI_truth, index.DB(x=spatial, cl = as.numeric(annotation$group))$DB)
      AIC_df = data.frame(sdimx = spatial$sdimx, sdimy = spatial$sdimy, truth = annotation$group, cl = label_df[,clus_name])
      AIC = c(AIC, spat_metric(AIC_df, comp="truth", method = "AIC"))
      AIC0 = c(AIC0, spat_metric(AIC_df, comp="null", method = "AIC"))
      VN = c(VN, spat_metric(AIC_df, method = "Var_norm"))
      DEA_metrics = spat_metric(AIC_df, DEA_method="gini", GO_object=GO_test, clust_labels=label_df[,clus_name], markers_truth=markers_truth, method="DEA") 
      TPR = c(TPR, DEA_metrics["TPR"])
      PPV = c(PPV, DEA_metrics["PPV"])
      FDR = c(FDR, DEA_metrics["FDR"])
      FM = c(FM, DEA_metrics["FM"])
      CSI = c(CSI, DEA_metrics["CSI"])
      ACC = c(ACC, DEA_metrics["ACC"])
      F1 = c(F1, DEA_metrics["F1"])
    }
    print(R_model_levels[t])
 }

  for(t in 1:length(NR_model_levels)){
    if(NR_model_levels[t] == "CIMLR"){
      clus_name = paste0(NR_model_levels[t],"_clus")
    }else{
      clus_name = paste0("LC_",NR_model_levels[t],"_clus")
    }
    if(NR_model_levels[t] == "HVG"){
      load(fs::path(save_path, "results", "cluster", "HVG",
                    paste(fname_attach, NR_model_levels[t], "dim", ndims, "k", num_nn, "LC", sep = "_"),ext = "RData"))
      GO = HVG_LC
      label_df = GO$cluster_label
    }else if(NR_model_levels[t] == "SNF"){
      load(fs::path(save_path, "results", "cluster", "SNF",
                    paste(fname_attach, NR_model_levels[t], "k", num_nn, "bin_LC", sep = "_"), ext="RData"))
      label_df = SNF_bin_LC
    }else if(NR_model_levels[t] == "WNN"){
      load(fs::path(save_path, "results", "cluster", "WNN",
                    paste(fname_attach, NR_model_levels[t], "dim", ndims, "k", num_nn, "bin_LC", sep = "_"), ext="RData"))
      label_df = WNN_bin_LC
    }else if(NR_model_levels[t] == "CIMLR"){
        load(fs::path(save_path, "results", "cluster", "CIMLR", 
                    paste(fname_attach, NR_model_levels[t], "dim", ndims, "bin_clus", sep = "_"), ext = "RData"))
        GO = CIMLR_bin_clus
        label_df = GO$cluster_label
    }else{
      load(fs::path(save_path, "results", "cluster", NR_model_levels[t],
                    paste(fname_attach, NR_model_levels[t], "dim", ndims, "k", num_nn, "bin_LC", sep = "_"),ext = "RData"))
      GO = eval(as.name(paste(NR_model_levels[t], "bin_LC", sep = "_")))
      print(names(GO))
      label_df = GO$cluster_label
      print(head(label_df))
    }

    #### aggregate values
    spat_prob = c(spat_prob, params_df$spat_prob[i])
    rep_id = c(rep_id, params_df$rep_id[i])
    model = c(model, NR_model_levels[t])
    type = c(type, NR_model_type[t])
    random_rep_id = c(random_rep_id, 0)

    DBI = c(DBI, index.DB(x=spatial, cl = as.numeric(label_df[,clus_name]))$DB)
    ARI = c(ARI, ARI(label_df[,clus_name], annotation$group))
    AMI = c(AMI, AMI(label_df[,clus_name], annotation$group))
    DBI_truth = c(DBI_truth, index.DB(x=spatial, cl = as.numeric(annotation$group))$DB)
    AIC_df = data.frame(sdimx = spatial$sdimx, sdimy = spatial$sdimy, truth = annotation$group, cl = label_df[,clus_name])
    AIC = c(AIC, spat_metric(AIC_df, comp="truth", method = "AIC"))
    AIC0 = c(AIC0, spat_metric(AIC_df, comp="null", method = "AIC"))
    VN = c(VN, spat_metric(AIC_df, method = "Var_norm"))
    DEA_metrics = spat_metric(AIC_df, DEA_method="gini", GO_object=GO_test, clust_labels=label_df[,clus_name], markers_truth=markers_truth, method="DEA")
    TPR = c(TPR, DEA_metrics["TPR"])
    PPV = c(PPV, DEA_metrics["PPV"])
    FDR = c(FDR, DEA_metrics["FDR"])
    FM = c(FM, DEA_metrics["FM"])
    CSI = c(CSI, DEA_metrics["CSI"])
    ACC = c(ACC, DEA_metrics["ACC"])
    F1 = c(F1, DEA_metrics["F1"])
    
    print(NR_model_levels[t])
  }
  print(i)
}

metrics = data.frame(spat_prob = spat_prob, rep_id = rep_id, random_rep_id = random_rep_id, model = model, type = type,
                         ARI = ARI, AMI = AMI, DBI = DBI, DBI_truth = DBI_truth, AIC = AIC, AIC0 = AIC0, VN = VN, TPR = TPR, PPV = PPV, FDR = FDR, FM = FM, CSI = CSI, ACC = ACC, F1 = F1)
#### save results
save(metrics, file = fs::path(save_path, "results", "cluster", paste("metrics_dim",ndims,"k",num_nn,sep="_"), ext = "RData"))
