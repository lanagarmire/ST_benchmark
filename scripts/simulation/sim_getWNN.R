rm(list = ls())
setwd("/home/liyijun/ST_benchmark_01082020_re_1/")
source("functions/combine_Seuratv4_02192021.R")
library(data.table)

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

#### load simulation parameters
params_df = read.csv(fs::path(data_path,"sim_params_df",ext="csv"),row.names = 1)
method_proc = read.csv("method_proc.csv", row.names = 1)
method_name = "WNN"
proc = method_proc[method_name, 1]

#### load dataset
task_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
fname_attach = paste("sim", "spat", params_df$spat_prob[task_id], "rep", params_df$rep_id[task_id], sep = "_")
load(fs::path(save_path, "data", fname_attach, ext = "RData"))
cl0 = length(unique(annotation$group))

####### load HVG
HVG_path = fs::path(save_path, "results", "HVG")
HVG_dat = read.csv(fs::path(HVG_path, paste(fname_attach, "HVG", proc, sep = "_"), ext = "csv"), row.names = 1)

####### load SG
SG_path = fs::path(save_path, "results", "SG")
SG_bin_dat = read.csv(fs::path(SG_path, paste(fname_attach, "SG_bin", proc, sep = "_"), ext = "csv"), row.names = 1)

dim(HVG_dat)
dim(SG_bin_dat)

##### get hidden dimension
if(ndims==0){
  ndims = read.csv(fs::path(data_path, data_ref, paste("thres",thres,sep="_"), "results", "scVI", paste(fname_attach, "hyperparams", sep = "_"), ext = "csv"))
  ndims = ndims$n_latent
}
print(ndims)

## run WNN 
WNN_bin = combine_Seuratv4(data1 = HVG_dat, data2 = SG_bin_dat,                           
                           res = 0.5, res_step = 0.1, max_iter=100,
                           data1_name = "HVG", data2_name = "SG",
                           data1_npcs = 100, data2_npcs=100,
                           data1_dims = ndims, data2_dims = ndims, n_neighbors = num_nn, smooth = as.logical(WNN_smooth),
                           ncl_truth = cl0, n.iter=-1)
table(WNN_bin$SO@meta.data$seurat_clusters)

### save WNN results
WNN_folder_name = "WNN"
WNN_data_name = paste(fname_attach, "WNN_dim", ndims, "k", num_nn, "bin", sep = "_")
WNN_save_path = fs::path(save_path,"results",WNN_folder_name)
if(!dir.exists(WNN_save_path)){
  dir.create(WNN_save_path, recursive = T)
}
save(WNN_bin, file = fs::path(WNN_save_path, WNN_data_name, ext = "RData"))

### save WNN clustering results
WNN_bin_LC = data.frame(cell_ID = rownames(WNN_bin$SO@meta.data),
                            LC_WNN_clus = WNN_bin$SO@meta.data$seurat_clusters) 
save_path_clus = fs::path(save_path, "results", "cluster", WNN_folder_name)
if(!dir.exists(save_path_clus)){
  dir.create(save_path_clus, recursive = T)
}

save(WNN_bin_LC,
     file = fs::path(save_path_clus, paste(fname_attach, method_name, "dim", ndims, "k", num_nn, "bin_LC", sep = "_"), ext="RData"))
