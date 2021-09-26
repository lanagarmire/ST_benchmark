rm(list=ls())
setwd("/home/liyijun/ST_benchmark_01082020_re_1/")
source("functions/run_cimlr_03232021.R")
library(data.table)

#### get arguments from bash scripts
args = commandArgs(trailingOnly = TRUE)
data_path = as.character(args[1])
thres = as.numeric(args[2])
data_ref = as.character(args[3])
diff_in_cov = as.numeric(args[4])
svg_pval_thres = as.numeric(args[5])
ndims = as.numeric(args[6])
scVI_error = as.character(args[7])
SNF_metric = as.character(args[8])
WNN_smooth = as.numeric(args[9])
CIMLR_k = as.numeric(args[10])
save_name = as.character(args[11])
n_random_reps = as.numeric(args[12])
cores_ratio = as.numeric(args[13])

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

#save_path = fs::path(data_path, data_ref, paste("thres",thres,sep="_"),"data")
if(!dir.exists(save_path)){
  dir.create(save_path, recursive = T)
}

#### load simulation parameters
params_df = read.csv(fs::path(data_path,"sim_params_df",ext="csv"),row.names = 1)
method_proc = read.csv("method_proc.csv", row.names = 1)
method_name = "CIMLR"
proc = method_proc[method_name, 1]

#### load dataset
task_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
fname_attach = paste("sim", "spat", params_df$spat_prob[task_id], "rep", params_df$rep_id[task_id], sep = "_")
load(fs::path(save_path, "data", fname_attach, ext = "RData"))
cl0 = length(unique(annotation$group))

####### load HVG
HVG_path = fs::path(save_path,"results", "HVG")
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

### run CIMLR
set.seed(12345)

CIMLR_bin_clus = run_cimlr(hvg_data = HVG_dat, svg_data = SG_bin_dat, ncl_truth = cl0, no.dim = ndims, cores.ratio=cores_ratio)
table(CIMLR_bin_clus$cluster_label$CIMLR_clus)

CIMLR_save_path = fs::path(save_path, "results", "cluster", method_name)
if(!dir.exists(CIMLR_save_path)){
  dir.create(CIMLR_save_path,recursive = T)
}
save(CIMLR_bin_clus, file = fs::path(CIMLR_save_path, paste(fname_attach, "CIMLR", "dim", ndims, "bin_clus", sep = "_"), ext="RData"))

