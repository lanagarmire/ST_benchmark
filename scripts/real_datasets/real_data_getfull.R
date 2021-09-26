rm(list=ls())
library(data.table)
library(Giotto)
setwd("/home/liyijun/ST_benchmark_01082020_re_1")

#### get arguments
args = commandArgs(trailingOnly = TRUE)
data_path = as.character(args[1])
data_ref = as.character(args[2])
diff_in_cov = as.numeric(args[3])
svg_pval_thres = as.numeric(args[4])
ndims = as.numeric(args[5])
scVI_error = as.character(args[6])
SNF_metric = as.character(args[7])
WNN_smooth = as.numeric(args[8])
CIMLR_k = as.numeric(args[9])
save_name = as.character(args[10])

### set parameters
if(save_name == "diff_in_cov"){
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


method_proc = read.csv("method_proc.csv",row.names = 1)
method_name = "full"
proc = method_proc[method_name,]
params_df = read.csv(fs::path(data_path,data_ref,"real_data_names",ext="csv"),row.names = 1)

#### load dataset
task_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
fname_attach = params_df$name[task_id]
load(fs::path(data_path, data_ref, "data", fname_attach, ext = "RData"))
cl0 = length(unique(annotation$group))

####### load HVG
HVG_path = fs::path(save_path, "results", "HVG")
HVG_dat = read.csv(fs::path(HVG_path, paste(fname_attach, "HVG", proc, sep = "_"), ext = "csv"), row.names = 1)

####### load SG
SG_path = fs::path(save_path, "results", "SG")
SG_bin_dat = read.csv(fs::path(SG_path, paste(fname_attach, "SG_bin", proc, sep = "_"), ext = "csv"), row.names = 1)

### full -- binspect
concatenation_bin = rbind(as.data.frame(HVG_dat), as.data.frame(SG_bin_dat))
dim(concatenation_bin)

save_full_path = fs::path(save_path, "results", "concatenation")
if(!dir.exists(save_full_path)){
  dir.create(save_full_path, recursive = T)
}
concatenation_data_name = paste(fname_attach, "concatenation_bin", sep = "_")
save(concatenation_bin, file = fs::path(save_full_path, concatenation_data_name, ext = "RData"))
fwrite(concatenation_bin, file = fs::path(save_full_path, concatenation_data_name, ext = "csv"))
