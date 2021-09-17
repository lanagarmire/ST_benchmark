#.libPaths( c( "/nfs/home/jingzhe/R/x86_64-pc-linux-gnu-library/4.0" , .libPaths() ) )
setwd('/home/liyijun/ST_benchmark_01082020') #change
library(MOFA2)
library(data.table)
library(dplyr)

#### get arguments from bash scripts
args = commandArgs(trailingOnly = TRUE)
data_path = as.character(args[1])
#data_path = "simulation/simulation_04302021"
thres = as.numeric(args[2])
#thres=0.6
data_ref = as.character(args[3])
#data_ref="ST_MOB1"
diff_in_cov = as.numeric(args[4])
svg_pval_thres = as.numeric(args[5])
#num_nn = as.numeric(args[6])
ndims = as.numeric(args[6])
scVI_error = as.character(args[7])
SNF_metric = as.character(args[8])
WNN_smooth = as.numeric(args[9])
CIMLR_k = as.numeric(args[10])
save_name = as.character(args[11])
n_random_reps = as.numeric(args[12])

### load simulation parameters
params_df = read.csv(fs::path(data_path,"sim_params_df",ext="csv"),row.names = 1)
method_proc = read.csv("method_proc.csv", row.names=1)
method_name = "MOFA+"

#### get dataname ##change
task_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
data_name = paste("sim", "spat", params_df$spat_prob[task_id], "rep", params_df$rep_id[task_id], sep = "_")

## get cell name reference ##change
if(save_name == "thres"){
  save_path = fs::path(data_path, data_ref, paste(save_name, format(thres,drop0Trailing=F),sep="_"))
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

file = fread(fs::path(save_path, "results", "SG", paste(data_name, "SG_bin", method_proc[method_name,1], sep="_"), ext ="csv"))

### get MOFA+ output from .py script
#change date
if(ndims == 0){
  ndims = read.csv(fs::path(save_path, "results", "scVI", paste(data_name, "hyperparams", sep = "_"), ext = "csv"))
  ndims = ndims$n_latent
}

seed_list = c(1:n_random_reps)+12345-1
for(i in 1:length(seed_list)){
  model_path = fs::path(save_path, "results", method_name, paste(data_name, method_name, "dim", ndims, "bin", seed_list[i], sep = "_"), ext = "hdf5") 
  model <- load_model(model_path, remove_inactive_factors = F)
  factors <- get_factors(model, factors = "all")
  factors <- factors$group0
  rownames(factors) = colnames(file)[2:ncol(file)]

  write.csv(factors, fs::path(save_path, "results", method_name, paste(data_name, method_name, "dim", ndims, "bin", seed_list[i], sep = "_"), ext = "csv"),quote=FALSE)
}
#change date
