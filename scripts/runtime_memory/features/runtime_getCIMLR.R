rm(list=ls())
setwd("/home/liyijun/ST_benchmark_01082020/")
library(data.table)

#### get arguments from bash scripts
args = commandArgs(trailingOnly = TRUE)
data_path = as.character(args[1])
ndims = as.numeric(args[2])
setwd(data_path)

#### load simulation parameters
params_df = read.csv(fs::path("data", "features", ext="csv"),row.names = 1)
method_proc = read.csv("/home/liyijun/ST_benchmark_01082020/method_proc.csv", row.names = 1)
method_name = "CIMLR"
proc = method_proc[method_name, 1]

#### load dataset
task_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
fname_attach = params_df$features[task_id]
rep = params_df$reps[task_id]
cl0 = 8

####### load HVG
hvg_data = read.csv(fs::path("data", paste("data_features", fname_attach, proc, 1, rep, sep = "_"), ext = "csv"), row.names = 1)

####### load SG
svg_data = read.csv(fs::path("data", paste("data_features", fname_attach, proc, 2, rep, sep = "_"), ext = "csv"), row.names = 1)

### run CIMLR
library(CIMLR)
dat_list = list(hvg_data, svg_data)
  
tic = Sys.time()
res = CIMLR(dat_list, cl0, no.dim = ndims, k=10, cores.ratio = 0)
toc = Sys.time()
elapse = toc-tic  
 
if(!dir.exists(fs::path("results", method_name))){
  dir.create(fs::path("results", method_name),recursive = T)
}

save(elapse, file = fs::path("results", method_name, paste(fname_attach, method_name,rep, "time", sep="_"), ext="RData"))

