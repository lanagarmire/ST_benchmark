rm(list=ls())
library(data.table)
library(Giotto)
setwd("/home/liyijun/ST_benchmark_01082020")

args = commandArgs(trailingOnly = TRUE)
data_path = as.character(args[1])
setwd(data_path)

#### load simulation parameters
method_proc = read.csv("/home/liyijun/ST_benchmark_01082020/method_proc.csv",row.names = 1)
method_name = "full"
proc = method_proc[method_name,]
params_df = read.csv(fs::path("data", "cells", ext="csv"), row.names = 1)

#### load dataset
task_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))*7
fname_attach = params_df$cells[task_id]
rep = params_df$reps[task_id]

####### load HVG
HVG_dat = read.csv(fs::path("data", paste("data_cells", fname_attach, proc, 1, rep, sep = "_"), ext = "csv"), row.names = 1)
colnames(HVG_dat) = paste("cells", 1:ncol(HVG_dat), sep="_")

####### load SG
SG_bin_dat = read.csv(fs::path("data", paste("data_cells", fname_attach, proc, 2, rep, sep = "_"), ext = "csv"), row.names = 1)
colnames(SG_bin_dat) = paste("cells", 1:ncol(SG_bin_dat), sep="_")

### full -- binspect
tic = Sys.time()
concatenation_bin = rbind(as.data.frame(HVG_dat), as.data.frame(SG_bin_dat))
toc = Sys.time()
elapse = toc-tic

if(!dir.exists(fs::path("results", "concatenation"))){
  dir.create(fs::path("results", "concatenation"), recursive = T)
}

save(elapse, file = fs::path("results", "concatenation", paste(fname_attach, "concatenation", rep, "time", sep="_"), ext="RData"))

