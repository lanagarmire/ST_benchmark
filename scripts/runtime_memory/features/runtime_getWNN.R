rm(list = ls())
setwd("/home/liyijun/ST_benchmark_01082020/")

library(data.table)
library(Seurat)
library(dplyr)

#### get arguments from bash scripts
args = commandArgs(trailingOnly = TRUE)
data_path = as.character(args[1])
ndims = as.numeric(args[2])

scale = F
center = F
sf = 6000

setwd(data_path)
params_df = read.csv(fs::path("data", "features", ext="csv"),row.names = 1)
method_proc = read.csv("/home/liyijun/ST_benchmark_01082020/method_proc.csv", row.names = 1)
method_name = "WNN"
proc = method_proc[method_name, 1]

task_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
data1 = read.csv(fs::path("data", paste("data_features", params_df$features[task_id], proc, 1, params_df$reps[task_id], sep="_"), ext="csv"), row.names=1)
#colnames(data1) = NULL
data1_name = "HVG"

data2 = read.csv(fs::path("data", paste("data_features", params_df$features[task_id], proc, 2, params_df$reps[task_id], sep="_"), ext="csv"), row.names=1)
#colnames(data2) = NULL
data2_name = "SG"

data1_npcs = 100
data2_npcs = 100
print(data1_npcs)
print(data2_npcs)

tic = Sys.time()
combine_SO = CreateSeuratObject(counts = data1, assay=data1_name)
data2_assay = CreateAssayObject(counts = data2)
combine_SO[[data2_name]] = data2_assay

#preprocess and dimension reduction
DefaultAssay(combine_SO) = data1_name
VariableFeatures(combine_SO) = rownames(data1)
combine_SO = NormalizeData(combine_SO, scale.factor = sf)%>%
  ScaleData(do.scale=scale, do.center=center) %>% 
  RunPCA(npcs = data1_npcs, reduction.name=paste0(data1_name,"_pca")) 

DefaultAssay(combine_SO) = data2_name
VariableFeatures(combine_SO) = rownames(data2)
combine_SO = NormalizeData(combine_SO, scale.factor = sf)%>%
  ScaleData(do.scale = scale, do.center = center) %>% 
  RunPCA(npcs = data2_npcs, reduction.name=paste0(data2_name,"_pca")) 

combine_SO = FindMultiModalNeighbors(combine_SO,
                                     reduction.list = list(paste0(data1_name,"_pca"),
                                                           paste0(data2_name,"_pca")),
                                     dims.list = list(1:ndims, 1:ndims),
                                     k.nn = 30, smooth = FALSE)
toc = Sys.time()
elapse = toc-tic

if(!dir.exists(fs::path("results","WNN"))){
  dir.create(fs::path("results","WNN"), recursive = T)
}

save(elapse, file = fs::path("results", "WNN", paste(params_df$features[task_id], "WNN", params_df$reps[task_id], "time", sep="_"), ext="RData"))

