#rm(list=ls())
setwd("/home/liyijun/ST_benchmark_01082020_re_1/")
source("functions/getCluster_Leiden_02222021.R")
library(Matrix)
library(data.table)
library(dplyr)

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

#save_path = fs::path(data_path, data_ref, paste("thres",thres,sep="_"),"data")
if(!dir.exists(save_path)){
  dir.create(save_path, recursive = T)
}

#### load simulation parameters
seed_list = c(1:n_random_reps)+12345-1
params_df = read.csv(fs::path(data_path,"sim_params_df",ext="csv"),row.names = 1)
method_proc = read.csv("method_proc.csv", row.names = 1)

#### load dataset
task_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
fname_attach = paste("sim", "spat", params_df$spat_prob[task_id], "rep", params_df$rep_id[task_id], sep = "_")
load(fs::path(save_path, "data", fname_attach, ext = "RData"))
cl0 = length(unique(annotation$group))

####### load HVG
HVG_path = fs::path(save_path, "results", "HVG")
HVG_dat = read.csv(fs::path(HVG_path, paste(fname_attach, "HVG", method_proc["HVG", 1], sep = "_"), ext = "csv"), row.names = 1)

####### load SG
SG_path = fs::path(save_path, "results", "SG")
SG_bin_dat = read.csv(fs::path(SG_path, paste(fname_attach, "SG_bin", method_proc["SG",1], sep = "_"), ext = "csv"), row.names = 1)

##### get hidden dimension
if(ndims==0){
  ndims = read.csv(fs::path(data_path, data_ref, paste("thres",thres,sep="_"), "results", "scVI", paste(fname_attach, "hyperparams", sep = "_"), ext = "csv"))
  ndims = ndims$n_latent
}

####### run clustering analysis
### HVG
HVG_LC = getCluster_Leiden(gene_expr=HVG_dat, 
                               spat_mat=spatial, 
                               data_type = "HVG", 
                               reduced = FALSE,
                               pc_dims = min(ndims, nrow(HVG_dat)), 	
                               ncl_truth = cl0, 	
                               res = 0.5, res_step = 0.1, max_iter = 100, num_iter = -1,
                               py_path = python_path, num_neighbors=num_nn, pb = T)
save_path_HVG = fs::path(save_path, "results", "cluster", "HVG")
if(!dir.exists(save_path_HVG)){
  dir.create(save_path_HVG, recursive = T)
}
save(HVG_LC, file = fs::path(save_path_HVG, paste(fname_attach, "HVG_dim", ndims, "k", num_nn, "LC", sep = "_"), ext="RData"))

##### SG
SG_bin_LC = getCluster_Leiden(gene_expr=SG_bin_dat, 
                                  spat_mat=spatial, 
                                  data_type = "SG", 
                                  reduced = FALSE,
                                  pc_dims = min(ndims, nrow(SG_bin_dat)), 	
                               ncl_truth = cl0, 	
                               res = 0.5, res_step = 0.1, max_iter = 100, num_iter = -1,
                                  py_path = python_path, num_neighbors = num_nn, pb = T)
save_path_SG = fs::path(save_path, "results", "cluster", "SG")
if(!dir.exists(save_path_SG)){
  dir.create(save_path_SG, recursive = T)
}
save(SG_bin_LC, file = fs::path(save_path_SG, paste(fname_attach, "SG_dim", ndims, "k", num_nn, "bin_LC", sep = "_"),ext="RData"))

### full
concatenation_data_name = paste(fname_attach, "concatenation_bin", sep = "_")
load(fs::path(save_path, "results", "concatenation", concatenation_data_name, ext="RData"))
concatenation_bin_LC = getCluster_Leiden(gene_expr=concatenation_bin, 
                                    spat_mat=spatial, 
                                    data_type = "concatenation", 
                                    reduced = FALSE,
                                    pc_dims = min(ndims,nrow(concatenation_bin)), 	
                               ncl_truth = cl0, 	
                               res = 0.5, res_step = 0.1, max_iter = 100, num_iter = -1,
                                    py_path = python_path, num_neighbors=num_nn, pb = T)
save_path_concatenation = fs::path(save_path, "results", "cluster", "concatenation")
if(!dir.exists(save_path_concatenation)){
  dir.create(save_path_concatenation, recursive = T)
}
print(save_path_concatenation)
save(concatenation_bin_LC,
     file = fs::path(save_path_concatenation, paste(fname_attach, "concatenation_dim", ndims, "k", num_nn, "bin_LC", sep = "_"), ext="RData"))

### MOFA+
for(i in 1:length(seed_list)){
  MOFAp_fname_attach = paste(fname_attach, "MOFA+_dim", ndims, "bin", seed_list[i], sep = "_")
  MOFAp_bin = read.csv(file = fs::path(save_path, "results", "MOFA+", MOFAp_fname_attach, ext="csv"), row.names = 1)
  dim(MOFAp_bin)
  dim(annotation)
  row.names(MOFAp_bin) = annotation$cell_ID #assign cell names
  MOFAp_bin_LC = getCluster_Leiden(gene_expr=t(MOFAp_bin), 
                                     spat_mat=spatial, 
                                     data_type = "MOFA+", 
                                     reduced = TRUE,
                                     pc_dims = ndims, 	
                               ncl_truth = cl0, 	
                               res = 0.5, res_step = 0.1, max_iter = 100, num_iter = -1,
                                     py_path = python_path, num_neighbors=num_nn, pb = T)
  save_path_MOFAp = fs::path(save_path, "results", "cluster", "MOFA+")
  if(!dir.exists(save_path_MOFAp)){
    dir.create(save_path_MOFAp, recursive = T)
  }
  save(MOFAp_bin_LC,
     file = fs::path(save_path_MOFAp, paste(fname_attach, "MOFA+_dim", ndims, "k", num_nn, "bin_LC", seed_list[i], sep = "_"), ext="RData"))
}

### scVI 
for(i in 1:length(seed_list)){
  scVI_fname_attach = paste(fname_attach,"scVI_dim",ndims, seed_list[i], sep="_")
  scVI_bin = read.csv(fs::path(save_path, "results", "scVI", scVI_fname_attach, ext="csv"),header=F)
  head(scVI_bin)
  head(annotation)
  row.names(scVI_bin) = annotation$cell_ID
  scVI_bin_LC = getCluster_Leiden(gene_expr=t(scVI_bin), 
                                     spat_mat=spatial, 
                                     data_type = "scVI", 
                                     reduced = TRUE,
                                     pc_dims = ndims,	
                               ncl_truth = cl0, 	
                               res = 0.5, res_step = 0.1, max_iter = 100, num_iter = -1,
                                     py_path = python_path, num_neighbors=num_nn, pb = T)
  save_path_scVI = fs::path(save_path, "results", "cluster", "scVI")
  if(!dir.exists(save_path_scVI)){
    dir.create(save_path_scVI, recursive = T)
  }
  save(scVI_bin_LC, 
     file = fs::path(save_path_scVI, paste(fname_attach, "scVI_sil_dim", ndims, "k", num_nn, "bin_LC", seed_list[i], sep = "_"), ext="RData"))
}

### SNF
save_path_SNF = fs::path(save_path, "results", "cluster", "SNF")
SNF_bin = read.csv(fs::path(save_path_SNF, paste(fname_attach, "SNF_k", num_nn, "bin_LC", sep = "_"), ext="csv"))
SNF_bin_LC = data.frame(cell_ID = annotation$cell_ID, LC_SNF_clus = SNF_bin$x) 
head(SNF_bin_LC)
table(SNF_bin_LC$LC_SNF_clus)
if(!dir.exists(save_path_SNF)){
  dir.create(save_path_SNF, recursive = T)
}
save(SNF_bin_LC, file = fs::path(save_path_SNF, paste(fname_attach, "SNF_k", num_nn, "bin_LC", sep = "_"), ext="RData"))
