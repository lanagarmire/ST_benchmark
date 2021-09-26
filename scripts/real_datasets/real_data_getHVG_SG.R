library(Giotto)
library(dplyr)
library(Matrix)
library(data.table)
setwd("/home/liyijun/ST_benchmark_01082020_re_1")

#### get arguments from bash scripts
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


params_df = read.csv(fs::path(data_path,data_ref,"real_data_names",ext="csv"),row.names = 1)

#### load dataset
task_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
fname_attach = params_df$name[task_id]
load(fs::path(data_path, data_ref, "data", fname_attach, ext = "RData"))

##### get all HVGs
source("functions/getHVG_11182020.R")
HVG_res = getHVG(raw_mat = raw_expr,
                 spat_mat = spatial,
                 wd = fs::path(save_path, "results", "HVG"),
                 fname = fname_attach, python_path = python_path, is_docker = F, pb=T,
                 exp_thres = 1, min_cells = 3, min_genes = floor(0.05*nrow(raw_expr)), diff_cov = diff_in_cov)
HVG_original = HVG_res$HVG

###### get all SG
source("functions/getSG_01182021.R")
SG_res = getSG(raw_mat = raw_expr, spat_mat = spatial,
               wd = fs::path(save_path, "results", "SG"), fname = fname_attach,
               python_path = python_path, is_docker=F, pb=T,
               exp_thres = 1, min_cells = 3, min_genes = floor(0.05*nrow(raw_expr)), annot_df = annotation, spat_gene_method = "binspect",
               max_dist_dly = "auto")
SG_original_df = SG_res$SG %>% filter(adj.p.value < svg_pval_thres)
SG_original = SG_original_df$genes

##### select unique HVG and SG
genes_SG_only_sig = setdiff(SG_original, HVG_original)
genes_HVG_only = setdiff(HVG_original, SG_original)
genes_HVG_SG = intersect(HVG_original, SG_original)

######## save HVG
HVG_path = fs::path(save_path, "results", "HVG")
if(!dir.exists(HVG_path)){
  dir.create(HVG_path, recursive=T)
}
HVG_GO = HVG_res$GO
HVG = genes_HVG_only
save(HVG_GO, HVG, file = fs::path(HVG_path,paste(fname_attach,"HVG",sep="_"),ext="RData"))
fwrite(data.frame(HVG_GO@norm_expr[HVG,]), file = fs::path(HVG_path, paste(fname_attach,"HVG_norm",sep="_"),ext="csv"),sep=",",row.names=T)
fwrite(data.frame(as.matrix(HVG_GO@raw_exprs[HVG,])), file = fs::path(HVG_path, paste(fname_attach,"HVG_raw",sep="_"),ext="csv"),sep=",",row.names=T)

######### save SG
SG_path = fs::path(save_path, "results", "SG")
if(!dir.exists(SG_path)){
  dir.create(SG_path, recursive = T)
}
SG_bin_GO = SG_res$GO
SG_bin = SG_res$SG
spat_genes_bin = genes_SG_only_sig
save(SG_bin_GO, SG_bin, spat_genes_bin, file = fs::path(SG_path, paste(fname_attach, "SG_bin", sep = "_"), ext="RData"))
fwrite(data.frame(SG_bin_GO@norm_expr[spat_genes_bin,]), file = fs::path(SG_path, paste(fname_attach, "SG_bin_norm", sep = "_"), ext = "csv"), sep=",", row.names = T)
fwrite(data.frame(as.matrix(SG_bin_GO@raw_exprs[spat_genes_bin,])), file = fs::path(SG_path, paste(fname_attach, "SG_bin_raw", sep = "_"), ext = "csv"), sep=",", row.names = T)

