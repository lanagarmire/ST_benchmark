#rm(list=ls())
library(Giotto)
library(dplyr)
library(Matrix)
library(data.table)
setwd("/home/liyijun/ST_benchmark_01082020")

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
#load(fs::path(data_path, data_ref, paste(data_ref, "gene_grp_df", sep="_"), ext="RData"))
params_df = read.csv(fs::path(data_path,"sim_params_df",ext="csv"),row.names = 1)

task_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
spat_prob = params_df$spat_prob[task_id]
rep_id = params_df$rep_id[task_id]
fname_attach = paste("sim", "spat", params_df$spat_prob[task_id], "rep", params_df$rep_id[task_id],sep = "_")
seed = as.integer(12345+task_id)

load(fs::path(data_path, "real_datasets", data_ref, ext="RData"))
### get HVG
source("functions/getHVG_11182020.R")
HVG_res = getHVG(raw_mat = raw_expr,
                 spat_mat = spatial,
                 wd = fs::path(data_path,data_ref),
                 fname = paste(data_ref, diff_in_cov, "HVG", sep = "_"), python_path = python_path, is_docker = F, pb=T,
                 exp_thres = 1, min_cells = 3, min_genes = floor(0.05*nrow(raw_expr)), diff_cov = diff_in_cov)
#removed no cells or genes
HVG_original = HVG_res$HVG

### only work with solely SVG
source("functions/getSG_01182021.R")
#load(fs::path("data_analysis/real_data/data",data_ref, ext = "RData"))
#load(fs::path("data_analysis/real_data/04302021/results/HVG", paste(data_ref, "HVG", sep = "_"), ext = "RData"))
#load(fs::path("data_analysis/real_data/04302021/results/SG", paste(data_ref, "SG_bin", sep = "_"), ext = "RData"))
SG_res = getSG(raw_mat = raw_expr, spat_mat = spatial,
               wd = fs::path(data_path, data_ref), fname = paste(data_ref, svg_pval_thres, "SG", sep="_"),
               python_path = python_path, is_docker=F, pb=T,
               exp_thres = 1, min_cells = 3, min_genes = floor(0.05*nrow(raw_expr)), annot_df = annotation, spat_gene_method = "binspect",
               max_dist_dly = "auto")

SG_original_df = SG_res$SG %>% filter(adj.p.value < svg_pval_thres)
SG_original = SG_original_df$genes
genes_SG_only_sig = setdiff(SG_original, HVG_original)
genes_HVG_only = setdiff(HVG_original, SG_original)
genes_HVG_SG = intersect(HVG_original, SG_original)

#genes_SG_only = spat_genes_bin[which(!(spat_genes_bin %in% HVG))]
#genes_SG_only_sig = SG_bin %>% filter(genes %in% genes_SG_only) %>% filter(adj.p.value < 0.05)
#genes_SG_only_sig = genes_SG_only_sig$genes
#genes_HVG_SG = intersect(HVG, spat_genes_bin)
#genes_HVG_only = HVG[which(!HVG %in% genes_HVG_SG)]

### compute gene proportions
cell_id = HVG_res$GO@cell_ID
gene_id = HVG_res$GO@gene_ID
grp_id = sort(unique(annotation$group))
expr_data = HVG_res$GO@norm_expr

gene_grp_id = rep(NA, length(gene_id))
grp_perc = matrix(NA, length(gene_grp_id), length(grp_id))

for(i in 1:nrow(grp_perc)){
  gene_name = gene_id[i]
  #subset the original normalized expression of the selected gene
  gene_vector = expr_data[rownames(expr_data) == gene_name, ]
  #sort in order of high to low expressed (normalized)
  sort_expr_gene = sort(gene_vector, decreasing = T)
  sort_expr_gene_non_zero = names(sort_expr_gene[sort_expr_gene != 0])

  for(j in 1:length(grp_id)){
    cell_id_subset_df = annotation %>% filter(group==grp_id[j])
    cell_id_subset = cell_id_subset_df$cell_ID
    n_grp = length(cell_id_subset)
    grp_perc[i,j] = length(intersect(sort_expr_gene_non_zero[1:n_grp], cell_id_subset))
  }
  grp_perc[i,] = grp_perc[i,]/sum(grp_perc[i,])
  gene_grp_id[i] = which.max(grp_perc[i,])
}

gene_grp_df = data.frame(gene_name = gene_id, gene_grp_id = gene_grp_id, stringsAsFactors = F)
gene_grp_df = cbind(gene_grp_df, grp_perc)
head(gene_grp_df)

## save gene group id
save_name = paste(data_ref, "gene_grp_df", sep="_")
save(gene_grp_df, file = fs::path(data_path, data_ref, save_name, ext="RData"))

#### set base parameters
cl0 = length(unique(annotation$group))
n_grps = table(annotation$group)

### assign gene pattern ids
clear = rep(NA, nrow(gene_grp_df))
for(i in 1:nrow(gene_grp_df)){
  if(max(gene_grp_df[i, as.character(sort(unique(annotation$group)))]) > thres){
    clear[i] = "yes"
  }else{
    clear[i] = "no"
  }
}
gene_grp_df$clear = clear 
head(gene_grp_df)
df2 = gene_grp_df%>%filter(gene_name %in% genes_SG_only_sig)
table(df2$clear)
df2_yes = df2%>%filter(clear=="yes")

#### simulate new patterns
GO_test = SG_res$GO
orig_raw = GO_test@raw_exprs
orig_norm = expr_data
new_raw = orig_raw
set.seed(seed)

for(j in 1:length(grp_id)){
  #select the pattern
  cell_id_subset_df = annotation %>% filter(cell_ID %in% cell_id) %>% filter(group==grp_id[j])
  cell_id_subset = cell_id_subset_df$cell_ID
  
  #subset the genes that belong to this pattern
  gene_grp_df_subset = df2_yes%>%filter(gene_grp_id==grp_id[j])
  gene_grp_df_subset = gene_grp_df_subset$gene_name
  
  if(length(gene_grp_df_subset)>0){
    for(i in 1:length(gene_grp_df_subset)){
      gene = gene_grp_df_subset[i]
      
      gobject = simulateOneGenePatternGiottoObject(GO_test, pattern_name = "pattern", pattern_cell_ids = cell_id_subset, 
                                                   gene_name = gene, spatial_prob = spat_prob, show_pattern = F)
      new_raw[row.names(new_raw)==gene,] = gobject@raw_exprs[row.names(gobject@raw_exprs)==gene,]
    }
  }
}

### create raw data (no need to filter)
myinst = createGiottoInstructions(python_path = python_path, is_docker = F)
GO = createGiottoObject(raw_exprs = new_raw, 
                        spatial_locs = GO_test@spatial_locs%>%select(c("sdimx","sdimy")), 
                        instructions = myinst)

### filter 
#GO = filterGiotto(gobject = GO,
#                  expression_threshold = 1, gene_det_in_min_cells = 1, min_det_genes_per_cell = 1,
#                  expression_values = c('raw'), verbose = T)

### normalize
GO = normalizeGiotto(gobject = GO, scalefactor = 6000)
norm_expr = GO@norm_expr

### plot new simulated patterns
for(i in 1:nrow(df2_yes)){
  gene = df2_yes$gene_name[i]
  spatGenePlot2D(GO, expression_values = 'norm', genes = gene, point_shape = 'border', point_border_stroke = 0.1,
                 show_network = F, network_color = 'lightgrey', point_size = 2.5, cow_n_col = 1, show_plot = F,
                 save_plot = T, save_param = list(save_dir = fs::path(save_path, "data"), 
                                                  save_folder = paste('simulated_spat',spat_prob,"rep",rep_id,sep="_"), 
                                                  save_name = paste0(gene,'_sim'),
                                                  base_width = 9, base_height = 7, units = 'cm'))
}

###### add annotations
GO = addStatistics(gobject = GO)
GO@cell_metadata = inner_join(GO@cell_metadata, 
                                annotation%>%filter(cell_ID %in% cell_id)%>%select(c("cell_ID","group")), by = "cell_ID")
annotation = GO@cell_metadata

#run dimension reduction: PCA, t-SNE, UMAP
GO = Giotto::runPCA(gobject = GO, expression_values = "normalized", genes_to_use = "all", 
                    scale_unit = F, center = F, method="factominer")
n_dims = 15
n_dims = min(dim(GO@dimension_reduction$cells$pca$pca$coordinates)[2], n_dims)
GO = Giotto::runUMAP(gobject = GO, dimensions_to_use = 1:n_dims, seed_number = seed)
GO = Giotto::runtSNE(gobject = GO, dimensions_to_use = 1:n_dims, seed_number = seed, check_duplicates=FALSE)

PCA = GO@dimension_reduction$cells$pca$pca$coordinates
UMAP = GO@dimension_reduction$cells$umap$umap$coordinates
tSNE = GO@dimension_reduction$cells$tsne$tsne$coordinates

raw_expr = new_raw
spatial = GO_test@spatial_locs%>%select(c("sdimx","sdimy"))
save(raw_expr, norm_expr, spatial, annotation, PCA, UMAP, tSNE, GO, 
     file = fs::path(save_path,"data", fname_attach, ext = "RData"))
cl0 = length(unique(annotation$group))
fwrite(data.frame(cl0 = cl0), file = fs::path(save_path,"data", paste(fname_attach, "truth", sep = "_"), ext = "csv"), sep=",", row.names = T)

#save plots
Giotto::plotPCA(GO,cell_color = "group",
                cell_color_code = Giotto:::getDistinctColors(cl0),
                point_size=2.5,save_param=list(save_dir = fs::path(save_path,"data"), save_folder = "PCA",
                                               save_name = paste0(fname_attach,"_PCA"),
                                               save_format = "png"),
                save_plot=TRUE, title="PCA of simulated data",
                legend_symbol_size = 2.5)
Giotto::plotUMAP(GO,cell_color = "group",
                 cell_color_code = Giotto:::getDistinctColors(cl0),
                 point_size=2.5,save_param=list(save_dir = fs::path(save_path,"data"), save_folder = "UMAP",
                                                save_name = paste0(fname_attach, "_UMAP"),
                                                save_format = "png"),
                 save_plot=TRUE, title="UMAP of simulated data",
                 legend_symbol_size = 2.5)
Giotto::plotTSNE(GO,cell_color = "group",
                 cell_color_code = Giotto:::getDistinctColors(cl0),
                 point_size=2.5,save_param=list(save_dir = fs::path(save_path,"data"), save_folder = "tSNE",
                                                save_name = paste0(fname_attach, "_tSNE"),
                                                save_format = "png"),
                 save_plot=TRUE, title="t-SNE of simulated data",
                 legend_symbol_size = 2.5)
Giotto::spatPlot(GO,cell_color = "group",
                 cell_color_code = Giotto:::getDistinctColors(cl0),
                 point_size=2.5,save_param=list(save_dir = fs::path(save_path,"data"), save_folder = "spat",
                                                save_name = paste0(fname_attach, "_spat"),
                                                save_format = "png"),
                 save_plot=TRUE, title="spatial coordinates of simulated data",
                 legend_symbol_size = 2.5)

####### save SG
#source("functions/getSG_01182021.R")
#sim_SG_path = fs::path(data_path, data_ref, paste("thres",thres,sep="_"), "results", "SG")
sim_SG_path = fs::path(save_path, "results", "SG")
if(!dir.exists(sim_SG_path)){
  dir.create(sim_SG_path, recursive = T)
}
#SG_bin_test = getSG(raw_mat = new_raw,
#                    spat_mat = spatial,
#                    annot_df = annotation,
#                    wd = sim_SG_path,
#                    python_path = python_path, is_docker = F, spat_gene_method = "binspect",
#                    fname =  paste(fname_attach, "SG_bin", sep = "_"), 
#                    exp_thres = 1, min_cells = 1, min_genes = 1, max_dist_dly = "auto", pb = T)
#SG_bin_GO = SG_bin_test$GO
#SG_bin = SG_bin_test$SG
#spat_genes_bin = df2_yes$gene_name
#save(SG_bin_GO, SG_bin, spat_genes_bin, file = fs::path(sim_SG_path, paste(fname_attach, "SG_bin", sep = "_"), ext="RData"))
fwrite(data.frame(GO@norm_expr[df2_yes$gene_name,]),
       file = fs::path(sim_SG_path, paste(fname_attach, "SG_bin_norm", sep = "_"), ext = "csv"), sep=",", row.names = T)
fwrite(data.frame(as.matrix(GO@raw_exprs[df2_yes$gene_name,])),
       file = fs::path(sim_SG_path, paste(fname_attach, "SG_bin_raw", sep = "_"), ext = "csv"), sep=",", row.names = T)

##### save HVG
#n_genes = length(df2_yes$gene_name)
#sim_HVG_path = fs::path(data_path, data_ref, paste("thres",thres,sep="_"), "results", "HVG")
sim_HVG_path = fs::path(save_path, "results", "HVG")
if(!dir.exists(sim_HVG_path)){
  dir.create(sim_HVG_path, recursive = T)
}
HVG = genes_HVG_only
#save(HVG_GO, HVG, file = fs::path(sim_HVG_path, paste(fname_attach, "HVG", sep = "_"), ext = "RData"))
fwrite(data.frame(HVG_res$GO@norm_expr[HVG,]),
       file = fs::path(sim_HVG_path, paste(fname_attach, "HVG_norm", sep = "_"), ext = "csv"), sep=",", row.names = T)
fwrite(data.frame(as.matrix(HVG_res$GO@raw_exprs[HVG,])),
       file = fs::path(sim_HVG_path, paste(fname_attach, "HVG_raw", sep = "_"), ext = "csv"), sep=",", row.names = T)

####### save SG
#source("functions/getSG_01182021.R")
#sim_SG_path = fs::path(data_path, data_ref, paste("thres",thres,sep="_"), "results", "SG")
#if(!dir.exists(sim_SG_path)){
#  dir.create(sim_SG_path, recursive = T)
#}
#SG_bin_test = getSG(raw_mat = new_raw,
#                    spat_mat = spatial,
#                    annot_df = annotation,
#                    wd = sim_SG_path,
#                    python_path = python_path, is_docker = F, spat_gene_method = "binspect",
#                        fname =  paste(fname_attach, "SG_bin", sep = "_"), min_cells = 1, min_genes = 1, max_dist_dly = "auto", pb = T)
#SG_bin_GO = SG_bin_test$GO
#SG_bin = SG_bin_test$SG
#spat_genes_bin = df2_yes$gene_name
#save(SG_bin_GO, SG_bin, spat_genes_bin, file = fs::path(sim_SG_path, paste(fname_attach, "SG_bin", sep = "_"), ext="RData"))
#fwrite(data.frame(GO@norm_expr[df2_yes$gene_name,]),
#       file = fs::path(sim_SG_path, paste(fname_attach, "SG_bin_norm", sep = "_"), ext = "csv"), sep=",", row.names = T)
#fwrite(data.frame(as.matrix(GO@raw_exprs[df2_yes$gene_name,])),
#       file = fs::path(sim_SG_path, paste(fname_attach, "SG_bin_raw", sep = "_"), ext = "csv"), sep=",", row.names = T)
