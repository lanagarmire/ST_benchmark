library(Giotto)
library(dplyr)
library(data.table)

setwd("/home/liyijun/ST_benchmark_01082020/")

### read external parameters
args = commandArgs(trailingOnly = TRUE)
data_path = as.character(args[1])
data_ref = as.character(args[2])

### set reference dataset
library(data.table)
library(Giotto)
load(fs::path(data_path, "real_datasets", data_ref, ext="RData"))

cell_id = colnames(GO@raw_exprs)
grp_id = sort(unique(annotation$group))
gene_id = rownames(GO@raw_exprs)
expr_data = GO@norm_expr

gene_grp_id = rep(NA, nrow(GO@raw_exprs))
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
