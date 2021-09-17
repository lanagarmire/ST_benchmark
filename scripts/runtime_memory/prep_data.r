#################### cells ##############################
rm(list=ls())
library(data.table)

load("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/Visium/data/Visium_HCere.RData")

raw_expr = as.matrix(raw_expr)
save_path = "/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/runtime_cells/data"
num_cells_vec = c(150, 300, 500, 1000, 1500, 3000, 4500)

seed_list = 12345 + c(1:5) - 1
n_cells = ncol(raw_expr)
set.seed(12345)
genes_idx1 = sample(nrow(raw_expr), size=500)
genes_idx2 = sample(nrow(raw_expr), size=500)

#save varying cell numbers
for(i in 1:length(num_cells_vec)){
  for(j in 1:length(seed_list)){
    set.seed(seed_list[j])
    #cells_idx1 = sample(n_cells, size=num_cells_vec[i])
    #cells_idx2 = sample(n_cells, size=num_cells_vec[i])
    cells_idx = sample(n_cells, size=num_cells_vec[i])
 
    raw_data_1 = raw_expr[genes_idx1, cells_idx]
    raw_data_2 = raw_expr[genes_idx2, cells_idx]

    norm_data_1 = norm_expr[genes_idx1, cells_idx]
    norm_data_2 = norm_expr[genes_idx2, cells_idx]

    fwrite(data.frame(raw_data_1), file = fs::path(save_path, paste("data", "cells", num_cells_vec[i], "raw", 1, seed_list[j], sep="_"), ext="csv"), row.names=T)
    fwrite(data.frame(raw_data_2), file = fs::path(save_path, paste("data", "cells", num_cells_vec[i], "raw", 2, seed_list[j], sep="_"), ext="csv"), row.names=T)

    fwrite(data.frame(norm_data_1), file = fs::path(save_path, paste("data", "cells", num_cells_vec[i], "norm", 1, seed_list[j], sep="_"), ext="csv"), row.names=T)
    fwrite(data.frame(norm_data_2), file = fs::path(save_path, paste("data", "cells", num_cells_vec[i], "norm", 2, seed_list[j], sep="_"), ext="csv"), row.names=T)
  }
}

#make the cell and feature index files
cells_df = expand.grid(num_cells_vec, seed_list)
colnames(cells_df) = c("cells", "reps")
fwrite(cells_df, file = fs::path(save_path, "cells", ext="csv"), row.names=T)



############################# features ############################
rm(list=ls())
library(data.table)

load("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/Visium/data/Visium_HCere.RData")

raw_expr = as.matrix(raw_expr)
save_path = "/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/runtime_features/data"
num_features_vec = c(150, 300, 500, 1000, 3000, 5000, 10000)

seed_list = 12345 + c(1:5) - 1
set.seed(12345)
cells_idx = sample(ncol(raw_expr), size=500)
n_features = nrow(raw_expr)
#save varying feature numbers
for(i in 1:length(num_features_vec)){
  for(j in 1:length(seed_list)){
    set.seed(seed_list[j])
    features_idx1 = sample(n_features, size=num_features_vec[i])
    features_idx2 = sample(n_features, size=num_features_vec[i])

    raw_data_1 = raw_expr[features_idx1, cells_idx]
    raw_data_2 = raw_expr[features_idx2, cells_idx]

    norm_data_1 = norm_expr[features_idx1, cells_idx]
    norm_data_2 = norm_expr[features_idx2, cells_idx]

    fwrite(data.frame(raw_data_1), file = fs::path(save_path, paste("data", "features", num_features_vec[i], "raw", 1, seed_list[j], sep="_"), ext="csv"), row.names=T)
    fwrite(data.frame(raw_data_2), file = fs::path(save_path, paste("data", "features", num_features_vec[i], "raw", 2, seed_list[j], sep="_"), ext="csv"), row.names=T)
 
    fwrite(data.frame(norm_data_1), file = fs::path(save_path, paste("data", "features", num_features_vec[i], "norm", 1, seed_list[j], sep="_"), ext="csv"), row.names=T)
    fwrite(data.frame(norm_data_2), file = fs::path(save_path, paste("data", "features", num_features_vec[i], "norm", 2, seed_list[j], sep="_"), ext="csv"), row.names=T)
  }
}

features_df = expand.grid(num_features_vec, seed_list)
colnames(features_df) = c("features", "reps")
#features_df = data.frame(features = num_features_vec)

fwrite(features_df, file = fs::path(save_path, "features", ext="csv"), row.names=T)
#fwrite(ndims_df, file = fs::path(save_path, "ndims", ext="csv"), row.names=T)
#fwrite(num_nn_df, file = fs::path(save_path, "num_nn", ext="csv"), row.names=T)

