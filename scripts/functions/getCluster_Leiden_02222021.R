##expr matrix: gene x cells
getCluster_Leiden = function(gene_expr, spat_mat, pb, data_type, pc_dims, py_path,reduced,
                             center = F, su = F, pca_method = "factominer", pcs=100,#PCA
                             num_neighbors=30, min_shared = 5, top_shared=3, #shared Nearest Network
                             ncl_truth = 0, res = 0.8, res_step = 0.4, seed = 1234, #Leiden clustering
                             max_iter = 50, num_iter=1000, partition_type="RBConfigurationVertexPartition"
){
  library(Giotto)
  #create Giotto object
  instrs = createGiottoInstructions(python_path = py_path)
  GO_test = createGiottoObject(raw_exprs = gene_expr, norm_expr = gene_expr, spatial_locs = spat_mat, 
                               instructions = instrs)
  if(reduced==TRUE){
    #integration tool includes dimension reduction
    GO_test@dimension_reduction$cells$pca$pca$coordinates = t(gene_expr)
  }else{
    #do PCA
    GO_test = runPCA(gobject = GO_test, 
                     expression_values="normalized",
                     reduction = "cells", 
                     genes_to_use=NULL, 
                     center = center, 
                     scale_unit = su, 
                     ncp = pcs,
                     method = pca_method,
                     verbose = pb)
  }
  #GO_test = createNearestNetwork(gobject=GO_test, 
  #                               type = "sNN", 
  #                               dim_reduction_to_use = "pca",
  #                               dim_reduction_name = "pca", 
  #                               dimensions_to_use = 1:pc_dims,
  #                               expression_values = "normalized", 
  #                               k = num_neighbors,
  #                               minimum_shared = min_shared, 
  #                               top_shared = top_shared, 
  #                               verbose = pb)
  if(ncl_truth != 0){
    res_init = res
    num_cl = 0
    iter0 = 1
    while(num_cl != ncl_truth & num_neighbors > 0 & iter0 < 15){
      if(num_cl != 0){
        if(num_cl > ncl_truth){
          if(num_neighbors <= 10){
            num_neighbors = num_neighbors + 1
          }else{
            num_neighbors = num_neighbors + 5
          }
        }else if(num_cl < ncl_truth){
          if(num_neighbors <= 10){
            num_neighbors = num_neighbors - 1
          }else{
            num_neighbors = num_neighbors - 5
          }
        }
      }
      GO_test = createNearestNetwork(gobject=GO_test,
                                 type = "sNN",
                                 dim_reduction_to_use = "pca",
                                 dim_reduction_name = "pca",
                                 dimensions_to_use = 1:pc_dims,
                                 expression_values = "normalized",
                                 k = num_neighbors,
                                 minimum_shared = min_shared,
                                 top_shared = top_shared,
                                 verbose = pb)

      iter0 = iter0 + 1
      iter = 1
      res = res_init
      while(num_cl != ncl_truth & iter < max_iter & res <= 1 & res > 9.1e-4){
        GO_test = doLeidenCluster(gobject = GO_test, 
                                resolution = res,
                                name = paste0("LC_",data_type,"_clus"),
                                n_iterations = num_iter,
                                partition_type = partition_type,
                                python_path = py_path,
                                set_seed = T,
                                seed_number = seed)
        num_cl = length(table(GO_test@cell_metadata[,2]))
        if(num_cl > ncl_truth){
          res = res - res_step/iter
        }else if(num_cl < ncl_truth){
          res = res + res_step/iter
        }
      #print(paste0(num_cl, " clusters; ", iter, " iterations; ", res, " resolution"))
        print(res)
        iter = iter + 1
      }
      print(num_neighbors)
    }  
  }else{
    GO_test = doLeidenCluster(gobject = GO_test, 
                              resolution = res,
                              name = paste0("LC_",data_type,"_clus"),
                              n_iterations = num_iter,
                              partition_type = partition_type,
                              python_path = py_path)
  }
  
  #save results
  label_df = as.data.frame(GO_test@cell_metadata)
  GO_ls = list(GO=GO_test, cluster_label = label_df, num_nn = num_neighbors, res = res)
  return(GO_ls)
}
