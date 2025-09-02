getCluster_SC3 = function(count_expr,norm_expr,celltype_df,cl0,SVM,data_type,d_region_min,d_region_max,
                          su=F,pcs,filter=FALSE,bio=FALSE,StepByStep=FALSE){
  library(SingleCellExperiment)
  library(SC3)
  library(scater)
  
  # remove cells with no variability
  non_zero_variance_cells = names(which(apply(norm_expr, 2, sd) != 0))
  print(paste0("non zero variance cells: ",length(non_zero_variance_cells)))
  
  count_expr = count_expr[,colnames(count_expr) %in% non_zero_variance_cells,drop=FALSE]
  norm_expr = norm_expr[,colnames(norm_expr) %in% non_zero_variance_cells,drop=FALSE]
  
  # create a SingleCellExperiment object
  sce = SingleCellExperiment(assays = list(counts = as.matrix(count_expr),logcounts = as.matrix(norm_expr)))
  rowData(sce)$feature_symbol = rownames(sce)
  if(!StepByStep){
    if(!SVM){
      if(length(non_zero_variance_cells)>2000){
        sce = sc3(sce, ks=cl0, gene_filter=filter, biology=bio, kmeans_nstart=50, svm_max=ncol(count_expr)+1000, d_region_min=d_region_min, d_region_max=d_region_max) #set the maximum cells to larger than the actual number of cells, to avoid triggering SVM
      }else{
        sce = sc3(sce, ks=cl0, gene_filter=filter, biology=bio, kmeans_nstart=100, svm_max=ncol(count_expr)+1000, d_region_min=d_region_min, d_region_max=d_region_max)
      }    
    }else{
      sce = sc3(sce, ks=cl0, gene_filter=filter, biology=bio,svm_max=10, svm_num_cells = round(ncol(count_expr)*0.2)) #set the maximum cells to a very small number so that SVM will be triggered
      #sce = sc3_run_svm(sce, ks = cl0)
    }
    individual_clus_res = NA
  }else{
    if(!SVM){
      if(length(non_zero_variance_cells)>2000){
        sce = sc3_prepare(sce, gene_filter=filter,kmeans_nstart=50, svm_max=ncol(count_expr)+1000, d_region_min=d_region_min, d_region_max=d_region_max) #set the maximum cells to larger than the actual number of cells, to avoid triggering SVM
      }else{
        sce = sc3_prepare(sce, gene_filter=filter,kmeans_nstart=100, svm_max=ncol(count_expr)+1000, d_region_min=d_region_min, d_region_max=d_region_max)
      }    
    }else{
      sce = sc3_prepare(sce, gene_filter=filter, svm_max=10, svm_num_cells = round(ncol(count_expr)*0.2)) #set the maximum cells to a very small number so that SVM will be triggered
      #sce = sc3_run_svm(sce, ks = cl0)
    }
    
    sce = sc3_calc_dists(sce)
    sce = sc3_calc_transfs(sce)
    sce = sc3_kmeans(sce, ks = cl0)
    individual_clus_res = metadata(sce)$sc3$kmeans
    sce = sc3_calc_consens(sce)
  }
  
  col_data = colData(sce)
  col_data = as.data.frame(col_data)
  col_data$cell_ID=rownames(col_data) 
  label_df = col_data
  colnames(label_df)[which(colnames(label_df)==paste("sc3",cl0,"clusters",sep="_"))]=paste("SC3",data_type,"clus",sep="_")

  ## colData(sce) = cbind(colData(sce),celltype_df%>%dplyr::select(-cell_ID))
  ## col_data = colData(sce)
  ## col_data = as.data.frame(col_data)
  ## label_df = data.frame(cell_ID = rownames(col_data), clus = col_data[, paste("sc3",cl0,"clusters",sep="_")])
  ## colnames(label_df)[2] = paste("SC3",data_type,"clus",sep="_")
  
  res_ls = list(SC3_obj = sce, cluster_label = label_df, individual_result = individual_clus_res, transf_ls = metadata(sce)$sc3$transformations)
  return(res_ls)
}