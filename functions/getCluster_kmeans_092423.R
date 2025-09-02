getCluster_kmeans = function(gene_expr, data_type, reduced, PCA_mat = NULL, pc_dims, dist = "euclidean", batch_size, batch,
                             su = F, pcs=100, cl0=0, iter_max=100, nstart=100, seed=1234){ #Kmeans clustering
  library(proxy)
  library(FactoMineR)
  library(Rfast)
  library(ClusterR)
  
  non_zero_variance_cells = names(which(apply(gene_expr, 2, sd) != 0))
  gene_expr = gene_expr[,colnames(gene_expr) %in% non_zero_variance_cells,drop=FALSE] ## these two lines ignorable when reduced = T
  
  ## perform dimension reduction
  if(reduced==TRUE){
    #PCA_df = t(gene_expr)
    PCA_df = PCA_mat[,1:pc_dims]
  }else{
    pca_res = PCA(t(gene_expr),scale.unit=su,ncp=pcs,graph=F)
    eigenvalues = pca_res$eig[, 1]
    loadings = sweep(pca_res$var$coord, 2, sqrt(eigenvalues[1:min(pcs,ncol(pca_res$var$coord))]),FUN = "/")
    rownames(loadings) = rownames(gene_expr)
    colnames(loadings) = paste0("Dim.", 1:ncol(loadings))
    coords = pca_res$ind$coord
    rownames(coords) = colnames(gene_expr)
    colnames(coords) = paste0("Dim.", 1:ncol(coords))
    
    PCA_df = coords[,1:pc_dims]
  }
  
  ## compute distance metrics
  if(dist == "euclidean"){
    celldist = proxy::dist(PCA_df,method=dist,parallel=4)
  }else if(dist == "pearson"){
    celldist = stats::as.dist(1-cora(t(PCA_df)))
  }else if(dist == "spearman"){
    PCA_df_ranked = t(apply(PCA_df,1,rank))
    celldist = stats::as.dist(1-cora(t(PCA_df_ranked)))
  }
  
  #do kmeans clustering 
  kmeans_opt = cl0
  set.seed(seed)

  if(!batch){
    kclusters = kmeans(x=celldist,centers=kmeans_opt,iter.max=iter_max,nstart=nstart)
    label_df = data.table::data.table(cell_ID = names(kclusters[["cluster"]]),
                                               name = kclusters[["cluster"]])
  }else{
    km_mini = MiniBatchKmeans(as.matrix(celldist),kmeans_opt,batch_size=round(batch_size*length(non_zero_variance_cells)))
    km_pred = predict_MBatchKMeans(as.matrix(celldist),km_mini$centroids)
    label_df = data.table::data.table(cell_ID = rownames(PCA_df), name = c(km_pred))
  }
  data.table::setnames(label_df, "name", paste0("KM_",data_type,"_",dist,"_clus"))
  label_df = as.data.frame(label_df)
  
  #save results
  if(!batch){
    res_ls = list(kmeans_obj = kclusters, celldist = celldist, cluster_label = label_df, kmeans_opt = kmeans_opt)
  }else{
    res_ls = list(kmeans_obj = km_mini, celldist = celldist, cluster_label = label_df, kmeans_opt = kmeans_opt)
  }
    
  return(res_ls)
}