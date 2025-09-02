getCluster_cellTree = function(gene_expr, ndims=NULL, cl0, data_type, method = "maptpx", filter=FALSE, su = FALSE, use_prob=TRUE){
  library(cellTree)
  
  ## remove non variance cells
  non_zero_variance_cells = names(which(apply(gene_expr, 2, sd) != 0))
  gene_expr = gene_expr[,colnames(gene_expr) %in% non_zero_variance_cells,drop=FALSE]
  
  ## run LDA
  if(use_prob){
    lda_res = compute.lda(gene_expr, method = method, log.scale = su, sd.filter = filter)
    celldist = get.cell.dist(lda_res)
    
    max_prob = apply(lda_res$omega,1,FUN=function(x){which.max(x)})
    label_df = data.frame(cell_ID = names(max_prob), clus = as.numeric(unname(max_prob)))
    colnames(label_df)[2] = paste("cellTree",data_type,"clus",sep="_")
  }else{
    lda_res = compute.lda(gene_expr, method = method, log.scale = su, sd.filter = filter)
    celldist = get.cell.dists(lda_res)
    HC = hclust(d=as.dist(celldist))
    memb=cutree(HC,k=cl0)
    
    label_df = data.frame(cell_ID = names(memb), clus = as.numeric(unname(memb)))
    colnames(label_df)[2] = paste("cellTree",data_type,"clus",sep="_")
  }
  
  
  res_ls = list(cellTree_obj = lda_res, cluster_label = label_df, celldist=celldist)
  return(res_ls)
}