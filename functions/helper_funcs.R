##########################
#### some helper functions
##########################

############
### extract variance percentages
############
compute_var_pc = function(gobject,name,reduction=c('cells', 'genes')){
  reduction = match.arg(reduction, c('cells', 'genes'))
  pca_obj = gobject@dimension_reduction[[reduction]]$pca[[name]]
  eigs = pca_obj$misc$eigenvalues
  var_expl = eigs/sum(eigs)*100
  var_expl_cum = cumsum(eigs)/sum(eigs)*100
  return(list(var_PC = var_expl, var_PC_cum = var_expl_cum))
}

PC_dist_heat = function(PC_dist,cls_order,cls_col,name){
  #PC_dist = as.matrix(dist(PC_mat[,1:nPC,drop=FALSE]))
  PC_dist_reordered = PC_dist[unlist(cls_order),unlist(cls_order)]
  fig = heatmap.2(PC_dist_reordered,Rowv=FALSE,Colv=FALSE,trace="none",dendrogram="none",main=name,col=colorRampPalette(brewer.pal(9,"Reds"))(30),
                  RowSideColors=unname(unlist(mapply(FUN=function(x,y){rep(x,y)},cls_col,sapply(cls_order,FUN=function(x){length(x)})))),
                  ColSideColors=unname(unlist(mapply(FUN=function(x,y){rep(x,y)},cls_col,sapply(cls_order,FUN=function(x){length(x)})))))
  #return(fig)
}

sNN_heat = function(gobject,N,cls_order,cls_col,pca_name,nPC,num_nn,name){
  gobject = createNearestNetwork(gobject,type="sNN",dim_reduction_to_use="pca",dim_reduction_name=pca_name,dimensions_to_use=1:nPC,k=num_nn,name=name)
  sNN=gobject@nn_network$sNN[[name]]$igraph
  wts_mat = as.matrix(sNN[1:N,1:N])
  wts_mat_reordered = wts_mat[unlist(cls_order),unlist(cls_order)]
  fig = heatmap.2(wts_mat_reordered,Rowv=FALSE,Colv=FALSE,trace="none",dendrogram="none",main=name, col=colorRampPalette(brewer.pal(9,"Reds"))(30),                 
                  RowSideColors=unname(unlist(mapply(FUN=function(x,y){rep(x,y)},cls_col,sapply(cls_order,FUN=function(x){length(x)})))),
                  ColSideColors=unname(unlist(mapply(FUN=function(x,y){rep(x,y)},cls_col,sapply(cls_order,FUN=function(x){length(x)})))))
}

compute_purity=function(gobject,ref,lab,matched){
  conf_mat = table(gobject@cell_metadata[[ref]],gobject@cell_metadata[[lab]])
  if(matched){
    purity_clus=diag(conf_mat)/colSums(conf_mat)
    purity_overall=sum(diag(conf_mat))/sum(conf_mat)
  }else{
    purity_clus=apply(conf_mat,2,FUN=function(x){max(x)/sum(x)})
    purity_overall=purity_clus*colSums(conf_mat)/sum(conf_mat)
  }
  purity_clus = data.frame(cluster=colnames(conf_mat),purity=purity_clus)
  return(list(confusion_matrix=conf_mat,purity=purity_clus,overall_purity=purity_overall))
}

compute_cluster_specific_F1 = function(gobject,ref,lab){
  conf_mat = table(gobject@cell_metadata[[ref]],gobject@cell_metadata[[lab]])
  n_clus = nrow(conf_mat)
  cluster_acc = lapply(1:n_clus,FUN=function(i){TP = conf_mat[i,i]; TN=sum(conf_mat[-i,-i]); FP=sum(conf_mat[i,-i]); FN=sum(conf_mat[-i,i]); precision=TP/(TP+FP); recall=TP/(TP+FN); F1=ifelse(precision+recall!=0,2*precision*recall/(precision+recall),0); return(c(TP,FP,TN,FN,precision,recall,F1))})
  cluster_acc_df = cbind(colnames(conf_mat),as.data.frame(Reduce(rbind,cluster_acc)))
  colnames(cluster_acc_df) = c("cluster","TP","FP","TN","FN","precision","recall","F1")
  weighted_F1 = sum(cluster_acc_df$F1*(cluster_acc_df$TP+cluster_acc_df$FN))/sum(conf_mat)
  return(list(confusion_matrix=conf_mat,F1=cluster_acc_df,weighted_F1=weighted_F1))
}

add_trailing_zero = function(x){
  if(x==as.integer(x)){
    ret = paste0(x,".0")
  }else{
    ret = as.character(x)
  }
  return(ret)
}

jaccard = function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

stablize_sct = function(raw_exprs,method="poisson",min_cells=1){
  n_cells = ncol(raw_exprs)
  clip.range = c(-sqrt(n_cells/30),sqrt(n_cells/30))
  cell.attr.df = data.frame(umi = colSums(raw_exprs), log_umi=log10(colSums(raw_exprs)))
  vst.out = vst(raw_exprs,,method=method,vst.flavor="v2",return_corrected_umi=T,return_gene_attr=T,return_cell_attr=T,cell_attr=cell.attr.df,n_cells=min(5000,ncol(raw_exprs)),min_cells=min_cells)
  scale.data = vst.out$y
  # clip the residuals
  scale.data[scale.data < clip.range[1]] = clip.range[1]
  scale.data[scale.data > clip.range[2]] = clip.range[2]
  # center
  scale.data = t(apply(scale.data,1,FUN=function(x){x-mean(x)}))
  vst.out$gene_attr$gene_ID=rownames(vst.out$gene_attr)
  vst.out$y = scale.data
  return(vst.out)
}

affinity_to_snn = function(aff_matrix,num_nn,rank_threshold=3,shared_threshold=5){
  ## convert to distance matrix
  dist_matrix = aff_matrix
  dist_matrix = (dist_matrix-min(dist_matrix))/diff(range(dist_matrix))
  diat_matrix = 1/(1+dist_matrix)
  diag(dist_matrix)=0
  
  ## make sNN
  cell_names = rownames(aff_matrix)
  names(cell_names) = 1:nrow(aff_matrix)
  nn_network = dbscan::kNN(x = as.dist(dist_matrix), k = num_nn, sort = TRUE)
  from = to = weight = distance = from_cell_ID = to_cell_ID = shared = NULL
  nn_network_dt = data.table::data.table(from = rep(1:nrow(nn_network$id), k=num_nn),
                                         to = as.vector(nn_network$id),
                                         weight = 1/(1 + as.vector(nn_network$dist)),
                                         distance = as.vector(nn_network$dist))
  nn_network_dt[, from_cell_ID := cell_names[from]]
  nn_network_dt[, to_cell_ID := cell_names[to]]
  snn_network = dbscan::sNN(x = nn_network, k = num_nn, kt = NULL)
  snn_network_dt = data.table::data.table(from = rep(1:nrow(snn_network$id), k=num_nn),
                                          to = as.vector(snn_network$id),
                                          weight = 1/(1 + as.vector(snn_network$dist)),
                                          distance = as.vector(snn_network$dist),
                                          shared = as.vector(snn_network$shared))
  snn_network_dt = snn_network_dt[stats::complete.cases(snn_network_dt)]
  snn_network_dt[, from_cell_ID := cell_names[from]]
  snn_network_dt[, to_cell_ID := cell_names[to]]
  data.table::setorder(snn_network_dt, from, -shared)
  snn_network_dt[, rank := 1:.N, by = from]
  snn_network_dt = snn_network_dt[rank <= rank_threshold | shared >= shared_threshold]
  all_index = unique(x = c(nn_network_dt$from_cell_ID, nn_network_dt$to_cell_ID))
  missing_indices = all_index[!all_index %in% unique(snn_network_dt$from)]
  nn_network_igraph = igraph::graph_from_data_frame(snn_network_dt[,.(from_cell_ID, to_cell_ID, weight, distance, shared, rank)], directed = TRUE, vertices = all_index)
  
  ##
  return(list(nn=nn_network_igraph,dist=dist_matrix))
}




subsetGiotto_v2 = function(gobject,cell_ids=NULL,gene_ids=NULL){
  if(!is.null(gobject@custom_expr)){
    if(is.null(cell_ids)){
      cell_ids = colnames(gobject@raw_exprs)
    }
    if(is.null(gene_ids)){
      gene_ids = rownames(gobject@raw_exprs)
    }
    old_custom_expr = gobject@custom_expr
    gobject@custom_expr=NULL
    gobject = subsetGiotto(gobject,cell_ids=cell_ids,gene_ids=gene_ids)
    gobject@custom_expr=old_custom_expr[rownames(old_custom_expr)%in%gene_ids,colnames(old_custom_expr)%in%cell_ids]
  }else{
    gobject=subsetGiotto(gobject,cell_ids=cell_ids,gene_ids=gene_ids)
  }
  return(gobject)
}

