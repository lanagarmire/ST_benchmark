run_cimlr = function(hvg_data, svg_data, ncl_truth, no.dim=NA, k=10, cores.ratio=1, sd = 12345){#normalized 
  library(CIMLR)
  
  dat_list = list(hvg_data, svg_data)

   set.seed(sd)
  res = CIMLR(dat_list, ncl_truth, no.dim = no.dim, k=k, cores.ratio = cores.ratio)
  
  cimlr_labels = res$y$cluster
  cell_ID = colnames(hvg_data)
  
  similarity_mat = res$S
  
  return(list(cluster_label = data.frame(cell_ID = cell_ID, CIMLR_clus =  cimlr_labels),
              similarity_matrix = similarity_mat))
}