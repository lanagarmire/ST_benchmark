######## function that gets the HVG
### raw_mat:genes x cell, spat_mat: cell x 2

getHVG = function(raw_mat, spat_mat, wd, fname, python_path, is_docker,#input 
                  exp_thres = 1, min_cells = 10, min_genes =10, #filter
                  sf = 6000, norm_m = "standard", ls_norm = T,  lognorm = T, scale_G = T, scale_C = T, #normalization
                  HVG_method = "cov_loess", diff_cov = 0.1, num_exp_grp = 20, z_thres = 1.5,  #hvg
                  #perc_cell_thres = 4, mean_exp_thres = 0.5, #extract HVG
                  pb = F 
                  ){
  library(Giotto)
  library(data.table)
  library(Matrix)
  #create Giotto object
  myinst = createGiottoInstructions(python_path = python_path, is_docker = is_docker)
  HVG_test = createGiottoObject(raw_exprs = raw_mat, 
                                spatial_locs = spat_mat, instructions = myinst)
  #filter
  HVG_test = filterGiotto(gobject = HVG_test, 
                          expression_threshold = exp_thres,
                          gene_det_in_min_cells = min_cells, ##filter genes 
                          min_det_genes_per_cell = min_genes,##filter cells
                          expression_values = c('raw'), 
                          verbose = pb)
  #normalize
  HVG_test = normalizeGiotto(gobject = HVG_test, 
                             scalefactor = sf,
                             norm_methods = norm_m, 
                             library_size_norm = ls_norm,
                             log_norm = lognorm, 
                             scale_genes = scale_G,
                             scale_cells = scale_C, 
                             verbose = pb)
  HVG_test = addStatistics(gobject = HVG_test)
  #get HVG
  HVG_test = calculateHVG(gobject = HVG_test,
                          method = HVG_method,
                          difference_in_cov = diff_cov,
                          expression_values = c('normalized'), 
                          nr_expression_groups = num_exp_grp, 
                          zscore_threshold = z_thres,
                          save_plot = T,
                          save_param = list(save_dir = wd, save_name = paste0(fname,"_HVG"), save_format = "png"))
  gene_metadata = fDataDT(HVG_test)
  #featgenes = gene_metadata[hvg == 'yes' & perc_cells > perc_cell_thres
  #                          & mean_expr_det > mean_exp_thres]$gene_ID
  featgenes = gene_metadata[hvg == 'yes']$gene_ID
  #save the information
  #setwd(wd)
  HVG_res = list(GO=HVG_test, HVG=featgenes)
  #save(HVG_res, file = "HVG_res.RData")
  return(HVG_res)
}
