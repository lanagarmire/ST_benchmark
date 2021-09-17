########## function that gets spatially variable genes

getSG = function(raw_mat, spat_mat, annot_df, wd, fname,top_spat_genes=6,seed_number=12345,pb = F,python_path, is_docker,#input 
                 exp_thres = 1, min_cells = 1, min_genes =1, #filter
                 sf = 6000, norm_m = "standard", ls_norm = T,  lognorm = T, scale_G = T, scale_C = T, #normalization
                 min_k = 0, max_dist_dly = "auto", method_sn = "Delaunay", method_dly = "deldir", #spatial network
                 Y=T, j=T, S=0, k=4, max_dist_knn = NULL, 
                 spat_gene_method #spatial gene exreaction method
                 #method_bin = "kmeans", ns = 3, max_iter = 10, perc_rk = 30, ft=T, #spatial genes
                 #hub=F, min_hub=3
                 ){
  library(Giotto)
  library(data.table)
  library(Matrix)
  #create Giotto object
  myinst = createGiottoInstructions(python_path = python_path, is_docker = is_docker)
  spat_test = createGiottoObject(raw_exprs = raw_mat, 
                                 spatial_locs = spat_mat, instructions = myinst)
  spat_test@cell_metadata$group = annot_df$group

  #filter
  spat_test = filterGiotto(gobject = spat_test, 
                           expression_threshold = exp_thres,
                           gene_det_in_min_cells = min_cells, ##filter genes
                           min_det_genes_per_cell = min_genes,##filter cells
                           expression_values = c('raw'), 
                           verbose = pb)
  #normalize
  spat_test = normalizeGiotto(gobject = spat_test, 
                              scalefactor = sf,
                              norm_methods = norm_m, 
                              library_size_norm = ls_norm,
                              log_norm = lognorm, 
                              scale_genes = scale_G,
                              scale_cells = scale_C, 
                              verbose = pb)
  spat_test = addStatistics(gobject = spat_test)
  #spat_test@cell_metadata$group = annot_df$group
  #create spatial network
  if(spat_gene_method == "binspect"){
    spat_test = createSpatialNetwork(gobject = spat_test, 
                                     method = method_sn, 
                                     delaunay_method = method_dly,
                                     minimum_k = min_k,
                                     maximum_distance_delaunay = max_dist_dly,
                                     Y=Y, j=j, S=S,  
                                     k=k,
                                     maximum_distance_knn=max_dist_knn,
                                     verbose=pb)
    Giotto::spatPlot(spat_test, show_network=T, network_color = "blue", cell_color = "group",
                     cell_color_code = Giotto:::getDistinctColors(length(table(annot_df$group))),
                     point_size=2.5,save_param=list(save_dir = wd,
                                                    save_name = paste0(fname,"_",spat_gene_method,"_spat_network"),
                                                    save_format = "png"),
                     save_plot=TRUE, title="spatial network",
                     legend_symbol_size = 2.5)
    #get spatial genes
    spatialgenes = binSpect(gobject = spat_test, expression_values = "normalized", set.seed = seed_number)
    spatGenePlot(spat_test, expression_values = "scaled", genes = spatialgenes[1:top_spat_genes]$genes,
                 point_shape="border", point_border_stroke = 0.1, show_network = F, 
                 network_color = 'lightgrey', point_size = 2.5, cow_n_col = 3,
                 save_param=list(save_dir = wd,
                                 save_name = paste0(fname,"_",spat_gene_method,"_top_sg"),
                                 base_width = 15,
                                 save_format = "png"),
                 save_plot=TRUE)
  }else if(spat_gene_method=="silhouetterank"){
    set.seed(seed_number)
    spatialgenes = silhouetteRank(gobject = spat_test, expression_values = "normalized")
    spatGenePlot(spat_test, expression_values = "scaled", genes = spatialgenes[1:top_spat_genes]$genes,
                 point_shape="border", point_border_stroke = 0.1, show_network = F, 
                 network_color = 'lightgrey', point_size = 2.5, cow_n_col = 3,
                 save_param=list(save_dir = wd,
                                 save_name = paste0(fname,"_",spat_gene_method,"_top_sg"),
                                 base_width = 15,
                                 save_format = "png"),
                 save_plot=TRUE)
  }else if(spat_gene_method=="SPARK"){
    set.seed(seed_number)
    spatialgenes = spark(gobject = spat_test) #raw
    spatGenePlot(spat_test, expression_values = "scaled", genes = spatialgenes[1:top_spat_genes]$genes,
                 point_shape="border", point_border_stroke = 0.1, show_network = F, 
                 network_color = 'lightgrey', point_size = 2.5, cow_n_col = 3,
                 save_param=list(save_dir = wd,
                                 save_name = paste0(fname,"_",spat_gene_method,"_top_sg"),
                                 base_width = 15,
                                 save_format = "png"),
                 save_plot=TRUE)
  }else if(spat_gene_method=="trendsceek"){
    set.seed(seed_number)
    spatialgenes = trendSceek(gobject = spat_test, expression_values = "normalized") 
    spatGenePlot(spat_test, expression_values = "scaled", genes = spatialgenes[1:top_spat_genes]$genes,
                 point_shape="border", point_border_stroke = 0.1, show_network = F, 
                 network_color = 'lightgrey', point_size = 2.5, cow_n_col = 3,
                 save_param=list(save_dir = wd,
                                 save_name = paste0(fname,"_",spat_gene_method,"_top_sg"),
                                 base_width = 15,
                                 save_format = "png"),
                 save_plot=TRUE)
  }else if(spat_gene_method=="spatialde"){
    set.seed(seed_number)
    spatialgenes = spatialDE(gobject = spat_test, expression_values = "raw",
                             save_plot = T, save_param = list(save_dir = wd,
                                                              save_name = paste0(fname,"_spatialDE"),
                                                              base_width = 15,
                                                              save_format = "png"))
  }
  #save the information
  #setwd(wd)
  SG_res = list(GO=spat_test, SG=spatialgenes)
  #save(SG_res, file = paste0("SG_res_",seed_number,".RData"))
  return(SG_res)
}
