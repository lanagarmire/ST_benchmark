#######################################
#### this script does preprocessing for downloaded datasets
#######################################

convertFactor = function(x) {
  # Check if the factor can be converted to numeric without warnings
  suppressWarnings(
    if(is.factor(x)){
      if (all(!is.na(as.numeric(levels(x))))) {
        # Convert factor to numeric
        return(as.numeric(levels(x))[x])
      } else {
        # Convert factor to character
        return(as.character(x))
      }
    }else{
      return(x)
    }
  )
}

## add space for additional HVG and SVG method
## add space for random non-genes

benchmark_preprocess = function(raw_exprs, spatial_locs, cell_metadata, scale_factor, cell_type_col, python_path,
                                remove_rare_CT, rare_CT_threshold, cell_filter, 
                                gene_filter, 
                                pearson_residual, 
                                geneset_cutoff,
                                save_path, 
                                run_svg, SVG_path, SVG_method = "SPARK",  
                                run_hvg, HVG_path, 
                                union_path, data_name, QC_path, random_gs_reps=10, random_gs_seed=12345){
  dir.create(QC_path,recursive=TRUE,showWarnings = FALSE)
  ## remove useless genes
  gene_IDs = rownames(raw_exprs)
  
  ## check mitochondrial genes, blank genes, and negativeprobes
  mito_gs = gene_IDs[grep(pattern="(^MT-)|(^mt-)",gene_IDs)]
  print(paste("removing",length(mito_gs),"mitochondrial genes",sep=" "))
  blank_gs = gene_IDs[grep(pattern="(^Blank)|(blank)",gene_IDs)]
  print(paste("removing",length(blank_gs),"blank genes",sep=" "))
  negprb_gs = gene_IDs[grep(gene_IDs,pattern="(^NegPrb)")]
  print(paste("removing",length(negprb_gs),"negative probes",sep=" "))

  ## remove rare cell types
  if(cell_type_col!="cell_type"){
    colnames(cell_metadata)[which(colnames(cell_metadata)==cell_type_col)]="cell_type"
  }
  cell_metadata$cell_type = convertFactor(cell_metadata$cell_type)
  if(remove_rare_CT){
    CT_tab = table(cell_metadata$cell_type)/sum(table(cell_metadata$cell_type))
    rare_CT = names(CT_tab[which(CT_tab < rare_CT_threshold)])
    rare_CT_ID = cell_metadata %>% filter(cell_type %in% rare_CT) %>% pull(cell_ID)
    print(CT_tab)
  }else{
    rare_CT_ID=c()
  }
  
  ## filter cells
  if(!is.null(cell_filter)){
    cell_filter_ID = colnames(raw_exprs)[which(colSums(raw_exprs>0)/nrow(raw_exprs) < cell_filter)]
  }else{
    cell_filter_ID = c()
  }
  print(paste("removing",length(cell_filter_ID),"cells/spots by custom filter",sep=" "))
  
  ## create Giotto object
  cell_IDs = colnames(raw_exprs)
  rownames(spatial_locs)=spatial_locs$cell_ID
  instr = createGiottoInstructions(python_path=python_path)
  GO_obj = createGiottoObject(raw_exprs=raw_exprs[setdiff(gene_IDs,c(mito_gs,blank_gs,negprb_gs)),setdiff(cell_IDs,union(rare_CT_ID,cell_filter_ID))],spatial_locs=spatial_locs[setdiff(cell_IDs,union(rare_CT_ID,cell_filter_ID)),c("sdimx","sdimy")],instructions=instr)
  GO_obj = addCellMetadata(gobject = GO_obj, new_metadata = cell_metadata %>% filter(cell_ID %in% setdiff(cell_IDs,union(rare_CT_ID,cell_filter_ID))), by_column = T, column_cell_ID = "cell_ID")
  
  ## basic filtering
  GO_obj = filterGiotto(GO_obj,expression_threshold=1,gene_det_in_min_cells = 1, min_det_genes_per_cell = 1,verbose=TRUE)
  GO_obj = normalizeGiotto(gobject = GO_obj, scalefactor = scale_factor)
  GO_obj = addStatistics(GO_obj)
  
  ## further gene filtering
  if(gene_filter=="by_pca"){
    percentages = seq(0.01,0.1,by=0.01)
    N_spots = nrow(GO_obj@cell_metadata)
    top_PC_percentages = c()
    for(i in 1:length(percentages)){
      GO_obj_filter = GO_obj
      GO_obj_filter = filterGiotto(GO_obj_filter,gene_det_in_min_cells=ceiling(percentages[i]*N_spots)+1,min_det_genes_per_cell = 1)
      GO_obj_filter = normalizeGiotto(gobject=GO_obj_filter,scalefactor=scale_factor)
      GO_obj_filter = runPCA(GO_obj_filter,name="pca_filter",genes_to_use="all",center=TRUE,scale_unit=TRUE)
      PC_var_ls = compute_var_pc(GO_obj_filter,name="pca_filter")
      PC_elbow = PCAtools::findElbowPoint(PC_var_ls$var_PC)
      top_PC_percentages = c(top_PC_percentages, PC_var_ls$var_PC_cum[PC_elbow])
    }
    filter_thres_idx = which.max(diff(top_PC_percentages))+1
    gs_filter_thres = percentages[filter_thres_idx]
  }else if(!is.null(gene_filter)){
    gs_filter_thres = gene_filter
  }else{
    gs_filter_thres = NULL
  }
  if(!is.null(gs_filter_thres)){
    gs_filtered = rownames(GO_obj@raw_exprs)[which(rowSums(GO_obj@raw_exprs>0)/ncol(GO_obj) < gs_filter_thres)]
  }else{
    gs_filtered = c()  
  }
  print(paste("removing",length(gs_filtered),"genes by custom filter",sep=" "))
  
  ## save QC metrics 
  dir.create(fs::path(QC_path,data_name),recursive=TRUE,showWarnings = FALSE)
  png(fs::path(QC_path,data_name,"perc_feat",ext="png"))
  hist(colSums(GO_obj@raw_exprs>0)/nrow(GO_obj@raw_exprs),breaks="FD",main="percentage of features detected")
  dev.off()
  
  png(fs::path(QC_path,data_name,"perc_cells",ext="png"))
  hist(rowSums(GO_obj@raw_exprs>0)/ncol(GO_obj@raw_exprs),breaks="FD",main="percentage of cells detected")
  dev.off()
  
  print("save QC metrics")
  QC_mets = list(mito_gs=mito_gs, blank_gs=blank_gs, negprob_gs=negprb_gs, N_spots_post_filter=N_spots, top_PC_percentages=data.frame(perc=percentages,PC_perc=top_PC_percentages), filter_thres=gs_filter_thres)
  saveRDS(QC_mets,file = fs::path(QC_path,data_name,ext="RDS"))
  
  ## SCT adjustment
  if(pearson_residual){
    print("run SCT adjustment")
    stablize_sct_res = stablize_sct(GO_obj@raw_exprs)
    GO_obj@custom_expr = stablize_sct_res$y
    GO_obj@gene_metadata = merge(GO_obj@gene_metadata, stablize_sct_res$gene_attr,by="gene_ID",sort=FALSE)
  }
  
  ## run SVG
  if(run_svg){
    print("run SVG")
    dir.create(SVG_path,recursive=T,showWarnings = F)
    if(SVG_method == "SPARK"){
        tryCatch(expr = {
        SG_sp = Giotto::spark(gobject = GO_obj, percentage = gs_filter_thres, min_count = 10, num_core = 8) 
        saveRDS(SG_sp,file=fs::path(SVG_path,paste(data_name,"rare_CT_removed",SVG_method,sep="_"),ext="RDS"))},
        error = function(e){
            if(gs_filter_thres < 0.1){
                print(paste("error at gene filter threshold ",gs_filter_thres,"! Increasing threshold by 0.01",sep=""))
                print(e)
                gs_filter_thres=gs_filter_thres+0.01
                SG_sp = Giotto::spark(gobject = GO_obj, percentage = gs_filter_thres, min_count = 10, num_core = 8) 
                saveRDS(SG_sp,file=fs::path(SVG_path,paste(data_name,"rare_CT_removed",SVG_method,sep="_"),ext="RDS"))
            }else{
                print(paste("error at gene filter threshold ",gs_filter_thres,"! Setting threshold at 0.1",sep=""))
                SG_sp = Giotto::spark(gobject = GO_obj, percentage = 0.1, min_count = 10, num_core = 8) 
                saveRDS(SG_sp,file=fs::path(SVG_path,paste(data_name,"rare_CT_removed",SVG_method,sep="_"),ext="RDS"))
            }
            })
    }else if(SVG_method == "spatialDE"){
        #if(pearson_residual){
        #    expr_values = GO_obj@custom_expr
        #}else{
        #    expr_values = GO_obj@norm_expr
        #}
        expr_values = GO_obj@norm_expr
        print(dim(expr_values))
        #expr_values = ifelse(pearson_residual,GO_obj@custom_expr,GO_obj@norm_expr)
        reticulate::use_python(required = T, python = python_path)
        reader_path = system.file("python", "SpatialDE_wrapper.py", package = 'Giotto')
        reticulate::source_python(file = reader_path)

        ## get spatial locations
        spatial_locs = as.data.frame(GO_obj@spatial_locs)
        rownames(spatial_locs) = spatial_locs$cell_ID
        spatial_locs = subset(spatial_locs, select = -cell_ID)

        print(dim(spatial_locs))
        SG_sde = Spatial_DE(as.data.frame(t(as.matrix(expr_values))), spatial_locs)
        SG_sp = SG_sde
        SG_sp = SG_sp$ms_results
        SG_sp = SG_sp %>% rename(adjusted_pvalue = qval, combined_pvalue = pval, genes = g)
        saveRDS(SG_sde,file=fs::path(SVG_path,paste(data_name,"rare_CT_removed",SVG_method,sep="_"),ext="RDS"))
    }
  }else{
    print("load SVG")
    SG_sp = readRDS(fs::path(SVG_path,paste(data_name,"rare_CT_removed",SVG_method,sep="_"),ext="RDS"))
    if(SVG_method=="spatialDE"){
        SG_sp = SG_sp$ms_results
        SG_sp = SG_sp %>% rename(adjusted_pvalue = qval, combined_pvalue = pval, genes = g)
    }else if(SVG_method=="SPARK" & "list" %in% class(SG_sp)){
      SG_sp = SG_sp$SVG_df
    }
  }

  print("set SVG levels and run PCA")
  SG_sp = SG_sp %>% arrange(adjusted_pvalue, combined_pvalue)
  SG_sp = SG_sp %>% filter(adjusted_pvalue < 0.05)

  SVG_sp_low = SG_sp %>% filter(adjusted_pvalue < 0.05) %>% pull(genes)
  SVG_sp_med = SG_sp %>% filter(adjusted_pvalue <= quantile(SG_sp %>% filter(adjusted_pvalue < 0.05) %>% pull(adjusted_pvalue),0.25)) %>% pull(genes)
  SVG_sp_high = SG_sp %>% filter(adjusted_pvalue <= quantile(SG_sp %>% filter(adjusted_pvalue < 0.05) %>% pull(adjusted_pvalue),0.5)) %>% pull(genes)

  GO_obj@gene_metadata[[paste(SVG_method,"svg_low",sep="_")]] = as.numeric(GO_obj@gene_ID %in% SVG_sp_low)
  GO_obj@gene_metadata[[paste(SVG_method,"svg_med",sep="_")]] = as.numeric(GO_obj@gene_ID %in% SVG_sp_med)
  GO_obj@gene_metadata[[paste(SVG_method,"svg_high",sep="_")]] = as.numeric(GO_obj@gene_ID %in% SVG_sp_high)
  
  ## PCA of SVGs
  pca_scale = pca_center = ifelse(pearson_residual,FALSE,TRUE)
  pca_values = ifelse(pearson_residual,"custom","normalized")
  svg_low_s = Sys.time()
  GO_obj = runPCA(GO_obj, expression_values = pca_values, name = paste("pca",SVG_method,"svg_low",sep="_"), genes_to_use = SVG_sp_low, center = pca_center, scale_unit = pca_scale, ncp = min(length(SVG_sp_low),100), method="irlba")
  svg_low_t = Sys.time()

  svg_med_s = Sys.time()
  GO_obj = runPCA(GO_obj, expression_values = pca_values, name = paste("pca",SVG_method,"svg_med",sep="_"), genes_to_use = SVG_sp_med, center = pca_center, scale_unit = pca_scale, ncp = min(length(SVG_sp_med),100), method="irlba")
  svg_med_t = Sys.time()
  
  svg_high_s = Sys.time()
  GO_obj = runPCA(GO_obj, expression_values = pca_values, name = paste("pca",SVG_method,"svg_high",sep="_"), genes_to_use = SVG_sp_high, center = pca_center, scale_unit = pca_scale, ncp = min(length(SVG_sp_high),100), method="irlba")
  svg_high_t = Sys.time()
  
  screePlot(GO_obj, show_plot=FALSE, name=paste("pca",SVG_method,"svg_low",sep="_"), method="irlba", genes_to_use = paste(SVG_method,"svg_low",sep="_"), scale_unit = pca_scale, center = pca_center, save_plot=TRUE, default_save_name=paste("screeplot",SVG_method,"svg_low",data_name,sep="_"),
            save_param=list(base_width = 18, base_height = 7, save_dir = SVG_path))
  screePlot(GO_obj, show_plot=FALSE, name=paste("pca",SVG_method,"svg_med",sep="_"), method="irlba", genes_to_use = paste(SVG_method,"svg_med",sep="_"), scale_unit = pca_scale, center = pca_center, save_plot=TRUE, default_save_name=paste("screeplot",SVG_method,"svg_med",data_name,sep="_"),
            save_param=list(base_width = 18, base_height = 7, save_dir = SVG_path))
  screePlot(GO_obj, show_plot=FALSE, name=paste("pca",SVG_method,"svg_high",sep="_"), method="irlba", genes_to_use = paste(SVG_method,"svg_high",sep="_"), scale_unit = pca_scale, center = pca_center, save_plot=TRUE, default_save_name=paste("screeplot",SVG_method,"svg_high",data_name,sep="_"),
            save_param=list(base_width = 18, base_height = 7, save_dir = SVG_path))

  ## run HVG
  HVG_method = ifelse(pearson_residual,"sct","loess")
  if(!dir.exists(HVG_path)){
    dir.create(HVG_path,recursive=TRUE,showWarnings=FALSE)
  }
  if(run_hvg){
    print("compute HVGs")
    if(pearson_residual){
      HVG_df = GO_obj@gene_metadata
      npool = HVG_df %>% filter(residual_variance > 1) %>% nrow()
      if(npool < geneset_cutoff){
        cutoffs = c(0.1,0.3,0.5)
      }else{
        cutoffs = c(0.5,0.7,0.9)
      }
      HVG_thres = quantile(HVG_df %>% filter(residual_variance > 1) %>% pull(residual_variance),cutoffs)
      HVG_thres = unname(HVG_thres)
      colnames(HVG_df)[which(colnames(HVG_df)=="residual_variance")] = "thres"
      
      saveRDS(HVG_df,fs::path(HVG_path,paste(data_name,"rare_CT_removed",HVG_method,sep="_"),ext="RDS")) 
    }else{
      HVG_df = calculateHVG(GO_obj,method="cov_loess",difference_in_cov=0.1,return_gobject=F,show_plot=F,return_plot=F,save_plot=F)
      npool = HVG_df %>% filter(cov_diff > 0) %>% nrow()
      if(npool < geneset_cutoff){
        cutoffs = c(0.1,0.3,0.5)
      }else{
        cutoffs = c(0.5,0.7,0.9)
      }
      HVG_thres = quantile(HVG_df %>% filter(cov_diff > 0) %>% pull(cov_diff),cutoffs)
      HVG_thres = unname(HVG_thres)
      colnames(HVG_df)[which(colnames(HVG_df)=="cov_diff")] = "thres"
      colnames(HVG_df)[which(colnames(HVG_df)=="genes")] = "gene_ID"
      
      saveRDS(HVG_df,fs::path(HVG_path,paste(data_name,"rare_CT_removed",HVG_method,sep="_"),ext="RDS"))
    }
  }else{
    HVG_df = readRDS(fs::path(HVG_path,paste(data_name,"rare_CT_removed",HVG_method,sep="_"),ext="RDS"))
    if(pearson_residual){
      colnames(HVG_df)[which(colnames(HVG_df)=="residual_variance")] = "thres"
      npool = HVG_df %>% filter(thres > 1) %>% nrow()
      if(npool < geneset_cutoff){
        cutoffs = c(0.1,0.3,0.5)
      }else{
        cutoffs = c(0.5,0.7,0.9)
      }
      HVG_thres = quantile(HVG_df %>% filter(thres > 1) %>% pull(thres),cutoffs)
      HVG_thres = unname(HVG_thres)
    }else{
      colnames(HVG_df)[which(colnames(HVG_df)=="cov_diff")] = "thres"
      colnames(HVG_df)[which(colnames(HVG_df)=="genes")] = "gene_ID"
      npool = HVG_df %>% filter(thres > 0) %>% nrow()
      if(npool < geneset_cutoff){
        cutoffs = c(0.1,0.3,0.5)
      }else{
        cutoffs = c(0.5,0.7,0.9)
      }
      HVG_thres = quantile(HVG_df %>% filter(thres > 0) %>% pull(thres),cutoffs)
      HVG_thres = unname(HVG_thres)
    }
  }
  
  print("set HVG levels and run PCA")
  ## low threshold
  HVG_low = setdiff(HVG_df %>% filter(thres >= HVG_thres[1]) %>% pull(gene_ID),gs_filtered)
  #GO_obj@gene_metadata = GO_obj@gene_metadata %>% mutate(hvg_low = ifelse(gene_ID %in% HVG_low, "yes", "no"))
  GO_obj@gene_metadata[[paste(HVG_method,"hvg_low",sep="_")]] = ifelse(GO_obj@gene_metadata$gene_ID %in% HVG_low, "yes", "no")
  dir.create(fs::path(HVG_path,"low"),showWarnings = FALSE, recursive=TRUE)
  saveRDS(HVG_low, fs::path(HVG_path,"low",paste(data_name,"rare_CT_removed",HVG_method,sep="_"),ext="RDS"))
  
  ## medium threshold
  HVG_med = setdiff(HVG_df %>% filter(thres >= HVG_thres[2]) %>% pull(gene_ID),gs_filtered)
  #GO_obj@gene_metadata = GO_obj@gene_metadata %>% mutate(hvg_med = ifelse(gene_ID %in% HVG_med, "yes", "no"))
  GO_obj@gene_metadata[[paste(HVG_method,"hvg_med",sep="_")]] = ifelse(GO_obj@gene_metadata$gene_ID %in% HVG_med, "yes", "no")
  dir.create(fs::path(HVG_path,"med"),showWarnings = FALSE, recursive=TRUE)
  saveRDS(HVG_med, fs::path(HVG_path,"med",paste(data_name,"rare_CT_removed",HVG_method,sep="_"),ext="RDS"))
  
  ## high threshold
  HVG_high = setdiff(HVG_df %>% filter(thres >= HVG_thres[3]) %>% pull(gene_ID),gs_filtered)
  #GO_obj@gene_metadata = GO_obj@gene_metadata %>% mutate(hvg_high = ifelse(gene_ID %in% HVG_high, "yes", "no"))
  GO_obj@gene_metadata[[paste(HVG_method,"hvg_high",sep="_")]] = ifelse(GO_obj@gene_metadata$gene_ID %in% HVG_high, "yes", "no")
  dir.create(fs::path(HVG_path,"high"),showWarnings = FALSE, recursive=TRUE)
  saveRDS(HVG_high, fs::path(HVG_path,"high",paste(data_name,"rare_CT_removed",HVG_method,sep="_"),ext="RDS"))
  
  ## PCA of HVGs
  hvg_low_s = Sys.time()
  GO_obj = runPCA(GO_obj, expression_values = pca_values, scale_unit = pca_scale, center = pca_center, genes_to_use = paste(HVG_method,"hvg_low",sep="_"), name = paste("pca",HVG_method,"hvg_low",sep="_"), ncp = min(100, length(HVG_low)), method="irlba")
  hvg_low_t = Sys.time()
  
  hvg_med_s = Sys.time()
  GO_obj = runPCA(GO_obj, expression_values = pca_values, scale_unit = pca_scale, center = pca_center, genes_to_use = paste(HVG_method,"hvg_med",sep="_"), name = paste("pca",HVG_method,"hvg_med",sep="_"), ncp = min(100, length(HVG_med)), method="irlba")
  hvg_med_t = Sys.time()

  hvg_high_s = Sys.time()
  GO_obj = runPCA(GO_obj, expression_values = pca_values, scale_unit = pca_scale, center = pca_center, genes_to_use = paste(HVG_method,"hvg_high",sep="_"), name = paste("pca",HVG_method,"hvg_high",sep="_"), ncp = min(100, length(HVG_high)), method="irlba")
  hvg_high_t = Sys.time()

  screePlot(GO_obj, show_plot=FALSE, name=paste("pca",HVG_method,"hvg_low",sep="_"), method="irlba", genes_to_use = paste(HVG_method,"hvg_low",sep="_"), scale_unit = pca_scale, center = pca_center, save_plot=TRUE, default_save_name=paste("screeplot",HVG_method,"hvg_low",data_name,sep="_"),
            save_param=list(base_width = 18, base_height = 7, save_dir = HVG_path))
  screePlot(GO_obj, show_plot=FALSE, name=paste("pca",HVG_method,"hvg_med",sep="_"), method="irlba", genes_to_use = paste(HVG_method,"hvg_med",sep="_"), scale_unit = pca_scale, center = pca_center, save_plot=TRUE, default_save_name=paste("screeplot",HVG_method,"hvg_med",data_name,sep="_"),
            save_param=list(base_width = 18, base_height = 7, save_dir = HVG_path))
  screePlot(GO_obj, show_plot=FALSE, name=paste("pca",HVG_method,"hvg_high",sep="_"), method="irlba", genes_to_use = paste(HVG_method,"hvg_high",sep="_"), scale_unit = pca_scale, center = pca_center, save_plot=TRUE, default_save_name=paste("screeplot",HVG_method,"hvg_high",data_name,sep="_"),
            save_param=list(base_width = 18, base_height = 7, save_dir = HVG_path))
  
  ## union gene sets
  print("set union gene sets and run PCA")
  dir.create(union_path,recursive=TRUE,showWarnings = FALSE)
  ## HVG_low, SVG_low
  union_LL = union(HVG_low,SVG_sp_low)
  union_LL_s = Sys.time()
  GO_obj = runPCA(GO_obj, expression_values = pca_values, name = paste("pca_union_LL",HVG_method,SVG_method,sep="_"), genes_to_use = union_LL, center= pca_center, scale = pca_scale, ncp = min(length(union_LL),100), method="irlba")
  union_LL_t = Sys.time()
  screePlot(GO_obj, show_plot=FALSE, name=paste("pca_union_LL",HVG_method,SVG_method,sep="_"), method="irlba", genes_to_use = union_LL, scale = pca_scale, center= pca_center, save_plot=TRUE, default_save_name=paste("screeplot_union_LL",HVG_method,SVG_method,data_name,sep="_"),
            save_param=list(base_width = 18, base_height = 7, save_dir = union_path))
  ## HVG_low, SVG_med
  union_LM = union(HVG_low,SVG_sp_med)
  union_LM_s = Sys.time()
  GO_obj = runPCA(GO_obj, expression_values = pca_values, name = paste("pca_union_LM",HVG_method,SVG_method,sep="_"), genes_to_use = union_LM, center= pca_center, scale = pca_scale, ncp = min(length(union_LM),100), method="irlba")
  union_LM_t = Sys.time()
  screePlot(GO_obj, show_plot=FALSE, name=paste("pca_union_LM",HVG_method,SVG_method,sep="_"), method="irlba", genes_to_use = union_LM, scale = pca_scale, center= pca_center, save_plot=TRUE, default_save_name=paste("screeplot_union_LM",HVG_method,SVG_method,data_name,sep="_"),
            save_param=list(base_width = 18, base_height = 7, save_dir = union_path))
  # HVG_low, SVG_high
  union_LH = union(HVG_low,SVG_sp_high)
  union_LH_s = Sys.time()
  GO_obj = runPCA(GO_obj, expression_values = pca_values, name = paste("pca_union_LH",HVG_method,SVG_method,sep="_"), genes_to_use = union_LH, center= pca_center, scale = pca_scale, ncp = min(length(union_LH),100), method="irlba")
  union_LH_t = Sys.time()
  screePlot(GO_obj, show_plot=FALSE, name=paste("pca_union_LH",HVG_method,SVG_method,sep="_"), method="irlba", genes_to_use = union_LH, scale = pca_scale, center= pca_center, save_plot=TRUE, default_save_name=paste("screeplot_union_LH",HVG_method,SVG_method,data_name,sep="_"),
            save_param=list(base_width = 18, base_height = 7, save_dir = union_path))
  # HVG_med, SVG_low
  union_ML = union(HVG_med,SVG_sp_low)
  union_ML_s = Sys.time()
  GO_obj = runPCA(GO_obj, expression_values = pca_values, name = paste("pca_union_ML",HVG_method,SVG_method,sep="_"), genes_to_use = union_ML, center= pca_center, scale = pca_scale, ncp = min(length(union_ML),100), method="irlba")
  union_ML_t = Sys.time()
  screePlot(GO_obj, show_plot=FALSE, name=paste("pca_union_ML",HVG_method,SVG_method,sep="_"), method="irlba", genes_to_use = union_ML, scale = pca_scale, center= pca_center, save_plot=TRUE, default_save_name=paste("screeplot_union_ML",HVG_method,SVG_method,data_name,sep="_"),
            save_param=list(base_width = 18, base_height = 7, save_dir = union_path))
  # HVG_med, SVG_med 
  union_MM = union(HVG_med,SVG_sp_med)
  union_MM_s = Sys.time()
  GO_obj = runPCA(GO_obj, expression_values = pca_values, name = paste("pca_union_MM",HVG_method,SVG_method,sep="_"), genes_to_use = union_MM, center= pca_center, scale = pca_scale, ncp = min(length(union_MM),100), method="irlba")
  union_MM_t = Sys.time()
  screePlot(GO_obj, show_plot=FALSE, name=paste("pca_union_MM",HVG_method,SVG_method,sep="_"), method="irlba", genes_to_use = union_MM, scale = pca_scale, center= pca_center, save_plot=TRUE, default_save_name=paste("screeplot_union_MM",HVG_method,SVG_method,data_name,sep="_"),
            save_param=list(base_width = 18, base_height = 7, save_dir = union_path))
  # HVG_med, SVG_high
  union_MH = union(HVG_med,SVG_sp_high)
  union_MH_s = Sys.time()
  GO_obj = runPCA(GO_obj, expression_values = pca_values, name = paste("pca_union_MH",HVG_method,SVG_method,sep="_"), genes_to_use = union_MH, center= pca_center, scale = pca_scale, ncp = min(length(union_MH),100), method="irlba")
  union_MH_t = Sys.time()
  screePlot(GO_obj, show_plot=FALSE, name=paste("pca_union_MH",HVG_method,SVG_method,sep="_"), method="irlba", genes_to_use = union_MH, scale = pca_scale, center= pca_center, save_plot=TRUE, default_save_name=paste("screeplot_union_MH",HVG_method,SVG_method,data_name,sep="_"),
            save_param=list(base_width = 18, base_height = 7, save_dir = union_path))
  # HVG_high, SVG_low
  union_HL = union(HVG_high,SVG_sp_low)
  union_HL_s = Sys.time()
  GO_obj = runPCA(GO_obj, expression_values = pca_values, name = paste("pca_union_HL",HVG_method,SVG_method,sep="_"), genes_to_use = union_HL, center= pca_center, scale = pca_scale, ncp = min(length(union_HL),100), method="irlba")
  union_HL_t = Sys.time()
  screePlot(GO_obj, show_plot=FALSE, name=paste("pca_union_HL",HVG_method,SVG_method,sep="_"), method="irlba", genes_to_use = union_HL, scale = pca_scale, center= pca_center, save_plot=TRUE, default_save_name=paste("screeplot_union_HL",HVG_method,SVG_method,data_name,sep="_"),
            save_param=list(base_width = 18, base_height = 7, save_dir = union_path))
  # HVG_high, SVG_med
  union_HM = union(HVG_high,SVG_sp_med)
  union_HM_s = Sys.time()
  GO_obj = runPCA(GO_obj, expression_values = pca_values, name = paste("pca_union_HM",HVG_method,SVG_method,sep="_"), genes_to_use = union_HM, center= pca_center, scale = pca_scale, ncp = min(length(union_HM),100), method="irlba")
  union_HM_t = Sys.time()
  screePlot(GO_obj, show_plot=FALSE, name=paste("pca_union_HM",HVG_method,SVG_method,sep="_"), method="irlba", genes_to_use = union_HM, scale = pca_scale, center= pca_center, save_plot=TRUE, default_save_name=paste("screeplot_union_HM",HVG_method,SVG_method,data_name,sep="_"),
            save_param=list(base_width = 18, base_height = 7, save_dir = union_path))
  # HVG_high, SVG_high
  union_HH = union(HVG_high,SVG_sp_high)
  union_HH_s = Sys.time()
  GO_obj = runPCA(GO_obj, expression_values = pca_values, name = paste("pca_union_HH",HVG_method,SVG_method,sep="_"), genes_to_use = union_HH, center= pca_center, scale = pca_scale, ncp = min(length(union_HH),100), method="irlba")
  union_HH_t = Sys.time()
  screePlot(GO_obj, show_plot=FALSE, name=paste("pca_union_HH",HVG_method,SVG_method,sep="_"), method="irlba", genes_to_use = union_HH, scale = pca_scale, center= pca_center, save_plot=TRUE, default_save_name=paste("screeplot_union_HH",HVG_method,SVG_method,data_name,sep="_"),
            save_param=list(base_width = 18, base_height = 7, save_dir = union_path))
  
  #### all
  print("run PCA for all genes")
  all_genes = GO_obj@gene_ID
  all_genes_s = Sys.time()
  GO_obj = runPCA(GO_obj, expression_values = pca_values, name = "pca_all", genes_to_use = all_genes, center= pca_center, scale = pca_scale, ncp = min(length(all_genes),100), method="irlba")
  all_genes_t = Sys.time()

  ### non_HVG_SVG_low
  if(!is.null(random_gs_reps)){
    non_union_LL = setdiff(all_genes,union_LL)
    if(length(non_union_LL) >= length(union_LL)){
        n_randoms_gs = length(union_LL)
        for(i in 1:random_gs_reps){
            set.seed(random_gs_seed+i)
            random_gs = sample(non_union_LL,size=n_randoms_gs,replace=FALSE)
            GO_obj@gene_metadata[[paste("non_union_LL",HVG_method,SVG_method,"rep",i,sep="_")]] = ifelse(GO_obj@gene_metadata$gene_ID %in% random_gs, "yes", "no")
            GO_obj = runPCA(GO_obj, expression_values = pca_values, name = paste("pca_non_union_LL",HVG_method,SVG_method,"rep",i,sep="_"), genes_to_use = non_union_LL, center= pca_center, scale = pca_scale, ncp = min(length(non_union_LL),100), method="irlba")
            screePlot(GO_obj, show_plot=FALSE, name=paste("non_union_LL",HVG_method,SVG_method,"rep",i,sep="_"), method="irlba", genes_to_use = non_union_LL, scale = pca_scale, center= pca_center, save_plot=TRUE, default_save_name=paste("screeplot_non_union_LL",HVG_method,SVG_method,"rep",i,data_name,sep="_"),
                        save_param=list(base_width = 18, base_height = 7, save_dir = union_path))
        }
    }
  }

  ##### compute entropy
  print("compute cell type entropy")
  GO_obj = baseline_entropy(GO_obj)
  
  ### save
  print("save preprocessed Giotto object for downstream")
  str(GO_obj)
  dir.create(save_path,recursive=T,showWarnings = F)
  end_times = c(hvg_low_t,hvg_med_t,hvg_high_t,svg_low_t,svg_med_t,svg_high_t,union_LL_t,union_LM_t,union_LH_t,union_ML_t,union_MM_t,union_MH_t,union_HL_t,union_HM_t,union_HH_t,all_genes_t)
  start_times = c(hvg_low_s,hvg_med_s,hvg_high_s,svg_low_s,svg_med_s,svg_high_s,union_LL_s,union_LM_s,union_LH_s,union_ML_s,union_MM_s,union_MH_s,union_HL_s,union_HM_s,union_HH_s,all_genes_s)
  runtime = as.numeric(difftime(end_times, start_times, units = "secs"))
  names(runtime) = c("hvg_low","hvg_med","hvg_high","svg_low","svg_med","svg_high","union_LL","union_LM","union_LH","union_ML","union_MM","union_MH","union_HL","union_HM","union_HH","all")
  GO_obj@parameters$runtime = runtime
  saveRDS(GO_obj,fs::path(save_path,paste(data_name,"rare_CT_removed_pre",sep="_"),ext="RDS"))
}

#########################
## debugging 
#########################
## library(dplyr)
## library(tidyr)
## library(Matrix)
## library(Giotto,lib.loc = "/sw/pkgs/med/Rgiotto/1.1.0")
## library(sctransform)
## 
## source("/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/functions/helper_funcs.R")
## 
## dat = readRDS("/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/GB_rev_4/real_data_analysis/data/single_cell/Xenium/raw/mouse_brain_ATN_counts.1.1.1.1.1.1.RDS")
## 
## raw_exprs = dat$raw_exprs 
## spatial_locs = dat$spatial_locs
## cell_metadata = dat$cell_metadata
## scale_factor = 6000
## cell_type_col = "SCT_snn_res.0.4"
## python_path = "/sw/pkgs/arc/python3.9-anaconda/2021.11/bin"
## remove_rare_CT = TRUE 
## rare_CT_threshold = 0.05
## cell_filter = 0.03 
## gene_filter = "by_pca"
## pearson_residual =  TRUE 
## save_path = "/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/GB_rev_4/real_data_analysis_trouble_shoot/SCT_adjusted/feature_selection"
## run_svg = TRUE
## SVG_path = "/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/GB_rev_4/real_data_analysis_trouble_shoot/SCT_adjusted/results/SVG"
## run_hvg = TRUE 
## HVG_path = "/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/GB_rev_4/real_data_analysis_trouble_shoot/SCT_adjusted/results/HVG" 
## union_path = "/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/GB_rev_4/real_data_analysis_trouble_shoot/SCT_adjusted/results/union"
## data_name = "mouse_brain_ATN_counts.1.1.1.1.1.1" 
## QC_path = "/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/GB_rev_4/real_data_analysis_trouble_shoot/SCT_adjusted/QC"