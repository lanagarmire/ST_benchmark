getCluster_monocle3_init = function(gene_expr,reduce = FALSE, metadata,red_method, pcs,red_dim_df, data_type, match_thres,round_thres, res_start,res_end,res_by,clus_method,seed_vec, CT_bin_size, max_clusters,
                               num_iter=1,num_neighbors=10, ncl_truth){
  library(monocle3)
  
  non_zero_variance_cells = names(which(apply(gene_expr, 2, sd) != 0))
  gene_expr = gene_expr[,colnames(gene_expr) %in% non_zero_variance_cells,drop=FALSE]
  
  cds = new_cell_data_set(gene_expr)
  cds = cds[,Matrix::colSums(exprs(cds)) != 0]
  if(reduce){
    cds = preprocess_cds(cds, num_dim = pcs, method="PCA", norm_method="none", scaling = su, verbose=pb)
  }else{
    colnames(red_dim_df) = paste0("PC",1:ncol(red_dim_df))
    reducedDims(cds)$PCA = red_dim_df[,1:pcs]
  }

  res_df = seed_df = c()
  n_round = 0
  n_clus_vec = c(0)
  n_matches=-1
  grid_counter=0

  no_match=TRUE
  while(no_match& (n_round < round_thres)){
    print(n_round)
    ## track resolution, seed, and number of clusters in the CURRENT round
    res_df_tmp = seed_df_tmp = n_clus_vec = c()
    print(c(res_start,res_end,res_by))
    
    ## make resolution values grid
    res_vec = setdiff(seq(res_start,res_end,by=res_by),c(0))
    
    ### run clustering on the given resolution grids
    for(i in 1:length(res_vec)){
      res = res_vec[i]
      for(j in 1:length(seed_vec)){ ## during initialization, we only run one rep for each resolution value
        seed = seed_vec[j]
        cds = cluster_cells(cds, reduction_method = red_method, k=num_neighbors, resolution=res, cluster_method=clus_method, num_iter=num_iter, random_seed=seed)
        #print(length(table(cds@clusters[[red_method]]$clusters)))
        
        ## Since monocle3 has a tendency of creating too many clusters, check if obtained number of clusters is reasonably small
        ## if so, save the number of clusters, resolution and seed values
        if(length(table(cds@clusters[[red_method]]$clusters)) < max_clusters){
          ## set name of clustering result
          col_name = paste("M3_clus",data_type,res,"seed",seed,sep="_")
          #print(length(table(cds@clusters[[red_method]]$clusters)))
          
          ## save clustering parameters
          res_df_tmp = c(res_df_tmp,res)
          seed_df_tmp = c(seed_df_tmp,seed)
          n_clus_vec = c(n_clus_vec,length(table(cds@clusters[[red_method]]$clusters)))
        }
      }
    }
    print(paste0("range of cluster numbers: ",range(n_clus_vec)))
    n_clus_df = data.frame(res = res_df_tmp, seed = seed_df_tmp, n_clus = n_clus_vec)
    
    ## obtain the resolution range that corresponds to n_CT +- bin_size
    ## closet number of clusters under n_CT - bin_size
    n_CT_lower = max(n_clus_df$n_clus[n_clus_df$n_clus <= n_CT-CT_bin_size])
    print(paste0("n_CT_lower: ",n_CT_lower))
    ## closet number of clusters above n_CT + bin_size
    n_CT_upper = min(n_clus_df$n_clus[n_clus_df$n_clus >= n_CT+CT_bin_size])
    print(paste0("n_CT_upper: ",n_CT_upper))
    
    ## update res_start, res_end, res_by
    res_start_cur = res_start
    res_end_cur = res_end
    if(is.infinite(n_CT_lower)){
      ## if no number of clusters below n_CT-bin_size, lower resolution
      res_start = res_start/2
    }else{
      ## obtain SMALLEST resolution corresponding to n_CT_lower
      res_start = min(n_clus_df %>% filter(n_clus == n_CT_lower) %>% pull(res))
    }
    print(paste0("new res_start= ",res_start))
    
    if(is.infinite(n_CT_upper)){
      ## if no number of clusters above n_CT+bin_size, increase resolution
      res_end=min(res_end*2,1)
    }else{
      ## obtain LARGEST resolution corresponding to n_CT_upper
      res_end = max(n_clus_df %>% filter(n_clus == n_CT_upper) %>% pull(res))
    }
    print(paste0("new res_end= ",res_end))
    
    ### check stopping criterion: obtain number of clusters range of n_CT+_bin_size
    if ((n_CT_lower == (n_CT-CT_bin_size)) & (n_CT_upper == (n_CT+CT_bin_size))) {
      no_match = FALSE
    } else if(res_start_cur == res_start & res_end_cur == res_end){
      ## if stoping criterion not satified AND res_start or res_end did not change, use a finer grid
      grid_counter = grid_counter + 1
      res_by = (res_end - res_start) / 10 / grid_counter
    }else{
      ## update res_by according to new res_start and res_end
      res_by = (res_end-res_start)/10
    }

    ## update round counter, record resolution values and seed values used at this round
    res_df = c(res_df, res_df_tmp)
    seed_df = c(seed_df, seed_df_tmp)
    n_round = n_round + 1
  }
  
  print(str(res_df))
  print(str(seed_df))
  return(list(res = res_df, seed = seed_df, n_rounds = n_round, res_start = res_start, res_end = res_end, M3_obj = cds))
}

getCluster_monocle3 = function(gene_expr,reduce = FALSE, metadata,red_method, pcs,red_dim_df, data_type, match_thres, round_thres,res_start,res_end,res_by,clus_method,seed_vec, max_clusters,
                               num_iter=1,num_neighbors=10, ncl_truth){
  library(monocle3)
  
  non_zero_variance_cells = names(which(apply(gene_expr, 2, sd) != 0))
  gene_expr = gene_expr[,colnames(gene_expr) %in% non_zero_variance_cells,drop=FALSE]
  
  cds = new_cell_data_set(gene_expr)
  cds = cds[,Matrix::colSums(exprs(cds)) != 0]
  if(reduce){
    cds = preprocess_cds(cds, num_dim = pcs, method="PCA", norm_method="none", scaling = su, verbose=pb)
  }else{
    colnames(red_dim_df) = paste0("PC",1:ncol(red_dim_df))
    reducedDims(cds)$PCA = red_dim_df[,1:pcs]
  }

  res_df = seed_df = c()
  n_round = 0
  n_clus_vec = c(0)
  n_matches=-1
  clus_res_df = metadata
  grid_counter = 0
  
  while(n_matches < match_thres & n_round < round_thres){
    print(n_round)
    ## track resolution, seed, and number of clusters in the CURRENT round
    res_df_tmp = seed_df_tmp = n_clus_vec = c()
    
    ## make resolution values grid
    print(c(res_start,res_end,res_by))
    res_vec = setdiff(seq(res_start,res_end,by=res_by),c(0))
    
    ########## run clustering for resolution grid, for each resolution value, run multiple (default: 10) times under different seeds
    for(i in 1:length(res_vec)){
      res = res_vec[i]
      for(j in 1:length(seed_vec)){
        seed = seed_vec[j]
        cds = cluster_cells(cds, reduction_method = red_method, k=num_neighbors, resolution=res, cluster_method=clus_method, num_iter=num_iter, random_seed=seed)
        
        ## Since monocle3 has a tendency of creating too many clusters, check if obtained number of clusters is reasonably small
        ## if so, save the number of clusters, resolution and seed values and clustering labels
        if(length(table(cds@clusters[[red_method]]$clusters)) < max_clusters){
          ## create data frame with clustering label assignment
          clus_res_df_tmp = data.frame(cell_ID = names(cds@clusters[[red_method]]$clusters), clus = as.numeric(cds@clusters[[red_method]]$clusters))
          col_name = paste("M3_clus",data_type,res,"seed",seed,sep="_")
          colnames(clus_res_df_tmp)[2] = col_name
          
          ## if clustering result is new, save it
          if(!(col_name %in% colnames(clus_res_df))){
            clus_res_df = merge(clus_res_df, clus_res_df_tmp, by="cell_ID", sort=FALSE)
            #print(length(table(cds@clusters[[red_method]]$clusters)))
            
            ## save clustering parameters
            res_df_tmp = c(res_df_tmp,res)
            seed_df_tmp = c(seed_df_tmp,seed)
            n_clus_vec = c(n_clus_vec,length(table(cds@clusters[[red_method]]$clusters)))
          }
        }
      }
    }
    print(paste0("range of cluster numbers: ",range(n_clus_vec)))
    n_clus_df = data.frame(res = res_df_tmp, seed = seed_df_tmp, n_clus = n_clus_vec)
    
    ######### gather the number of clusters corresponding each grid
    n_matches = n_clus_df %>% mutate(match=n_clus==n_CT) %>% group_by(res) %>% summarise(n=sum(match)) %>% pull(n) %>% max()
    print(paste0("number of matches: ",n_matches))
    
    ### if there weren't enough matches, then discard the clsutering results from this round
    if(n_matches < match_thres){
      print("removing non-match results")
      clus_res_df = metadata
    }else{
      res_df = c(res_df, res_df_tmp)
      seed_df = c(seed_df, seed_df_tmp)
    }
    
    ## obtain the resolution range that corresponds to n_CT +- bin_size
    ## closet number of clusters under n_CT
    n_CT_lower = max(n_clus_df$n_clus[n_clus_df$n_clus < n_CT])
    print(paste0("n_CT_lower: ",n_CT_lower))
    ## closet number of clusters above n_CT
    n_CT_upper = min(n_clus_df$n_clus[n_clus_df$n_clus > n_CT])
    print(paste0("n_CT_upper: ",n_CT_upper))
    
    ## record current values of res_start and res_end
    res_start_cur = res_start
    res_end_cur = res_end
    if(is.infinite(n_CT_lower)){
      ## if no number of clusters below n_CT, lower resolution
      res_start = res_start / 2
    }else{
      ## obtain new res_start value
      ## gather all resolution values corresponding to n_CT_lower
      res_start_values = n_clus_df %>% filter(n_clus == n_CT_lower) %>% pull(res)
      ## pick the resolution with highest frequency (since each were run multiple times)
      res_start_tab = table(res_start_values)
      print("res start values matching n_CT_lower")
      print(res_start_tab)
      res_start = as.numeric(names(res_start_tab)[which.max(res_start_tab)])
    }
    print(paste0("new res_start= ",res_start))
    
    if(is.infinite(n_CT_upper)){
      ## if no number of clusters above n_CT, increase resolution
      res_end=min(res_end*2,1)
    }else{
      ## obtain new res_end value
      ## gather all resolution values corresponding to n_CT_upper
      res_end_values = n_clus_df %>% filter(n_clus == n_CT_upper) %>% pull(res)
      ## pick the resolution with highest frequency (since each were run multiple times)
      res_end_tab = table(res_end_values)
      print("res end values matching n_CT_upper")
      print(res_end_tab)
      res_end = as.numeric(names(res_end_tab)[which.max(res_end_tab)])
    }
    print(paste0("new res_end= ",res_end))
    
    ## if the new res_start and res_end are equal, expand by in each direction by res_by
    if(res_start == res_end){
      res_start = res_start - res_by
      res_end = res_end + res_by
    }
    
    if(res_start_cur == res_start & res_end_cur == res_end){
      ## if res_start or res_end did not change, use a finer grid
      grid_counter = grid_counter + 1
      res_by = (res_end - res_start) / 10 / grid_counter
    }else{
      ## update res_by according to new res_start and res_end
      res_by = (res_end-res_start)/10
    }
    
    ## update round counter, record resolution values and seed values used at this round
    n_round = n_round + 1
  }
  
  print(str(res_df))
  print(str(seed_df))
  print(str(clus_res_df))
  return(list(res = res_df, seed = seed_df, n_rounds = n_round, clus_res_df = clus_res_df, M3_obj = cds))
}