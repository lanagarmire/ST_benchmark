Leiden_clus_res_grid_search_init = function(n_CT, match_thres,round_thres,res_start,res_end,res_by,GO,seed_vec,n_iter,geneset,nn_network_name,CT_bin_size){
  res_df = seed_df = c()
  n_round = 0
  n_clus_vec = c(0)
  n_matches=-1
  grid_counter=0
  
  no_match=TRUE
  while(no_match & (n_round < round_thres)){
    print(n_round)
    res_df_tmp = seed_df_tmp = c()
    
    ## make resolution values grid
    res_vec = setdiff(seq(res_start,res_end,by=res_by),c(0))
    print(paste0("resolution: ",res_vec))
    
    for(i in 1:length(res_vec)){
      res = res_vec[i]
      for(j in 1:length(seed_vec)){
        seed = seed_vec[j]
        GO = doLeidenCluster(gobject = GO, resolution = res, n_iterations = n_iter, seed_number = seed, name = paste("leiden_clus",geneset,res,"seed",seed,sep="_"), network_name = nn_network_name)
        
        res_df_tmp = c(res_df_tmp,res)
        seed_df_tmp = c(seed_df_tmp,seed)
      }
    }
    n_clus_vec = sapply(paste("leiden_clus",geneset,res_df_tmp,"seed",seed_df_tmp,sep="_"),FUN=function(i){
      mets = GO@cell_metadata%>%pull(i)%>%table()%>%length()})
    print(paste0("range of cluster numbers: ",range(n_clus_vec)))
    
    n_clus_df = data.frame(res = res_df_tmp, seed = seed_df_tmp, n_clus = n_clus_vec)
    n_matches = n_clus_df %>% mutate(match=n_clus==n_CT) %>% group_by(res) %>% summarise(n=sum(match)) %>% pull(n) %>% max()
    print(paste0("number of matches: ",n_matches))
       if(n_matches > 0){
            no_match=FALSE
          }
    
    n_CT_lower = max(n_clus_df$n_clus[n_clus_df$n_clus <= n_CT-CT_bin_size])
    print(paste0("n_CT_lower: ",n_CT_lower))
    
    n_CT_upper = min(n_clus_df$n_clus[n_clus_df$n_clus >= n_CT+CT_bin_size])
    print(paste0("n_CT_upper: ",n_CT_upper))
    
    res_start_cur = res_start
    res_end_cur = res_end
    if(is.infinite(n_CT_lower)){
      res_start = res_start/2
      #res_start_nCT = -1
    }else{
      res_start = min(n_clus_df %>% dplyr::filter(n_clus == n_CT_lower) %>% pull(res))
      #res_start_nCT = n_clus_df %>% dplyr::filter(res==res_start) %>% pull(n_clus)
    }
    print(paste0("new res_start= ",res_start))
    
    if(is.infinite(n_CT_upper)){
      ## if no number of clusters above n_CT+bin_size, increase resolution
      res_end=min(res_end*2,1)
    }else{
      ## obtain LARGEST resolution corresponding to n_CT_upper
      res_end = max(n_clus_df %>% dplyr::filter(n_clus == n_CT_upper) %>% pull(res))
    }
    print(paste0("new res_end= ",res_end))
    
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
  return(list(GO=GO,res=res_df,seed=seed_df,n_rounds=n_round,res_start = res_start, res_end = res_end))
}



Leiden_clus_res_grid_search = function(n_CT, match_thres,round_thres,res_start,res_end,res_by,GO,seed_vec,n_iter,geneset,nn_network_name){
  res_df = seed_df = c()
  n_round = 0
  n_clus_vec = c(0)
  n_matches=-1
  orig_meta = GO@cell_metadata
  grid_counter = 0
  
  while((n_matches < match_thres) & (n_round < round_thres)){
    print(n_round)
    res_df_tmp = seed_df_tmp = c()
    res_vec = setdiff(sort(unique(c(res_start,res_end,seq(res_start,res_end,by=res_by)))),res_df)
    
    print(head(GO@cell_metadata))
    print(paste0("resolution: ",res_vec))
    for(i in 1:length(res_vec)){
      res = res_vec[i]
      for(j in 1:length(seed_vec)){
        seed = seed_vec[j]
        clus_name = paste("leiden_clus",geneset,res,"seed",seed,sep="_")
        print(clus_name)
        if(!(clus_name %in% colnames(GO@cell_metadata))){
          print("clustering")
          GO = doLeidenCluster(gobject = GO, resolution = res, n_iterations = n_iter, seed_number = seed, name = clus_name, network_name = nn_network_name)
        }
        res_df_tmp = c(res_df_tmp,res)
        seed_df_tmp = c(seed_df_tmp,seed)
      }
    }
    print("compute n_clus")
    n_clus_vec = sapply(paste("leiden_clus",geneset,res_df_tmp,"seed",seed_df_tmp,sep="_"),FUN=function(i){
      mets = GO@cell_metadata%>%pull(i)%>%table()%>%length()})
    print(paste0("range of cluster numbers: ",range(n_clus_vec)))
    n_clus_df = data.frame(res = res_df_tmp, seed = seed_df_tmp, n_clus = n_clus_vec)
    n_matches = n_clus_df %>% mutate(match=n_clus==n_CT) %>% group_by(res) %>% summarise(n=sum(match)) %>% pull(n) %>% max()
    print(paste0("number of matches: ",n_matches))
    
    ### if there weren't enough matches, then discard the clsutering results from this round
    if(n_matches < match_thres){
      print("removing non-match results")
      GO@cell_metadata = orig_meta
    }else{
      res_df = c(res_df, res_df_tmp)
      seed_df = c(seed_df, seed_df_tmp)
    }
    
    n_CT_lower = max(n_clus_df$n_clus[n_clus_df$n_clus < n_CT])
    n_CT_upper = min(n_clus_df$n_clus[n_clus_df$n_clus > n_CT])
    ## record current values of res_start and res_end
    res_start_cur = res_start
    res_end_cur = res_end
    if(is.infinite(n_CT_lower)){
      ## if no number of clusters below n_CT, lower resolution
      res_start = res_start / 2
    }else{
      ## obtain new res_start value
      ## gather all resolution values corresponding to n_CT_lower
      res_start_values = n_clus_df %>% dplyr::filter(n_clus == n_CT_lower) %>% pull(res)
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
      res_end_values = n_clus_df %>% dplyr::filter(n_clus == n_CT_upper) %>% pull(res)
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
    
    n_round = n_round + 1
  }
  return(list(GO=GO,res=res_df,seed=seed_df,n_rounds=n_round))
}


## if this still doesn't work then only keep results with matches over threshold
