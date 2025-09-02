#########################
### compute baseline entropy
#########################

baseline_entropy = function(gobject){
  spat_df = base::merge(as.data.frame(gobject@spatial_locs),as.data.frame(gobject@cell_metadata)%>%dplyr::select(c(cell_ID,cell_type)), by="cell_ID",sort=FALSE)
  dist_df = spat_df %>% dplyr::select(c(sdimx,sdimy)) %>% dist() %>% as.matrix()
  width = as.numeric(quantile(dist_df[upper.tri(dist_df,diag=FALSE)],0.01))
  entro = c()
  for(i in 1:(dim(spat_df)[1])){
    pt_x = spat_df$sdimx[i]
    pt_y = spat_df$sdimy[i]
    label_df_tmp = spat_df %>% dplyr::filter(dplyr::between(sdimx, left = pt_x-width, right = pt_x+width)) %>%
      dplyr::filter(dplyr::between(sdimy, pt_y-width, pt_y+width))
    entro = c(entro, Entropy(table(as.character(label_df_tmp$cell_type))))
  }
  zero_entro = which(entro==0)
  eps=1e-4
  entro[zero_entro] = eps/length(zero_entro)
  
  gobject@cell_metadata$entro = entro
  return(gobject)
}
