spat_metrics = function(spat_df,clus_res,clus_col_name,width){
  #df = inner_join(spat_df,clus_res$GO@cell_metadata %>% select(c("cell_ID",clus_col_name)))
  df = inner_join(spat_df,clus_res %>% dplyr::select(c("cell_ID",clus_col_name)))
  local_ari = c()
    local_ari_names=c()
  for(i in 1:(dim(df)[1])){
    pt_x = df$sdimx[i]
    pt_y = df$sdimy[i]
    label_df_tmp = df %>% 
      filter(sdimx >= pt_x-width & sdimx <= pt_x+width) %>%
      filter(sdimy >= pt_y-width & sdimy <= pt_y+width)
      #filter(dplyr::between(sdimx, left = pt_x-width, right = pt_x+width)) %>%
      #filter(dplyr::between(sdimy, pt_y-width, pt_y+width))
    if(length(unique(label_df_tmp$cell_type))==1 & length(unique(label_df_tmp[[clus_col_name]]))==1){
      tmp = 1
    }else{
      tmp = AMI(label_df_tmp$cell_type,label_df_tmp[[clus_col_name]])
    }
    local_ari = c(local_ari, tmp)
       local_ari_names = c(local_ari_names, df$cell_ID[i])
  }
    names(local_ari) = local_ari_names
  df$entro = df$entro/sum(df$entro)
  
  AIC_df = data.frame(sdimx = df$sdimx, sdimy = df$sdimy, truth = df$cell_type, cl = df[[clus_col_name]])
  return(list(spat_ami = local_ari, SC = sum(df$entro * local_ari) ,VN = spat_metric(AIC_df, method = "Var_norm")))
}
