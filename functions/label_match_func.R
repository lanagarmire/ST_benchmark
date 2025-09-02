###################################
##### label matching 
###################################

label_match = function(clus_df,exprs_mat,method){
  print("cluster level expression")
  print(table(clus_df$clus))
  unique_clus = sort(unique(clus_df$clus))
  avg_expr_clus = sapply(unique_clus,FUN=function(x){
    cells = clus_df%>%pull(cell_ID)
    cells = cells[which(clus_df%>%pull(clus) == x)]
    cells_idx = which(colnames(exprs_mat) %in% cells)
    print(str(cells_idx))
    rowMeans(exprs_mat[,cells_idx,drop=FALSE])})
  
  print("cell type level expression")
  unique_grp = sort(unique(clus_df$group))
  avg_expr_grp = sapply(unique_grp,FUN=function(x){
    cells = clus_df%>%pull(cell_ID)
    cells = cells[which(clus_df%>%pull(group) == x)]
    cells_idx = which(colnames(exprs_mat) %in% cells)
    rowMeans(exprs_mat[,cells_idx,drop=FALSE])})
  
  print("compute cluster - cell type correlation")
  n_clus = length(unique_grp)
  grp_clus_tau = matrix(NA,n_clus,n_clus)
  for(row in 1:n_clus){
    for(col in 1:n_clus){
      grp_clus_tau[row,col] = cor.test(x=avg_expr_grp[,row],y=avg_expr_clus[,col],method="kendall")$estimate
    }
  }
  
  GS_match = galeShapley.marriageMarket(t(grp_clus_tau),grp_clus_tau) #proposal: group, engagement: clus
  print(n_clus)
  print(GS_match)
  
  args_ls = list(.x=as.numeric(as.factor(clus_df%>%pull(clus))))
  args_ls[2:(n_clus+1)] = NA
  names(args_ls)[2:(n_clus+1)] = as.character(c(GS_match$proposal))
  ## args_ls[2:(n_clus+1)] = 1:n_clus
  names(unique_grp)=as.character(1:n_clus)
  args_ls[2:(n_clus+1)] = unname(unique_grp[as.character(1:n_clus)])
  
  clus_df[,"clus"] = do.call(recode,args_ls)
  rename_col_idx = which(colnames(clus_df)=="clus")
  colnames(clus_df)[rename_col_idx] = method

  res = list(grp_clus_tau=grp_clus_tau,GS_match=GS_match,clus_df=clus_df)
  
  return(res)
}

#########################
## notes
#########################
#GS_match$proposal: row indices stand for the ordered labels in group, and the actual values stand for the matched labels in clus