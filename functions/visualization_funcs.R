
perm_test_median = function(group1,group2,paired) {
  if(paired){
    non_equal_idx = (group1!=group2)
    group1 = group1[non_equal_idx]
    group2 = group2[non_equal_idx]
    observed_stat = median(group1-group2,na.rm=TRUE)
  }else{
    observed_stat = median(group1)-median(group2)
  }
  test = wilcox.test(group1,group2,paired=paired,conf.int=TRUE,exact=TRUE,correct=FALSE)
  
  return(list(observed_stat = observed_stat, p_value = test$p.value, conf_int = test$conf.int))
}

make_boxplots_with_permutation_pval=function(data_df,metric,metric_name,panel_var,nGS,GS_test,GS_levels,
                                             cols,legend_symbol_size,label_rel_size,save_path,
                                             platform_width=12.16,platform_height=5.36){
  ### make pvalue boxplot
  library(ggplot2)
  dir.create(fs::path(save_path),recursive=TRUE,showWarnings = FALSE)
  p = data_df %>% 
    ggboxplot(x="GeneSet",y=metric,fill="GeneSet") + theme_bw() + scale_fill_manual(values=cols) + labs(fill="Gene Set") + ylab(metric_name) + xlab("") +
      theme(legend.text=element_text(size=rel(label_rel_size)), legend.title=element_text(size=rel(label_rel_size)),
            axis.text.x=element_text(size=rel(label_rel_size),angle=60,vjust=1,hjust=1), axis.text.y=element_text(size=rel(label_rel_size)), 
            axis.title.x=element_text(size=rel(label_rel_size)), axis.title.y=element_text(size=rel(label_rel_size))) + 
      guides(fill = guide_legend(override.aes = list(size = legend_symbol_size))) +
      facet_wrap(vars(!!sym(panel_var)),nrow=1) + 
      theme(strip.background = element_rect(fill="white"),strip.text.x = element_text(size = rel(label_rel_size)))
    
  ## compute pvalues
  xlocs = c(1:nGS)
  names(xlocs) = GS_levels
  print(xlocs)
  pvalue_df = NULL
  panel_vec = unique(data_df[[panel_var]])
  for(i in 1:length(panel_vec)){
    plat = panel_vec[i]
    dat = data_df %>% filter(!!sym(panel_var)==plat) %>% select(all_of(c("data",metric,"GeneSet"))) %>% pivot_wider(names_from=GeneSet,values_from=!!sym(metric))
    ymax = max(data_df %>% filter(!!sym(panel_var)==plat) %>% pull(!!sym(metric)),na.rm=TRUE)
    step = 0.12*diff(range(data_df %>% filter(!!sym(panel_var)==plat) %>% pull(!!sym(metric)),na.rm=TRUE))
    ypos = ymax
    
    for(j in 1:length(GS_test)){
      med_perm_test = perm_test_median(group1=dat[[GS_test[[j]][1]]],group2=dat[[GS_test[[j]][2]]],paired=T)
      print(med_perm_test$p_value)
      if(med_perm_test$p_value <= 0.1){
        pvalue_df = bind_rows(pvalue_df,data.frame(.y.=metric,group1=GS_test[[j]][1],group2=GS_test[[j]][2],p=med_perm_test$p_value,method="paired wilcox",type=plat,xmin=xlocs[GS_test[[j]][1]],xmax=xlocs[GS_test[[j]][2]],y.position=ypos+j*step))
      }
    }
  }

  if(!is.null(pvalue_df)){
    coltype_idx = which(colnames(pvalue_df)=="type")
    colnames(pvalue_df)[coltype_idx] = panel_var
    pvalue_df = pvalue_df %>% add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),symbols = c("***", "**", "*", ".",""))
    
    ## add in p-values
    p = p + stat_pvalue_manual(pvalue_df, inherit.aes = FALSE, hide.ns = TRUE,label = "{p.signif}")
  }
  
  ## save boxplot
  ggsave(plot=p, filename=fs::path(save_path,paste(metric,"med_platform_boxplot_pval",sep="_"),ext="png"),width=platform_width,height=platform_height,dpi=300,unit="in")
  return(pvalue_df)
}


perm_ci_median <- function(group1, group2, R=10000, alpha=0.05, n_grids=300) {
  non_equal_idx = (group1!=group2)
  group1 = group1[non_equal_idx]
  group2 = group2[non_equal_idx]
  
  test = wilcox.test(group1,group2,paired=TRUE,conf.int=TRUE,exact=TRUE,conf.level=0.9)
  return(c(test$estimate,test$conf.int))  
}

make_forestplot_with_effectsize_permutation_CI = function(data_df,clus_method,cols,metric,platform_vec,nGS,GS_test, GS_levels,n_perm,save_path,platform_width=12.16,platform_height=5.36){
  library(ggplot2)
  colnames(data_df)[which(colnames(data_df)==metric)]="met"
  dir.create(fs::path(save_path,clus_method),recursive=TRUE,showWarnings = FALSE)
  
  ## compute CIs
  forest_df = NULL
  for(i in 1:length(platform_vec)){
    plat = platform_vec[i]
    print(plat)
    dat = data_df %>% filter(clustering_method==clus_method, platform==plat) %>% select(data,met,GeneSet) %>% spread(key=GeneSet,value=met)
    
    for(j in 1:length(GS_test)){
      print(j)
      res = perm_ci_median(group1=dat[[GS_test[[j]][1]]], group2=dat[[GS_test[[j]][2]]], R=n_perm, alpha=0.05, n_grids=300) 
      forest_df = bind_rows(forest_df,data.frame(metric=metric,labeltext=paste(GS_test[[j]][1],"vs",GS_test[[j]][2],sep=" "),mean=res[1],est=format(res[1],scientific=TRUE,digits=3),high_prop=round(mean(dat[[GS_test[[j]][1]]]>=dat[[GS_test[[j]][2]]]),2),lower=res[2],upper=res[3],group=plat,method="paired_permutation"))
    }
  }
  
  ## draw forest plot
  for(i in 1:length(platform_vec)){
    plat = platform_vec[i]
    #png(fs::path(save_path,clus_method,paste(clus_method,metric,"med_platform_forest_CI",sep="_"),ext="png"),width=platform_width,height=platform_height,res=300,unit="in")
    fp=forest_df |> filter(group==plat) |> 
      forestplot(labeltext = c(labeltext, high_prop), xlab = "difference in median", title=plat, lwd.zero=2) |>
      fp_add_header(high_prop = "Proportion") |>
      fp_set_style(box = cols[1],
                   line = gpar(col=cols[2],lwd=2),
                   summary = cols[3],
                   zero = gpar(col="red",lty=3),
                   txt_gp = fpTxtGp(label = list(gpar(fontfamily = "", cex=2),
                                                 gpar(fontfamily = "", cex=2)),
                                    ticks = gpar(fontfamily = "", cex = 1.5),
                                    xlab  = gpar(fontfamily = "", cex = 1.5))) 
    
    png(fs::path(save_path,clus_method,paste(clus_method,metric,plat,"med_platform_forest_CI",sep="_"),ext="png"),width=platform_width,height=platform_height,res=300,unit="in")
    print(fp)
    dev.off()
  }
  
  return(forest_df)
}
