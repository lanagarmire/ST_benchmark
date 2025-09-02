spat_metric = function(label_df, comp=NULL, width=NULL, DEA_method=NULL, GO_object=NULL, clust_labels=NULL, markers_truth = NULL, method){
  if(method == "AIC"){
    library(spatstat)
    if(comp == "truth"){
      est_range = ripras(label_df$sdimx, label_df$sdimy,shape = "rectangle")
      
      dat_spat_truth = ppp(label_df$sdimx, label_df$sdimy,
                           xrange = est_range$xrange, yrange = est_range$yrange,
                           marks = factor(paste0("cl_",label_df$truth)))
      
      dat_spat_cl = ppp(label_df$sdimx, label_df$sdimy,
                        xrange = est_range$xrange, yrange = est_range$yrange,
                        marks = factor(paste0("cl_",label_df$cl)))
      
      fit_marks_truth = ppm(dat_spat_truth, ~marks)
      fit_marks_cl = ppm(dat_spat_cl, ~marks)
    }else if(comp == "null"){
      est_range = ripras(label_df$sdimx, label_df$sdimy,shape = "rectangle")
      
      dat_spat_cl = ppp(label_df$sdimx, label_df$sdimy,
                        xrange = est_range$xrange, yrange = est_range$yrange,
                        marks = factor(paste0("cl_",label_df$cl)))
      
      fit_marks_truth = ppm(dat_spat_cl, ~1)
      fit_marks_cl = ppm(dat_spat_cl, ~marks)
    }
    return(AIC(fit_marks_truth)-AIC(fit_marks_cl))
  }else if(method == "DBI"){
    library(clusterSim)
    #compute DBI
    DBI = index.DB(x = label_df[,c("sdimx","sdimy")],
                   cl = label_df$cl)
    return(DBI)
  }else if(method == "Var_norm"){
    #obtain group pairwise means
    library(dplyr)
    truth_labels = unique(label_df$truth)
    cl_labels = unique(label_df$cl)
    truth_avg_dist = cl_avg_dist = rep(NA,length(truth_labels))
    for(i in 1:length(truth_labels)){
      tmp1 = label_df %>%
        filter(truth == truth_labels[i]) %>%
        dplyr::select(sdimx,sdimy)%>%
        dist()
      if(length(tmp1) == 0){
        truth_avg_dist[i] = 0
      }else{
        truth_avg_dist[i] = mean(tmp1)
      }
      
      tmp2 = label_df %>%
        filter(cl == cl_labels[i]) %>%
        dplyr::select(sdimx,sdimy)%>%
        dist()
      if(length(tmp2) == 0){
        cl_avg_dist[i] = 0
      }else{
        cl_avg_dist[i] = mean(tmp2)
      }
    }
    #obtain density
    f = min(c(truth_avg_dist, cl_avg_dist))
    t = max(c(truth_avg_dist, cl_avg_dist))
    y_truth = density(truth_avg_dist, from = -10*f, to = 10*t)
    y_cl = density(cl_avg_dist, from = -10*f, to = 10*t)
    dx = diff(y_truth$x)[1]
    vn = 1-0.5*sum(abs(y_truth$y-y_cl$y))*dx
    
    return(vn)
  }else if(method == "spat_MI"){
    library(aricode)
    spat_MI = rep(NA,dim(label_df)[1])
    for(i in 1:(dim(label_df)[1])){
      pt_x = label_df$sdimx[i]
      pt_y = label_df$sdimy[i]
      label_df_tmp = label_df %>% filter(between(sdimx, left = pt_x-width, right = pt_x+width)) %>%
        filter(between(sdimy, pt_y-width, pt_y+width))
      spat_MI[i] = AMI(label_df_tmp$truth, label_df_tmp$cl)
    }
    label_df$spat_MI = spat_MI
    
    return(label_df)
  }else if(method == "wt_MI"){
    #compute weight matrix
    dist_mat = label_df %>% select(sdimx, sdimy) %>% dist()
    wt_mat = as.matrix(1.0/(dist_mat/sum(dist_mat)))
    #compute pxy, px, py
    tot = dim(label_df)[1]
    tab_xy = table(label_df$truth, label_df$cl)
    tab_truth = table(label_df$truth)
    tab_cl = table(label_df$cl)
    #compute weighted mutual information
    wtMI = 0
    for(i in 1:tot){
      for(j in 1:tot){
        pxy = tab_xy[as.character(label_df$truth[i]), as.character(label_df$cl[j])]/tot
        if(pxy != 0){
          px = tab_truth[as.character(label_df$truth[i])]/tot
          py = tab_cl[as.character(label_df$cl[i])]/tot
          #wtMI = wtMI + wt_mat[i,j]*pxy*log(pxy/(px*py))
          wtMI = wtMI + pxy*log(pxy/(px*py))
          print(pxy*log(pxy/(px*py)))
        }
      }
    }
    return(wtMI)
  }else if(method == "DEA"){
    library(Giotto)
    GO_object@cell_metadata$clust = clust_labels
    markers_df = findMarkers_one_vs_all(GO_object, method = DEA_method, min_genes = 0, cluster_column = "clust")
    markers = markers_df$genes

    #metric
    P = markers_truth
    TP = intersect(P, markers)
    FN = setdiff(P, TP)
    N = setdiff(GO_object@gene_ID, P)
    TN = intersect(N, setdiff(GO_object@gene_ID, markers))
    FP = setdiff(N, TN)

    TPR = length(TP)/length(P)
    PPV = length(TP)/(length(TP)+length(FP))
    FDR = 1-PPV
    FM = sqrt(PPV*TPR)
    CSI = length(TP)/(length(TP)+length(FN)+length(FP))
    ACC = (length(TP)+length(TN))/(length(TP)+length(TN)+length(FP)+length(FN))
    F1 = 2*length(TP)/(2*length(TP)+length(FP)+length(FN))

    return(c("TPR"=TPR, "PPV"=PPV, "FDR"=FDR, "FM"=FM, "CSI"=CSI, "ACC"=ACC, "F1"=F1)) 
  }
}

