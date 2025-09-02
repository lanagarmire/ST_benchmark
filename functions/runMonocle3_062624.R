############################################
######## building master function for the entire pipeline
############################################

## retain the distance matrices 
## retain the transformed matrices
## cluster concensus plot
### cluster stability

runMonocle3 = function(gobject, n_CT, geneset_vec, num_nn, seed_vec, nPCs, pearson_residuals, kNN_dist_metric,clus_method,
                     init_match_thres, init_round_thres, init_res_start, init_res_end, init_res_length,
                     match_thres, round_thres, res_length,
                     metrics_vec, save_dir, run_clustering, run_label_matching, run_diagnostics, data_name, fname_tail,
                     pt_size){
  for(GeneSet in geneset_vec){
    M3_path = fs::path(save_dir,GeneSet)
    dir.create(M3_path,recursive=TRUE,showWarnings = TRUE)
    
    gs = retrieve_genesets(gobject,GeneSet)
    ngs = length(gs)
    
    if(nPCs=="elbow"){
      var_ls = compute_var_pc(gobject,name=paste("pca",GeneSet,sep="_"))
      elbow = PCAtools::findElbowPoint(var_ls$var_PC)
      ndims = elbow
    }else{
      ndims = nPCs
    }
    
    if(run_clustering){
      print(paste("running Monocle3 on",GeneSet,sep=" "))
      
      ## run Monocle3
      Monocle3_res_init = getCluster_monocle3_init(gene_expr=gobject@raw_exprs,reduce = FALSE,red_method="PCA",pcs=min(ndims,ngs),red_dim_df=gobject@dimension_reduction$cells$pca[[paste("pca",GeneSet,sep="_")]]$coordinates,
                                                   kNN_dist_metric = kNN_dist_metric, data_type=GeneSet, 
                                                   match_thres=init_match_thres, round_thres=init_round_thres, res_start=init_res_start, res_end=init_res_end, res_by=(init_res_end-init_res_start)/init_res_length,
                                                   clus_method=clus_method, seed_vec=seed_vec[1], CT_bin_size=2,num_iter=1,num_neighbors=num_nn, n_CT=n_CT)
      Monocle3_res = getCluster_monocle3(gene_expr=gobject@raw_exprs,reduce = FALSE, metadata=gobject@cell_metadata, red_method="PCA", pcs=min(ndims,ngs), red_dim_df=gobject@dimension_reduction$cells$pca[[paste("pca",GeneSet,sep="_")]]$coordinates,
                                         kNN_dist_metric = kNN_dist_metric, data_type=GeneSet,
                                         match_thres=match_thres, round_thres=round_thres, res_start=Monocle3_res_init$res_start, res_end=Monocle3_res_init$res_end, res_by=(Monocle3_res_init$res_end-Monocle3_res_init$res_start)/res_length, 
                                         clus_method=clus_method, seed_vec=seed_vec, num_iter=1,num_neighbors=num_nn, n_CT=n_CT)
      saveRDS(Monocle3_res, file = fs::path(M3_path,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3",sep="_"),fname_tail),ext="RDS"))
    }else{
      print(paste("loading existing Monocle3 clustering results on",GeneSet,sep=" "))
      Monocle3_res = readRDS(file = fs::path(M3_path,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3",sep="_"),fname_tail),ext="RDS"))
    }
    
    if(run_label_matching){
      print("obtaining majority Monocle3 label")
      ## cluster_ID_df = as.data.frame(Monocle3_res$GO@cell_metadata)
      ## rownames(cluster_ID_df) = Monocle3_res$GO@cell_metadata$cell_ID
      cluster_ID_df = as.data.frame(Monocle3_res$clus_res_df[,which(grepl(colnames(Monocle3_res$clus_res_df),pattern="M3_clus")),with=FALSE])
      rownames(cluster_ID_df) = Monocle3_res$clus_res_df$cell_ID
      n_clus_vec = sapply(paste("M3_clus",GeneSet,Monocle3_res$res,"seed",Monocle3_res$seed,sep="_"),FUN=function(i){mets = Monocle3_res$clus_res_df%>%pull(i)%>%table()%>%length()})
      match_ID = which(n_clus_vec == n_CT)
      ## cluster_ID_df = cluster_ID_df[,unname(match_ID),with=FALSE]
      cluster_ID_df = cluster_ID_df[,match_ID]
      ## cluster_ID_df = as.data.frame(cluster_ID_df)
      
      if(length(match_ID)==1){
        similarity_score = c(1)
      }else{
        similarity_score = sapply(1:ncol(cluster_ID_df),FUN=function(x){tmp_df=cluster_ID_df[,-x]; labs=cluster_ID_df[,x]; mean(apply(tmp_df,2,FUN=function(y){AMI(y,labs)}))})
      }
      
      cluster_ID_df_tmp = cluster_ID_df[,which.max(similarity_score),drop=FALSE]
      cluster_ID_df_tmp$cell_ID=rownames(cluster_ID_df_tmp)
      label_df = merge(cluster_ID_df_tmp,Monocle3_res$clus_res_df %>% dplyr::select(cell_ID,cell_type),by="cell_ID",sort=FALSE)
      ## label_df = data.frame(cell_ID=gobject@cell_ID, clus=cluster_ID_df[,which.max(similarity_score)], cell_type=gobject@cell_metadata$cell_type)
      clus_pick = colnames(cluster_ID_df)[which.max(similarity_score)]
      colnames(label_df)[2] = paste("M3",GeneSet,"clus",sep="_")
      print(str(label_df))
      
      print("matching labels")
      colnames(label_df) = c("cell_ID","clus","group")
      if(pearson_residuals){
        match_res = label_match(label_df,gobject@custom_expr,paste("M3",GeneSet,"clus",sep="_"))
        ## match_res = label_match(label_df,gobject@norm_expr,paste("Monocle3",GeneSet,"clus",sep="_"))
      }else{
        match_res = label_match(label_df,gobject@norm_expr,paste("M3",GeneSet,"clus",sep="_"))
      }
      
      label_df = match_res$clus_df
      colnames(label_df)[which(colnames(label_df)=="group")] = "cell_type"
      print(head(label_df))
      
      majority_label = list(similarity=similarity_score, matched_labels=cluster_ID_df, label=label_df, match_res=match_res)
      saveRDS(majority_label, file=fs::path(M3_path,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3","majority",sep="_"),fname_tail),ext="RDS"))
    }else{
      print("loading existing majority Monocle3 clustering label")
      majority_label = readRDS(file=fs::path(M3_path,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3","majority",sep="_"),fname_tail),ext="RDS"))
      label_df = majority_label$label
      print(head(label_df))
    }
    gobject@cell_metadata = merge(gobject@cell_metadata,label_df[,c(1,2)],by="cell_ID",sort=FALSE)
    
    ## compute metrics
    if(length(metrics_vec)!=0){
      print("computing metrics")
      metrics=data.frame(AMI=double(),n_clus=double(),ASC=double(),CHI=double(),SC=double(),VN=double(),GeneSet=character())
      metrics[1,]=NA
      
      if("supervised" %in% metrics_vec){
        res_AMI = AMI(label_df[,paste0("M3_",GeneSet,"_clus")], label_df[,"cell_type"])
        metrics$AMI = res_AMI
      }
      if("unsupervised" %in% metrics_vec){
        res_unsup_ls = compute_silhouette_v3(as.data.frame(label_df[,c(1,2)]),Monocle3_res,"monocle3",M3_dist=kNN_dist_metric)
        metrics$CHI = res_unsup_ls$CHI
        metrics$ASC = res_unsup_ls$silhouette
        metrics$DBI = res_unsup_ls$DBI
        metrics$dunn = res_unsup_ls$dunn
        metrics$pearsongamma = res_unsup_ls$pearsongamma
        metrics$ASC_scaled = res_unsup_ls$silhouette_scaled
      }
      if("spatial metrics" %in% metrics_vec){
        spat_df = merge(gobject@spatial_locs, gobject@cell_metadata%>%dplyr::select(c(cell_ID,cell_type,entro)), by="cell_ID",sort=FALSE)
        dist_df = spat_df %>% dplyr::select(c(sdimx,sdimy)) %>% dist() %>% as.matrix()
        width = as.numeric(quantile(dist_df[upper.tri(dist_df,diag=FALSE)],0.01))
        
        spat_met_ls = spat_metrics(spat_df,label_df,paste0("M3_",GeneSet,"_clus"), width)
        spat_ami_df = spat_df
        spat_ami_df$spat_ami = spat_met_ls$spat_ami
        metrics$SC = spat_met_ls$SC
        metrics$VN = spat_met_ls$VN
        metrics$mean_spat_ami = mean(spat_met_ls$spat_ami)
      }
      if("F1" %in% metrics_vec){
        F1_ls = compute_cluster_specific_F1(gobject,ref="cell_type",lab=paste("M3",GeneSet,"clus",sep="_"))
        metrics$weighted_F1=F1_ls$weighted_F1
        gobject@cell_metadata$mis_classified = gobject@cell_metadata$cell_type != gobject@cell_metadata[[paste("M3",GeneSet,"clus",sep="_")]]
      }
      metrics$n_clus=n_CT
      metrics$GeneSet=GeneSet
      metrics$ndims = ndims
      metrics$PC_elbow = nPCs=="elbow"
      metrics$knn_dist = kNN_dist_metric
      
      metrics_ls = list(scalar_metrics=metrics,spat_ami=spat_ami_df,cluster_F1=F1_ls)
      
      #saveRDS(metrics,file=fs::path(M3_path,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3","metrics",sep="_"),fname_tail),ext="RDS"))
      saveRDS(metrics_ls,file=fs::path(M3_path,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3","metrics",sep="_"),fname_tail),ext="RDS"))
    }
    
    ## compute diagnostics
    if(run_diagnostics){
      ## PCA and sNN plots
      print("running diagnostics")
      dir.create(fs::path(M3_path,"diagnostics",data_name),recursive=T,showWarnings = F)
      #pdf(fs::path(M3_path,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3","diagnostics",sep="_"),fname_tail),ext="pdf"))
      
      ## PCA heatmaps
      cols = pal_ucscgb(alpha=1)(26)
      PCA_reoder_ls = lapply(unique(gobject@cell_metadata$cell_type),FUN=function(x){idx=which(gobject@cell_metadata$cell_type==x);names(idx)=rep(x,length(idx));return(idx)})
      dist_mat = as.matrix(dist(gobject@dimension_reduction$cells$pca[[paste("pca",GeneSet,sep="_")]][["coordinates"]][,1:min(ndims,ngs),drop=FALSE]))
      png(fs::path(M3_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3","PCA_heatmap",sep="_"),fname_tail),ext="png"))
      PC_dist_heat(dist_mat,PCA_reoder_ls,cols[1:n_CT],name=paste("pca",GeneSet,sep="_"))
      dev.off()
      
      ## label comparison on tSNE
      gobject_tmp = gobject
      tsne_pca_mat = gobject_tmp@dimension_reduction$cells$pca[[paste("pca",GeneSet,sep="_")]][["coordinates"]][,1:min(ndims,ngs)]
      nCells_tsne_pca_mat = nrow(tsne_pca_mat)
      ## remove duplicates for the sake of visualization
      tsne_pca_mat = as.data.frame(tsne_pca_mat) %>% distinct() %>% as.matrix()
      if(nrow(tsne_pca_mat) < nCells_tsne_pca_mat){
        print(paste("removing",-nrow(tsne_pca_mat)+nCells_tsne_pca_mat,"duplicate row(s) in GeneSet PCA before computing tSNE.",sep=" "))
        gobject_tmp = subsetGiotto(gobject,cell_ids=rownames(tsne_pca_mat))
      }
      ## gobject_tmp@dimension_reduction$cells$pca[[paste("pca",GeneSet,sep="_")]][["coordinates"]] = tsne_pca_mat
      ## gobject_tmp = subsetGiotto(gobject,cell_ids=rownames(tsne_pca_mat))
      gobject_tmp = runtSNE(gobject_tmp,dim_reduction_name=paste("pca",GeneSet,sep="_"),dimensions_to_use=1:min(ndims,ngs),name=paste("tsne",GeneSet,ndims,sep="_"))
      
      cols_code = cols[1:n_CT]
      names(cols_code)=unique(gobject_tmp@cell_metadata$cell_type)
      dir.create(fs::path(M3_path,"diagnostics",data_name),recursive=TRUE,showWarnings = TRUE)
      
      tsne_df = gobject_tmp@dimension_reduction$cells$tsne[[paste("tsne",GeneSet,ndims,sep="_")]]$coordinates %>% as.data.frame()
      tsne_df = tsne_df %>% mutate(cell_ID=rownames(tsne_df))
      tsne_df = merge(tsne_df, gobject_tmp@cell_metadata, by="cell_ID",sort=FALSE)
      tsne_df[[paste0("M3_",GeneSet,"_clus")]] = factor(tsne_df[[paste0("M3_",GeneSet,"_clus")]])
      tsne_df$cell_type = factor(tsne_df$cell_type)
      
      p_tsne_cell_type = tsne_df %>% ggplot(aes_string(x="Dim.1",y="Dim.2",color="cell_type")) + geom_point(size=pt_size) + scale_color_manual(values=cols_code) + guides(color = guide_legend(title="",override.aes = list(size = pt_size+2))) + ggtitle("cell type") + theme_classic()
      p_tsne_clus = tsne_df %>% ggplot(aes_string(x="Dim.1",y="Dim.2",color=paste0("M3_",GeneSet,"_clus"))) + geom_point(size=pt_size) + scale_color_manual(values=cols_code) + guides(color = guide_legend(title="",override.aes = list(size = pt_size+2)))+ ggtitle(paste0("M3_",GeneSet,"_clus")) + theme_classic() 
      p_tsne_clus_highlights = tsne_df %>% ggplot(aes_string(x="Dim.1",y="Dim.2",color=paste0("M3_",GeneSet,"_clus"))) + geom_point(size=pt_size,aes(alpha = ifelse(gobject_tmp@cell_metadata$mis_classified == 1, 1, 0.1))) + 
        scale_color_manual(values=cols_code) + scale_alpha_identity() + guides(color = guide_legend(title="",override.aes = list(size = pt_size+2)))+ ggtitle(paste0("M3_",GeneSet,"_clus")) + theme_classic() 
      
      ggsave(plot=p_tsne_cell_type, filename = fs::path(M3_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3","tsne_cell_type",sep="_"),fname_tail),ext="png"),width=7,height=7,units="in",dpi=300)
      ggsave(plot=p_tsne_clus, filename = fs::path(M3_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3","tsne",GeneSet,sep="_"),fname_tail),ext="png"),width=7,height=7,units="in",dpi=300)
      ggsave(plot=p_tsne_clus_highlights, filename = fs::path(M3_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3","tsne",GeneSet,"highlight",sep="_"),fname_tail),ext="png"),width=7,height=7,units="in",dpi=300)
      
      ## label comparison on spatial coordinates
      spatial_df = merge(gobject_tmp@spatial_locs,gobject@cell_metadata,by="cell_ID",merge=FALSE)
      spatial_df[[paste0("M3_",GeneSet,"_clus")]] = factor(spatial_df[[paste0("M3_",GeneSet,"_clus")]])
      spatial_df$cell_type = factor(spatial_df$cell_type)
      
      p_spat_cell_type = spatial_df %>% ggplot(aes_string(x="sdimx",y="sdimy",color="cell_type")) + geom_point(size=pt_size) + scale_color_manual(values=cols_code) + guides(color = guide_legend(title="",override.aes = list(size = pt_size+2))) + ggtitle("cell type") + theme_classic()
      p_spat_clus = spatial_df %>% ggplot(aes_string(x="sdimx",y="sdimy",color=paste0("M3_",GeneSet,"_clus"))) + geom_point(size=pt_size) + scale_color_manual(values=cols_code) + guides(color = guide_legend(title="",override.aes = list(size = pt_size+2)))+ ggtitle(paste0("M3_",GeneSet,"_clus")) + theme_classic() 
      p_spat_clus_highlights = spatial_df %>% ggplot(aes_string(x="sdimx",y="sdimy",color=paste0("M3_",GeneSet,"_clus"))) + geom_point(size=pt_size,aes(alpha = ifelse(gobject_tmp@cell_metadata$mis_classified == 1, 1, 0.1))) + 
        scale_color_manual(values=cols_code) + scale_alpha_identity() + guides(color = guide_legend(title="",override.aes = list(size = pt_size+2)))+ ggtitle(paste0("M3_",GeneSet,"_clus")) + theme_classic() 
      
      ggsave(plot=p_spat_cell_type, filename = fs::path(M3_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3","spat_cell_type",sep="_"),fname_tail),ext="png"),width=6,height=6,units="in",dpi=300)
      ggsave(plot=p_spat_clus, filename = fs::path(M3_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3","spat",GeneSet,sep="_"),fname_tail),ext="png"),width=6,height=6,units="in",dpi=300)
      ggsave(plot=p_spat_clus_highlights, filename = fs::path(M3_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3","spat",GeneSet,"highlight",sep="_"),fname_tail),ext="png"),width=6,height=6,units="in",dpi=300)
      
      ## F1
      F1_heat = F1_ls$confusion_matrix %>% as.data.frame() %>% ggplot(aes(Var1,Var2,fill=Freq))+geom_tile()+xlab("ground truth")+ylab(paste0("Monocle3: ",GeneSet))+scale_fill_gradient(low="white", high="blue")+theme_classic()
      F1_bar = F1_ls$F1 %>% mutate(cluster=factor(cluster))%>% ggplot(aes(x=cluster,y=F1,fill=cluster))+geom_bar(stat="identity")+scale_fill_manual(values=cols_code)+ylim(c(0,1))+theme_classic()
      recall_bar = F1_ls$F1 %>% mutate(cluster=factor(cluster))%>% ggplot(aes(x=cluster,y=recall,fill=cluster))+geom_bar(stat="identity")+scale_fill_manual(values=cols_code)+ylim(c(0,1))+theme_classic()
      precision_bar = F1_ls$F1 %>% mutate(cluster=factor(cluster))%>% ggplot(aes(x=cluster,y=precision,fill=cluster))+geom_bar(stat="identity")+scale_fill_manual(values=cols_code)+ylim(c(0,1))+theme_classic()
      
      ggsave(plot=F1_heat, filename = fs::path(M3_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3",GeneSet,"confusion_matched",sep="_"),fname_tail),ext="png"),width=6,height=6,units="in",dpi=300)
      ggsave(plot=F1_bar, filename = fs::path(M3_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3",GeneSet,"cluster_F1",sep="_"),fname_tail),ext="png"),width=6,height=6,units="in",dpi=300)
      ggsave(plot=precision_bar, filename = fs::path(M3_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3",GeneSet,"cluster_precision",sep="_"),fname_tail),ext="png"),width=6,height=6,units="in",dpi=300)
      ggsave(plot=recall_bar, filename = fs::path(M3_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3",GeneSet,"cluster_recall",sep="_"),fname_tail),ext="png"),width=6,height=6,units="in",dpi=300)
      
      ## spatial AMI
      p_spat_ami = metrics_ls$spat_ami %>% ggplot(aes(x=sdimx,y=sdimy,color=spat_ami)) + geom_point(size=pt_size) + scale_color_gradient(low="white",high="blue") + theme_classic()
      ggsave(plot=p_spat_ami, filename = fs::path(M3_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3",GeneSet,"spatial_ami",sep="_"),fname_tail),ext="png"),width=6,height=6,units="in",dpi=300)
      
      p_spat_ami_rev = metrics_ls$spat_ami %>% ggplot(aes(x=sdimx,y=sdimy,color=spat_ami)) + geom_point(size=pt_size) + scale_color_gradient(low="blue",high="white") + theme_classic()
      ggsave(plot=p_spat_ami_rev, filename = fs::path(M3_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"Monocle3",GeneSet,"spatial_ami_rev",sep="_"),fname_tail),ext="png"),width=6,height=6,units="in",dpi=300)
    }
  }
}

#############################
### testing script
### module load Bioinformatics
### module load Rmonocle3
#############################
## library(monocle3)
## library(dplyr)
## library(matchingR)
## library(aricode)
## library(fpc)
## library(clusterSim)
## library(ggsci)
## library(ggplot2)
## library(gplots)
## library(RColorBrewer)
## library(Giotto,lib.loc = "/sw/pkgs/med/Rgiotto/1.1.0")
## 
## source("/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/functions/retrieve_genesets.R")
## source("/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/functions/M3_clus_res_grid_search.R")
## source("/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/functions/label_match_func.R")
## source("/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/functions/compute_silhouette_v3.R")
## source("/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/functions/spat_metric_01142024.R")
## source("/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/functions/spat_metric_03192021.R")
## source("/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/functions/helper_funcs.R")
## 
## data_name = "ST_MOB1"
## dat = readRDS(fs::path("/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/GB_rev_4/real_data_analysis_trouble_shoot/SCT_adjusted/feature_selection",paste(data_name,"rare_CT_removed_pre",sep="_"),ext="RDS"))
## 
## n_CT = length(unique(dat@cell_metadata$cell_type))
## 
## gobject=dat
## geneset_vec=c("hvg_low")
## GeneSet=geneset_vec[1]
## num_nn=15
## seed_vec=12345+c(0:9)
## nPCs="elbow" 
## pearson_residuals=TRUE
## kNN_dist_metric="euclidean"
## clus_method="leiden"
## init_match_thres=1 
## init_round_thres=5 
## init_res_start=0.01 
## init_res_end=1
## init_res_length=10
## match_thres=3
## round_thres=3 
## res_length=100
## metrics_vec=c("supervised","unsupervised","spatial metrics","purity") 
## save_dir="/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/GB_rev_4/real_data_analysis_trouble_shoot/SCT_adjusted/results/cluster/M3"
## run_clustering=TRUE 
## run_label_matching=TRUE
## run_diagnostics=TRUE
## fname_tail=""
## pt_size=1.5
