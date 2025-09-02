############################################
######## building master function for the entire pipeline
############################################

runKmeans = function(gobject, n_CT, pearson_residuals, batch=FALSE, batchsize, geneset_vec, dist_vec, nPCs, metrics_vec, save_dir, run_clustering, run_label_matching, run_diagnostics, data_name, fname_tail, pt_size){
  for(GeneSet in geneset_vec){
    for(dist in dist_vec){
      KM_path = fs::path(save_dir,GeneSet,dist)
      dir.create(KM_path,recursive=T,showWarnings = F)
      
      print(paste("running kmeans on",GeneSet,"using",dist,sep=" "))
      ## run kmeans
      PCA_mat = gobject@dimension_reduction$cells$pca[[paste("pca",GeneSet,sep="_")]][["coordinates"]]
      if(nPCs=="elbow"){
        PC_var_ls = compute_var_pc(gobject,name=paste("pca",GeneSet,sep="_"))
        PC_elbow = PCAtools::findElbowPoint(PC_var_ls$var_PC)
        ndims = PC_elbow
      }else{
        ndims = nPCs
      }
      
      if(run_clustering){        
        KM_res = getCluster_kmeans(gobject@norm_expr, GeneSet, reduced=T, PCA_mat = PCA_mat, pc_dims=min(ndims,ncol(PCA_mat)), dist = dist, cl0=n_CT, batch = batch, batch_size = batchsize)
        saveRDS(KM_res,file=fs::path(KM_path,paste0(paste(data_name,"rare_CT_removed",GeneSet,"kmeans_res",sep="_"),fname_tail),ext="RDS"))
      }else{
        print("loading existing kmeans results")
        KM_res = readRDS(fs::path(KM_path,paste0(paste(data_name,"rare_CT_removed",GeneSet,"kmeans_res",sep="_"),fname_tail),ext="RDS"))
      }
      
      ## extract results
      ## label matching
      if(run_label_matching){
        label_df = merge(KM_res$cluster_label, gobject@cell_metadata %>% dplyr::select(cell_ID,cell_type),by="cell_ID",sort=F)
        print(str(label_df))
        
        print("matching labels")
        colnames(label_df) = c("cell_ID","clus","group")
        if(pearson_residuals){
          match_res = label_match(label_df,gobject@custom_expr,paste("KM",GeneSet,dist,"clus",sep="_"))
        }else{
          match_res = label_match(label_df,gobject@norm_expr,paste("KM",GeneSet,dist,"clus",sep="_"))
        }
        label_df = match_res$clus_df
        colnames(label_df)[which(colnames(label_df)=="group")] = "cell_type"
        print(head(label_df))
        
        majority_label = list(label=label_df, match_res=match_res)
        saveRDS(majority_label, file=fs::path(KM_path,paste0(paste(data_name,"rare_CT_removed",GeneSet,"kmeans_labels",sep="_"),fname_tail),ext="RDS"))
      }else{
        print("loading existing kmeans clustering label")
        majority_label = readRDS(file=fs::path(KM_path,paste0(paste(data_name,"rare_CT_removed",GeneSet,"kmeans_labels",sep="_"),fname_tail),ext="RDS"))
        label_df = majority_label$label
      }
      
      ## compute metrics
      if(length(metrics_vec)!=0){
        print("computing clustering metrics")
        metrics=data.frame(AMI=double(),n_clus=double(),ASC=double(),CHI=double(),SC=double(),VN=double(),GeneSet=character(),dist=double())
        metrics[1,]=NA
        
        if("supervised" %in% metrics_vec){
          res_AMI = AMI(label_df[,paste0("KM_",GeneSet,"_",dist,"_clus")], label_df[,"cell_type"])
          metrics$AMI = res_AMI
        }
        if("unsupervised" %in% metrics_vec){
          res_unsup_ls = compute_silhouette_v3(clus_df=as.data.frame(label_df[,c(1,2)]),clus_res=KM_res,clus_method="kmeans",ndims=ndims,red_dim=PCA_mat)
          metrics$CHI = res_unsup_ls$CHI
          metrics$ASC = res_unsup_ls$silhouette
          metrics$DBI = res_unsup_ls$DBI
        metrics$ASC_scaled = res_unsup_ls$silhouette_scaled
        metrics$dunn = res_unsup_ls$dunn
        metrics$pearsongamma = res_unsup_ls$pearsongamma
        }
        if("spatial metrics" %in% metrics_vec){
          spat_df = merge(gobject@spatial_locs, gobject@cell_metadata%>%dplyr::select(c(cell_ID,cell_type,entro)), by="cell_ID",sort=FALSE)
          dist_df = spat_df %>% dplyr::select(c(sdimx,sdimy)) %>% dist() %>% as.matrix()
          width = as.numeric(quantile(dist_df[upper.tri(dist_df,diag=FALSE)],0.01))
          
          spat_met_ls = spat_metrics(spat_df,label_df,paste0("KM_",GeneSet,"_",dist,"_clus"), width)
          spat_ami_df = spat_df
          spat_ami_df$spat_ami = spat_met_ls$spat_ami
          metrics$SC = spat_met_ls$SC
          metrics$VN = spat_met_ls$VN
          metrics$mean_spat_ami = mean(spat_met_ls$spat_ami)
        }
        if("F1" %in% metrics_vec){
          gobject@cell_metadata = merge(gobject@cell_metadata,label_df[,c(1,2)],by="cell_ID",sort=FALSE)
          F1_ls = compute_cluster_specific_F1(gobject,ref="cell_type",lab=paste("KM",GeneSet,dist,"clus",sep="_"))
          metrics$weighted_F1=F1_ls$weighted_F1
          gobject@cell_metadata$mis_classified = gobject@cell_metadata$cell_type != gobject@cell_metadata[[paste("KM",GeneSet,dist,"clus",sep="_")]]
        }
        
        metrics$n_clus=n_CT
        metrics$GeneSet=GeneSet
        metrics$dist=dist
        metrics$ndims = ndims
        metrics$PC_elbow = nPCs=="elbow"
        
        metrics_ls = list(scalar_metrics=metrics,spat_ami=spat_ami_df,cluster_F1=F1_ls)
        
        saveRDS(metrics_ls,file=fs::path(KM_path,paste0(paste(data_name,"rare_CT_removed",GeneSet,"kmeans","metrics",sep="_"),fname_tail),ext="RDS"))
      }
      
      ## compute diagnostics
      if(run_diagnostics){
        print("running diagnostics")
        dir.create(fs::path(KM_path,"diagnostics",data_name),recursive=TRUE,showWarnings = TRUE)
        
        ## compare organized PCA plots
        cols = pal_ucscgb(alpha=1)(26)
        cols_code = cols[1:n_CT]
        names(cols_code)=unique(gobject@cell_metadata$cell_type)
        
        PCA_reoder_ls = lapply(unique(gobject@cell_metadata$cell_type),FUN=function(x){idx=which(gobject@cell_metadata$cell_type==x);names(idx)=rep(x,length(idx));return(idx)})
        dist_mat = as.matrix(KM_res$celldist)
        png(fs::path(KM_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"kmeans","diagnostics",sep="_"),fname_tail),ext="png"))
        PC_dist_heat(dist_mat,PCA_reoder_ls,cols_code,name=paste("pca",GeneSet,sep="_"))
        dev.off()
        
        ## label comparison on tSNE
        gobject_tmp = gobject
        tsne_pca_mat = gobject_tmp@dimension_reduction$cells$pca[[paste("pca",GeneSet,sep="_")]][["coordinates"]][,1:min(ndims,ncol(PCA_mat))]
        nCells_tsne_pca_mat = nrow(tsne_pca_mat)
        ## remove duplicates for the sake of visualization
        tsne_pca_mat = as.data.frame(tsne_pca_mat) %>% distinct() %>% as.matrix()
        if(nrow(tsne_pca_mat) < nCells_tsne_pca_mat){
          print(paste("removing",-nrow(tsne_pca_mat)+nCells_tsne_pca_mat,"duplicate row(s) in GeneSet PCA before computing tSNE.",sep=" "))
          gobject_tmp = subsetGiotto(gobject,cell_ids=rownames(tsne_pca_mat))
        }
        ## gobject_tmp@dimension_reduction$cells$pca[[paste("pca",GeneSet,sep="_")]][["coordinates"]] = tsne_pca_mat
        ## gobject_tmp = subsetGiotto(gobject,cell_ids=rownames(tsne_pca_mat))
        gobject_tmp = runtSNE(gobject_tmp,dim_reduction_name=paste("pca",GeneSet,sep="_"),dimensions_to_use=1:min(ndims,ncol(PCA_mat)),name=paste("tsne",GeneSet,ndims,sep="_"))
        
        tsne_df = gobject_tmp@dimension_reduction$cells$tsne[[paste("tsne",GeneSet,ndims,sep="_")]]$coordinates %>% as.data.frame()
        tsne_df = tsne_df %>% mutate(cell_ID=rownames(tsne_df))
        tsne_df = merge(tsne_df, gobject_tmp@cell_metadata, by="cell_ID",sort=FALSE)
        tsne_df[[paste0("KM_",GeneSet,"_",dist,"_clus")]] = factor(tsne_df[[paste0("KM_",GeneSet,"_",dist,"_clus")]])
        tsne_df$cell_type = factor(tsne_df$cell_type)
        
        p_tsne_cell_type = tsne_df %>% ggplot(aes_string(x="Dim.1",y="Dim.2",color="cell_type")) + geom_point(size=pt_size) + scale_color_manual(values=cols_code) + guides(color = guide_legend(title="",override.aes = list(size = pt_size+2))) + ggtitle("cell type") + theme_classic()
        p_tsne_clus = tsne_df %>% ggplot(aes_string(x="Dim.1",y="Dim.2",color=paste0("KM_",GeneSet,"_",dist,"_clus"))) + geom_point(size=pt_size) + scale_color_manual(values=cols_code) + guides(color = guide_legend(title="",override.aes = list(size = pt_size+2)))+ ggtitle(paste0("KM_",GeneSet,"_",dist,"_clus")) + theme_classic() 
        p_tsne_clus_highlights = tsne_df %>% ggplot(aes_string(x="Dim.1",y="Dim.2",color=paste0("KM_",GeneSet,"_",dist,"_clus"))) + geom_point(size=pt_size,aes(alpha = ifelse(gobject_tmp@cell_metadata$mis_classified == 1, 1, 0.1))) + 
          scale_color_manual(values=cols_code) + scale_alpha_identity() + guides(color = guide_legend(title="",override.aes = list(size = pt_size+2)))+ ggtitle(paste0("KM_",GeneSet,"_",dist,"_clus")) + theme_classic() 
        
        ggsave(plot=p_tsne_cell_type, filename = fs::path(KM_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"kmeans","tsne_cell_type",sep="_"),fname_tail),ext="png"),width=7,height=7,units="in",dpi=300)
        ggsave(plot=p_tsne_clus, filename = fs::path(KM_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"kmeans","tsne",GeneSet,sep="_"),fname_tail),ext="png"),width=7,height=7,units="in",dpi=300)
        ggsave(plot=p_tsne_clus_highlights, filename = fs::path(KM_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"kmeans","tsne",GeneSet,"highlight",sep="_"),fname_tail),ext="png"),width=7,height=7,units="in",dpi=300)
        
        ## label comparison on spatial coordinates
        spatial_df = merge(gobject_tmp@spatial_locs,gobject@cell_metadata,by="cell_ID",merge=FALSE)
        spatial_df[[paste0("KM_",GeneSet,"_",dist,"_clus")]] = factor(spatial_df[[paste0("KM_",GeneSet,"_",dist,"_clus")]])
        spatial_df$cell_type = factor(spatial_df$cell_type)
        
        p_spat_cell_type = spatial_df %>% ggplot(aes_string(x="sdimx",y="sdimy",color="cell_type")) + geom_point(size=pt_size) + scale_color_manual(values=cols_code) + guides(color = guide_legend(title="",override.aes = list(size = pt_size+2))) + ggtitle("cell type") + theme_classic()
        p_spat_clus = spatial_df %>% ggplot(aes_string(x="sdimx",y="sdimy",color=paste0("KM_",GeneSet,"_",dist,"_clus"))) + geom_point(size=pt_size) + scale_color_manual(values=cols_code) + guides(color = guide_legend(title="",override.aes = list(size = pt_size+2)))+ ggtitle(paste0("KM_",GeneSet,"_",dist,"_clus")) + theme_classic() 
        p_spat_clus_highlights = spatial_df %>% ggplot(aes_string(x="sdimx",y="sdimy",color=paste0("KM_",GeneSet,"_",dist,"_clus"))) + geom_point(size=pt_size,aes(alpha = ifelse(gobject_tmp@cell_metadata$mis_classified == 1, 1, 0.1))) + 
          scale_color_manual(values=cols_code) + scale_alpha_identity() + guides(color = guide_legend(title="",override.aes = list(size = pt_size+2)))+ ggtitle(paste0("KM_",GeneSet,"_",dist,"_clus")) + theme_classic() 
        
        ggsave(plot=p_spat_cell_type, filename = fs::path(KM_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"kmeans","spat_cell_type",sep="_"),fname_tail),ext="png"),width=6,height=6,units="in",dpi=300)
        ggsave(plot=p_spat_clus, filename = fs::path(KM_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"kmeans","spat",GeneSet,sep="_"),fname_tail),ext="png"),width=6,height=6,units="in",dpi=300)
        ggsave(plot=p_spat_clus_highlights, filename = fs::path(KM_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"kmeans","spat",GeneSet,"highlight",sep="_"),fname_tail),ext="png"),width=6,height=6,units="in",dpi=300)
        
        ## spatial AMI
        p_spat_ami = metrics_ls$spat_ami %>% ggplot(aes(x=sdimx,y=sdimy,color=spat_ami)) + geom_point(size=pt_size) + scale_colour_gradientn(colours = terrain.colors(10)) + theme_classic()
        ggsave(plot=p_spat_ami, filename = fs::path(KM_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"kmeans",GeneSet,"spatial_ami",sep="_"),fname_tail),ext="png"),width=6,height=6,units="in",dpi=300)
        
        p_spat_ami_rev = metrics_ls$spat_ami %>% ggplot(aes(x=sdimx,y=sdimy,color=spat_ami)) + geom_point(size=pt_size) + scale_colour_gradientn(colours = terrain.colors(10)) + theme_classic()
        ggsave(plot=p_spat_ami_rev, filename = fs::path(KM_path,"diagnostics",data_name,paste0(paste(data_name,"rare_CT_removed",GeneSet,"kmeans",GeneSet,"spatial_ami_rev",sep="_"),fname_tail),ext="png"),width=6,height=6,units="in",dpi=300)
      }
    }
  }
}

#############################
### testing script
### module load R/4.2.0
### module load gsl
#############################
## rm(list=ls())
## 
## library(data.table)
## library(dplyr)
## library(tidyr)
## library(Matrix)
## library(matchingR)
## library(proxy)
## library(Rfast)
## library(aricode)
## library(DescTools)
## library(clusterSim)
## library(fpc)
## library(RColorBrewer)
## library(ggsci)
## library(gplots)
## library(ggplot2)
## library(Giotto,lib.loc = "/sw/pkgs/med/Rgiotto/1.1.0")
## 
## source("/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/functions/compute_silhouette_v3.R")
## source("/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/functions/spat_metric_03192021.R")
## source("/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/functions/spat_metric_01142024.R")
## source("/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/functions/getCluster_kmeans_092423.R")
## source("/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/functions/helper_funcs.R")
## source("/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/functions/label_match_func.R")
## 
## data_name = "ST_MOB1"
## dat = readRDS(fs::path("/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/GB_rev_4/real_data_analysis_trouble_shoot/SCT_adjusted/feature_selection",paste(data_name,"rare_CT_removed_pre",sep="_"),ext="RDS"))
## 
## n_CT = length(unique(dat@cell_metadata$cell_type))
## gobject = dat 
## batch=FALSE
## batchsize=NULL
## geneset_vec=c("hvg_low")
## pearson_residuals=TRUE
## GeneSet=geneset_vec[1]
## dist_vec=c("pearson") 
## dist=dist_vec[1]
## nPCs="elbow"
## metrics_vec=c("supervised","unsupervised","spatial metrics","F1")
## save_dir="/nfs/turbo/umms-lgarmire/liyijun/ST_benchmark_01082020_revision/GB_rev_4/real_data_analysis_trouble_shoot/SCT_adjusted/results/cluster/kmeans"
## run_clustering=TRUE
## run_label_matching=TRUE
## run_diagnostics=TRUE
## fname_tail=""
## pt_size=1.5
  