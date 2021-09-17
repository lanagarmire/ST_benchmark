rm(list=ls())

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

setwd("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/runtime_features")
if(!dir.exists("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig5")){
  dir.create("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig5")
}

features_df = read.csv("data/features.csv", row.names = 1)
n_random_reps = 5
seed_list = 12345+c(1:n_random_reps)-1

features_df = features_df %>% filter(!(features %in% c(150, 300)))

models = c("CIMLR", "concatenation", "MOFA+", "scVI", "SNF", "WNN")
model_labels = c("CIMLR", "concatenation", "MOFA+", "scVI", "SNF", "Seurat v4")
num_features = time = rep = method = c()

######################## runtime ####################
for(i in 1:nrow(features_df)){
  for(j in 1:length(models)){
    num_features = c(num_features, features_df$features[i])
    rep = c(rep, features_df$reps[i])
    method = c(method, models[j])  

    if(models[j] %in% c("CIMLR", "concatenation", "WNN")){
      load(fs::path("results", models[j], paste(features_df$features[i], models[j], features_df$reps[i], "time", sep="_"), ext="RData"))
      time = c(time, elapse)
    }else if(models[j] %in% c("SNF")){
      time_df = read.csv(file = fs::path("results", models[j], paste(features_df$features[i], models[j], features_df$reps[i], "time", sep="_"), ext="csv"), header = FALSE)
      time = c(time, time_df[1,1])
    }else if(models[j] %in% c("scVI", "MOFA+")){
      time_df = read.csv(file = fs::path("results", models[j], paste(features_df$features[i], models[j], "time", features_df$reps[i], sep="_"), ext="csv"), header = FALSE)
      time = c(time, time_df[1,1])      
    }
  }
}

runtime_features_df = data.frame(num_features = num_features, rep = rep, method = method, runtime = time)
runtime_features_df = runtime_features_df %>% dplyr::group_by(num_features, method) %>% dplyr::summarise(runtime = mean(runtime))
runtime_features_df = runtime_features_df %>% dplyr::group_by(num_features) %>% dplyr::mutate(runtime_rank = rank(runtime, ties.method = "first"))

manual_col_pallete0 = c("WNN" = "aquamarine3", "concatenation" = "coral4", "SG"= "black", "MOFA+" = "blue", "CIMLR" = "blueviolet", "scVI" = "darkgoldenrod1", "SNF" = "darkgreen")

## line plot
runtime_features_df %>% mutate(method = factor(method, levels=models)) %>%
ggplot(aes(x=num_features, y=runtime, color = method))+
       geom_line() + geom_point()+
       scale_color_manual(values = manual_col_pallete0, breaks = models, labels = model_labels)+
       ylim(0, 350) +
       labs(x = "number of features", y="", color = "model", title = "")+
       theme_bw()+
       theme(axis.text.x=element_text(size=13), legend.text=element_text(size=13), legend.title = element_text(size = 15), legend.position = "bottom")
ggsave(fs::path("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig5", "runtime_features", ext="png"), width = 6, height = 6)

## heatmap
runtime_features_df %>% mutate(method = factor(method, levels = models)) %>%
ggplot(aes(x=method, y = factor(num_features, levels = c(10000, 5000, 3000, 1000, 500)), fill=runtime))+
       geom_tile()+
       geom_text(aes(label = runtime_rank), size = 6)+
       scale_fill_viridis(discrete=F, option="C", begin=0.25, direction = 1)+
       scale_x_discrete(position="top")+
       scale_y_discrete(labels = c("10000"="10000 features", "5000"="5000 features", "3000"="3000 features", "1000"="1000 features", "500"="500 features"))+
       guides(colour=F)+
       labs(x="", y="", fill="")+
       theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size=12), legend.position = "none")
ggsave(fs::path("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig6", "runtime_features_heatmap", ext="png"), width = 4.83, height = 2.32) #change this



################## memory ##################
memory_features_df = data.frame(num_features = unique(features_df$features),
				scVI = c(804.864, 882.688, 901.12, 1396.252, 2398.208),
				SNF = c(0, 0, 0, 0, 0),
				MOFAp = c(0,	0,	0,	1177.6,	2170.88),
				concatenation = c(0, 0, 0, 0, 0),
				WNN = c(0,	0,	0,	0,	128.31),
				CIMLR = c(660.654,	660.378,	688.984,	718.464,	851.916))
colnames(memory_features_df)[4] = "MOFA+"

memory_features_df = gather(memory_features_df, method, memory, scVI:CIMLR)
memory_features_df = memory_features_df %>% dplyr::group_by(num_features) %>% dplyr::mutate(memory_rank = rank(memory, ties.method = "first"))

## line plots
memory_features_df %>% mutate(method = factor(method, levels=models)) %>%
  ggplot(aes(x=num_features, y=memory, group = method, color = method))+
  geom_line() + geom_point()+
  scale_color_manual(values = manual_col_pallete0, breaks = models, labels = model_labels)+
  ylim(0, 35000)+
  labs(x = "number of features", y="", color = "model", title = "")+
  theme_bw()+
  theme(axis.text.x=element_text(size=13), legend.text=element_text(size=13), legend.title = element_text(size = 15), legend.position = "bottom")
ggsave(fs::path("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig5", "memory_features", ext="png"), width = 6, height = 6)

## heatmaps
memory_features_df %>% mutate(method = factor(method, levels = models)) %>%
  ggplot(aes(x = method, y = factor(num_features, levels = c(10000, 5000, 3000,	1000, 500)), fill = memory)) +
  geom_tile() +
  geom_text(aes(label = memory_rank), size = 6)+
  scale_fill_viridis(discrete = F, option = "C", begin=0.25, direction=1)+
  scale_x_discrete(position = "top")+
  scale_y_discrete(labels = c("10000"="10000 features", "5000"="5000 features", "3000"="3000 features", "1000"="1000 features", "500"="500 features"))+
  guides(colour=F)+
  labs(x="", y="", fill = "")+
  theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size=12), legend.position="none")
ggsave(fs::path("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig6", "memory_features_heatmap", ext="png"), width = 4.83, height = 2.32) #change this
