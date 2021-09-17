rm(list=ls())

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

setwd("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/runtime_cells")
if(!dir.exists("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig5")){
  dir.create("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig5")
}

cells_df = read.csv("data/cells.csv", row.names = 1)
#cells_df = cells_df[-1,,drop=FALSE]
n_random_reps = 5
seed_list = 12345+c(1:n_random_reps)-1

cells_df = cells_df %>% filter(!(cells %in% c(150, 300)))

models = c("CIMLR", "concatenation", "MOFA+", "scVI", "SNF", "WNN")
model_labels = c("CIMLR", "concatenation", "MOFA+", "scVI", "SNF", "Seurat v4")
num_cells = time = rep = method = c()

###################### runtime #######################
for(i in 1:nrow(cells_df)){
  for(j in 1:length(models)){
    num_cells = c(num_cells, cells_df$cells[i])
    rep = c(rep, cells_df$reps[i])
    method = c(method, models[j]) 
 
    if(models[j] %in% c("CIMLR", "concatenation", "WNN")){
      load(fs::path("results", models[j], paste(cells_df$cells[i], models[j], cells_df$reps[i], "time", sep="_"), ext="RData"))
      time = c(time, elapse)
    }else if(models[j] %in% c("MOFA+", "SNF")){
      time_df = read.csv(file = fs::path("results", models[j], paste(cells_df$cells[i], models[j], cells_df$reps[i], "time", sep="_"), ext="csv"), header = FALSE)
      time = c(time, time_df[1,1])
    }else if(models[j] %in% c("scVI")){
      time_df = read.csv(file = fs::path("results", models[j], paste(cells_df$cells[i], models[j], "time", cells_df$reps[i], sep="_"), ext="csv"), header = FALSE)
      time = c(time, time_df[1,1])
    }
  }
}

runtime_cells_df = data.frame(num_cells = num_cells, rep = rep, method = method, runtime = time)
runtime_cells_df = runtime_cells_df %>% dplyr::group_by(num_cells, method) %>% dplyr::summarise(runtime = mean(runtime))
runtime_cells_df = runtime_cells_df %>% dplyr::group_by(num_cells) %>% dplyr::mutate(runtime_rank = rank(runtime, ties.method = "first"))

manual_col_pallete0 = c("WNN" = "aquamarine3", "concatenation" = "coral4", "SG"= "black", "MOFA+" = "blue", "CIMLR" = "blueviolet", "scVI" = "darkgoldenrod1", "SNF" = "darkgreen")

## line plot
runtime_cells_df %>% mutate(method = factor(method, levels=models)) %>%
ggplot(aes(x=num_cells, y=runtime, group = method, color = method))+
       geom_line() + geom_point()+
       scale_color_manual(values = manual_col_pallete0, breaks = models, labels = model_labels)+
       ylim(0, 350) +
       labs(x = "number of cells", y="", color = "model", title = "")+
       theme_bw()+
       theme(axis.text.x=element_text(size=13), legend.text=element_text(size=13), legend.title = element_text(size = 15), legend.position = "bottom")
ggsave(fs::path("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig5", "runtime_cells", ext="png"), width = 6, height = 6)

## heatmap
runtime_cells_df %>% mutate(method = factor(method, levels = models)) %>%
ggplot(aes(x = method, y = factor(num_cells, levels = c(4500, 3000, 1500, 1000, 500)), fill = runtime)) +
	geom_tile() +
	geom_text(aes(label = runtime_rank), size = 6)+
	scale_fill_viridis(discrete = F, option = "C", begin=0.25, direction=1)+
	scale_x_discrete(position = "top")+
	scale_y_discrete(labels = c("4500"="4500 cells", "3000"="3000 cells", "1500"="1500 cells", "1000"="1000 cells", "500"="500 cells"))+
	guides(colour=F)+
	labs(x="", y="", fill = "")+
	theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size=12), legend.position="none")
ggsave(fs::path("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig6", "runtime_cells_heatmap", ext="png"), width = 4.46, height = 2.32) #change this



############################ memory  ############################
memory_cells_df = data.frame(num_cells = unique(cells_df$cells), 
			     scVI = c(2205.696,	2209.792,	2224.128,	2252.8,	2291.712),
			     MOFAp = c(0,	0,	0,	493.6148,	1075.2), 
                             WNN = c(0,	0,	0,	0,	460.724),
			     SNF = c(0,	0,	0,	1120.256,	2373.632), 
                             CIMLR = c(662.006,	2105.344,	4392.96,	14974.976,	34777.088),
                             concatenation = c(0.00, 0.00, 0.00, 0.00, 0.00))
colnames(memory_cells_df)[3] = "MOFA+"

memory_cells_df = gather(memory_cells_df, method, memory, scVI:concatenation)
memory_cells_df = memory_cells_df %>% dplyr::group_by(num_cells) %>% dplyr::mutate(memory_rank = rank(memory, ties.method = "first"))

## line plots
memory_cells_df %>% mutate(method = factor(method, levels=models)) %>%
ggplot(aes(x=num_cells, y=memory, group = method, color = method))+
       geom_line() + geom_point()+
       scale_color_manual(values = manual_col_pallete0, breaks = models, labels = model_labels)+
       ylim(0, 35000) +
       labs(x = "number of cells", y="", color = "model", title = "")+
       theme_bw()+
       theme(axis.text.x=element_text(size=13), legend.text=element_text(size=13), legend.title = element_text(size = 15), legend.position = "bottom")
ggsave(fs::path("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig5", "memory_cells", ext="png"), width = 6, height = 6)

## heatmaps
memory_cells_df %>% mutate(method = factor(method, levels = models)) %>%
ggplot(aes(x = method, y = factor(num_cells, levels = c(4500, 3000, 1500, 1000, 500)), fill = memory)) +
	geom_tile() +
	geom_text(aes(label = memory_rank), size = 6)+
	scale_fill_viridis(discrete = F, option = "C", begin=0.25, direction=1)+
	scale_x_discrete(position = "top")+
	scale_y_discrete(labels = c("4500"="4500 cells", "3000"="3000 cells", "1500"="1500 cells", "1000"="1000 cells", "500"="500 cells"))+
	guides(colour=F)+
	labs(x="", y="", fill = "")+
	theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size=12), legend.position="none")
ggsave(fs::path("/home/liyijun/ST_benchmark_01082020/data_analysis/real_data/vis_manuscript/Fig6", "memory_cells_heatmap", ext="png"), width = 4.46, height = 2.32) #change this
