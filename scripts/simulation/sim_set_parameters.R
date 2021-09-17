rm(list=ls())

setwd("/home/liyijun/ST_benchmark_01082020")
library(data.table)

args = commandArgs(trailingOnly = TRUE)
data_path = as.character(args[1])
spat_l = as.numeric(args[2])
spat_r = as.numeric(args[3])
n_spat_prob = as.numeric(args[4])
n_reps = as.numeric(args[5])

spat_prob_vec = seq(spat_l, spat_r, length.out = n_spat_prob)
rep_id = c(1:n_reps)

params_df = expand.grid(spat_prob = spat_prob_vec, rep_id = rep_id)
head(params_df)

write.csv(params_df, file = fs::path(data_path, "sim_params_df", ext="csv"))
