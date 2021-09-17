# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 15:45:07 2021

@author: stefs
"""

import numpy as np
import pandas as pd
#import plotnine as p9

import scanpy as sc
import scvi

import anndata

import itertools

import sklearn.metrics as skm
import argparse
import os
import faulthandler
import torch
import time

faulthandler.enable()
#from tqdm.notebook import tqdm
from tqdm import tqdm

sc.set_figure_params(figsize=(4, 4))

def train_model(file_name_svg, file_name_hvg, n_latent, desired_num_clusters, seed, num_nn = 30, default_resolution = 0.8, n_hidden = 128, dropout = 0.1, gene_likelihood = "zinb", latent_distribution = "normal", error_measure = "reconstruction"):

  raw_svg_data = pd.read_csv(file_name_svg)
  raw_hvg_data = pd.read_csv(file_name_hvg)

  if "Unnamed: 0" in raw_hvg_data.columns:
    raw_hvg_data = raw_hvg_data.drop(columns = ["Unnamed: 0"])

  if "Unnamed: 0" in raw_svg_data.columns:
    raw_svg_data = raw_svg_data.drop(columns = ["Unnamed: 0"])

  combined_data = np.concatenate((raw_hvg_data.values,raw_svg_data.values))
  adata = anndata.AnnData(X = combined_data.T)

  torch.cuda.empty_cache()
  scvi.settings.seed = seed

  scvi.data.setup_anndata(adata)
  vae = scvi.model.SCVI(adata, n_hidden = n_hidden, n_latent = n_latent, gene_likelihood = gene_likelihood, latent_distribution = latent_distribution,dropout_rate = dropout)
  vae.train()
  adata.obsm["X_scVI"] = vae.get_latent_representation()
  adata.obsm["X_normalized_scVI"] = vae.get_normalized_expression()

  if error_measure == "silhouette":
    num_clusters = None
    resolution = default_resolution

    number_its = 1
    while num_clusters != desired_num_clusters and number_its < 50:

      sc.pp.neighbors(adata, n_neighbors=num_nn, use_rep="X_scVI")
      sc.tl.umap(adata, min_dist=0.2)
      sc.tl.leiden(adata, key_added="leiden_scVI",resolution=resolution, random_state=1234)

      num_clusters = len(set(adata.obs["leiden_scVI"]))

      if num_clusters > desired_num_clusters:
        resolution -= 0.4 / number_its

      if num_clusters < desired_num_clusters:
        resolution += 0.4 / number_its

      print("In the loop, resolution is ",resolution," and the number of clusters is ",num_clusters)

      number_its += 1

    metric = skm.silhouette_score(adata.obsm["X_scVI"],adata.obs["leiden_scVI"])

    print("Silhouette score is ",metric)

  if error_measure == "reconstruction":
    metric = -1.0*vae.get_reconstruction_error()["reconstruction_loss"]

  return adata,metric


def grid_search(file_name_svg, file_name_hvg, desired_num_clusters, seed, list_latent = [], list_hidden = [], list_dropouts = [0.2], gene_likelihood = ["zinb"], latent_distribution = ["normal"], default_resolution = 0.8, num_nn = 30, error_measure="reconstruction"):

  best_metric = -1.0
  best_adata = None
  best_nlatent = None
  best_nhidden = None
  best_gene_likely = None
  best_latent_dist = None

  comb_list = [list_latent,list_hidden,list_dropouts,gene_likelihood,latent_distribution]

  for n_latent,n_hidden,dropout,gene_likely,latent_dist in tqdm(list(itertools.product(*comb_list))):

    print("Initializing the model with n_latent = ",n_latent," n_hidden = ",n_hidden," n_dropout = ",dropout," gene_likelihood = ",gene_likely," and latent_distribution = ",latent_dist)
    adata, metric = train_model(file_name_svg = file_name_svg, file_name_hvg = file_name_hvg, num_nn = num_nn, default_resolution = default_resolution, n_hidden = n_hidden, n_latent = n_latent, dropout = dropout, gene_likelihood = gene_likely, latent_distribution = latent_dist, desired_num_clusters = desired_num_clusters, error_measure = error_measure, seed = seed)
    if metric > best_metric:
      best_metric = metric
      best_adata = adata
      best_nlatent = n_latent
      best_nhidden = n_hidden
      best_dropout = dropout
      best_gene_likely = gene_likely
      best_latent_dist = latent_dist

  return best_metric, best_adata, best_nlatent, best_nhidden, best_dropout,best_gene_likely,best_latent_dist
  


def main(job_id, data_path, ndims, scVI_error):

    print("The job ID is ",job_id)

    'set working directory'
    os.chdir("/home/liyijun/ST_benchmark_01082020")
    save_path = data_path

    'load parameters'
    params_name = f"{save_path}/data/cells.csv"
    params_df = pd.read_csv(params_name, dtype = str).drop(columns="Unnamed: 0")
    method_proc = pd.read_csv('/home/liyijun/ST_benchmark_01082020/method_proc.csv', index_col = 0)
    method_name = "scVI"
    
    data_name = params_df.iloc[job_id]['cells']
    rep = params_df.iloc[job_id]["reps"]
    num_clusters = 8
    hvg_name = f"{save_path}/data/data_cells_{data_name}_{method_proc.at[method_name,'proc']}_1_{rep}.csv"
    svg_name = f"{save_path}/data/data_cells_{data_name}_{method_proc.at[method_name,'proc']}_2_{rep}.csv"

    if ndims == 0:
      ndims = [5,10,15,20,25]
    else:
      ndims = [ndims]

    if not os.path.exists(f"{save_path}/results/{method_name}"):
      os.makedirs(f"{save_path}/results/{method_name}", exist_ok=True)

    'stochastic runs'
    sd = 12345
    tic = time.perf_counter()
    metric, adata, nlatent, nhidden, dropout, best_gene_likely, best_latent_dist = grid_search(list_latent = ndims, list_hidden = [128], list_dropouts = [0.1], file_name_svg = svg_name, file_name_hvg = hvg_name, desired_num_clusters = num_clusters, num_nn=30, error_measure = scVI_error, seed = sd)
    toc = time.perf_counter()
    elapse = toc-tic
    with open(f"{save_path}/results/{method_name}/{data_name}_{method_name}_time_{rep}.csv", 'w') as f:      
      f.write("%d" % elapse)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='provide arguments')
    parser.add_argument('job_id', type=int, help='')
    parser.add_argument('data_path', type=str, help='')
    parser.add_argument('ndims', type=int, help='')
    parser.add_argument('scVI_error', type=str, help='')
    args = parser.parse_args()
    #main(args.job_id-1, args.data_path, args.ndims, args.scVI_error)
    main(args.job_id*7-1, args.data_path, args.ndims, args.scVI_error)
