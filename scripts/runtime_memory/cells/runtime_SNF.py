import pandas as pd
import leidenalg as la
import igraph as ig
import snf
import sklearn
import numpy as np
import parser
import argparse
import os
import time

def do_snf(hvg_file_name, svg_file_name, save_name, num_nn=30, metric="sqeuclidean"):

  hvg_df = pd.read_csv(hvg_file_name)
  svg_df = pd.read_csv(svg_file_name)

  if "Unnamed: 0" in hvg_df.columns:
    hvg_df = hvg_df.drop(columns = ["Unnamed: 0"])

  if "Unnamed: 0" in svg_df.columns:
    svg_df = svg_df.drop(columns = ["Unnamed: 0"])

  cell_names = hvg_df.columns.values

  hvg_mat = hvg_df.values.T
  svg_mat = svg_df.values.T

  tic = time.perf_counter()
  hvg_affinity = snf.make_affinity(hvg_mat,K=num_nn,metric = metric)
  svg_affinity = snf.make_affinity(svg_mat,K=num_nn,metric = metric)

  fused = snf.snf(hvg_affinity,svg_affinity,K=num_nn)
  toc = time.perf_counter()
  elapse = toc-tic
  with open(save_name, 'w') as f:      
    f.write("%d" % elapse)

  return fused


def main(job_id, data_path):
  os.chdir("/home/liyijun/ST_benchmark_01082020")
  os.chdir(data_path)

  'load parameters'
  params_name = f"data/cells.csv"
  params_df = pd.read_csv(params_name, dtype = str).drop(columns="Unnamed: 0")
  method_proc = pd.read_csv('/home/liyijun/ST_benchmark_01082020/method_proc.csv', index_col = 0)
  method_name = "SNF"

  data_name = params_df.iloc[job_id]['cells']
  rep = params_df.iloc[job_id]['reps']
  hvg_name = f"data/data_cells_{data_name}_{method_proc.at[method_name,'proc']}_1_{rep}.csv"
  svg_name = f"data/data_cells_{data_name}_{method_proc.at[method_name,'proc']}_2_{rep}.csv"

  if not os.path.exists(f"results/{method_name}"):
    os.makedirs(f"results/{method_name}", exist_ok=True)

  n_neighbors = 30    
  save_name = f"results/{method_name}/{data_name}_{method_name}_{rep}_time.csv"
  fused_mat = do_snf(hvg_file_name = hvg_name,svg_file_name = svg_name, save_name = save_name, num_nn = n_neighbors)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='provide arguments')
    parser.add_argument('job_id', type=int, help='')
    parser.add_argument('data_path', type=str, help='')
    args = parser.parse_args()
    #main(args.job_id-1, args.data_path)
    main(args.job_id*7-1, args.data_path)
