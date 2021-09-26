import pandas as pd
import leidenalg as la
import igraph as ig
import snf
import sklearn
import numpy as np
import parser
import argparse
import os

def do_snf(hvg_file_name, svg_file_name, fuse_save_name, num_nn=30, metric="sqeuclidean"):

  hvg_df = pd.read_csv(hvg_file_name)
  svg_df = pd.read_csv(svg_file_name)

  if "Unnamed: 0" in hvg_df.columns:
    hvg_df = hvg_df.drop(columns = ["Unnamed: 0"])

  if "Unnamed: 0" in svg_df.columns:
    svg_df = svg_df.drop(columns = ["Unnamed: 0"])

  cell_names = hvg_df.columns.values

  hvg_mat = hvg_df.values.T
  svg_mat = svg_df.values.T

  hvg_affinity = snf.make_affinity(hvg_mat,K=num_nn,metric = metric)
  svg_affinity = snf.make_affinity(svg_mat,K=num_nn,metric = metric)

  fused = snf.snf(hvg_affinity,svg_affinity,K=num_nn)

  pd.DataFrame(data=fused,columns = cell_names).to_csv(fuse_save_name,index = False)

  return fused



def do_leiden_clustering(fused_mat, desired_num_clusters, label_save_name, current_num_clusters):

  g = ig.Graph.Weighted_Adjacency(fused_mat.tolist())

  resolution = 0.5
  num_clusters = current_num_clusters
  num_iters = 1

  num_cells = fused_mat.shape[0]

  while num_clusters != desired_num_clusters and num_iters < 100 and resolution <= 1 and resolution > 0:

    partition = la.find_partition(g,la.RBConfigurationVertexPartition,weights = g.es["weight"],resolution_parameter = resolution,seed=1234, n_iterations=-1)

    num_clusters = len(partition)

    if num_clusters > desired_num_clusters:
        resolution -= 0.1/num_iters
    elif num_clusters < desired_num_clusters:
        resolution += 0.1/num_iters

    num_iters += 1

  cluster_labels = np.zeros(num_cells)

  for i in range(len(partition)):
    cluster_members = partition[i]
    for member in cluster_members:
        cluster_labels[member] = i

  pd.Series(data = cluster_labels,name = "x").to_csv(label_save_name,index=False)

  return cluster_labels, partition



def main(job_id, data_path, data_ref, diff_in_cov, svg_pval_thres, num_nn, ndims, scVI_error, SNF_metric, WNN_smooth, CIMLR_k, save_name):
  os.chdir("/home/liyijun/ST_benchmark_01082020_re_1")
  os.chdir(data_path)

  'load parameters'
  params_name = f"{data_ref}/real_data_names.csv"
  params_df = pd.read_csv(params_name, dtype = str).drop(columns="Unnamed: 0")
  method_proc = pd.read_csv('/home/liyijun/ST_benchmark_01082020_re_1/method_proc.csv', index_col = 0)
  method_name = "SNF"
  
  if save_name == "diff_in_cov":
      save_path = f"{data_ref}/{save_name}_{diff_in_cov}"
  elif save_name == "svg_pval_thres":
      save_path = f"{data_ref}/{save_name}_{svg_pval_thres}"
  elif save_name == "ndims":
      save_path = f"{data_ref}/{save_name}_{ndims}"
  elif save_name == "scVI_error":
      save_path = f"{data_ref}/{save_name}_{scVI_error}"
  elif save_name == "SNF_metric":
      save_path = f"{data_ref}/{save_name}_{SNF_metric}"
  elif save_name == "WNN_smooth":
      save_path = f"{data_ref}/{save_name}_{WNN_smooth}"
  elif save_name == "CIMLR_k":
      save_path = f"{data_ref}/{save_name}_{CIMLR_k}"

  data_name = params_df.iloc[job_id]['name']
  num_clusters = int(params_df.iloc[job_id]['cl0'])
  hvg_name = f"{save_path}/results/HVG/{data_name}_HVG_{method_proc.at[method_name,'proc']}.csv"
  svg_name = f"{save_path}/results/SG/{data_name}_SG_bin_{method_proc.at[method_name,'proc']}.csv"

  if not os.path.exists(f"{save_path}/results/cluster/{method_name}"):
    os.makedirs(f"{save_path}/results/cluster/{method_name}", exist_ok=True)

  if not os.path.exists(f"{save_path}/results/{method_name}"):
    os.makedirs(f"{save_path}/results/{method_name}", exist_ok=True)

  ncl = 0
  n_neighbors = num_nn
  iter0 = 1
  while ncl != num_clusters and n_neighbors > 0 and iter0 < 15:
    if ncl != 0:
      if ncl > num_clusters:
          if n_neighbors <= 10:
            n_neighbors = n_neighbors + 1
          else:
            n_neighbors = n_neighbors + 5          
      elif ncl < num_clusters:
          if n_neighbors <= 10:
            n_neighbors = n_neighbors - 1
          else:
            n_neighbors = n_neighbors - 5
    
    print(n_neighbors)      
    fused_mat = do_snf(hvg_file_name = hvg_name,svg_file_name = svg_name, num_nn = n_neighbors, fuse_save_name = f"{save_path}/results/{method_name}/{data_name}_SNF.csv")
    labels, partition = do_leiden_clustering(fused_mat, label_save_name = f"{save_path}/results/cluster/{method_name}/{data_name}_SNF_bin_LC.csv",desired_num_clusters = num_clusters, current_num_clusters = ncl)
    ncl = len(partition)
    iter0 = iter0 + 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='provide arguments')
    parser.add_argument('job_id', type=int, help='')
    parser.add_argument('data_path', type=str, help='')
    parser.add_argument('data_ref', type=str, help='')
    parser.add_argument('diff_in_cov', type=float, help='')
    parser.add_argument('svg_pval_thres', type=float, help='')
    parser.add_argument('num_nn',type=int,help='')
    parser.add_argument('ndims',type=int, help='')
    parser.add_argument('scVI_error', type=str, help='')
    parser.add_argument('SNF_metric', type=str, help='')
    parser.add_argument('WNN_smooth', type=int, help='')
    parser.add_argument('CIMLR_k', type=int, help='')
    parser.add_argument('save_name', type=str, help='')
    args = parser.parse_args()
    main(args.job_id-1, args.data_path, args.data_ref, args.diff_in_cov, args.svg_pval_thres, args.num_nn, args.ndims, args.scVI_error, args.SNF_metric, args.WNN_smooth, args.CIMLR_k, args.save_name)
