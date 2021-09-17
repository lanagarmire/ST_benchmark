from mofapy2.run.entry_point import entry_point
import pandas as pd
import numpy as np
import math 
import argparse

### load job ID
parser = argparse.ArgumentParser(description='provide arguments')
parser.add_argument('job_id', type=int, help='')
parser.add_argument('data_path', type=str, help='')
parser.add_argument('thres',type=float,help='')
parser.add_argument('data_ref', type=str, help='')
parser.add_argument('diff_in_cov', type=float, help='')
parser.add_argument('svg_pval_thres', type=float, help='')
parser.add_argument('ndims',type=int, help='')
parser.add_argument('scVI_error', type=str, help='')
parser.add_argument('SNF_metric', type=str, help='')
parser.add_argument('WNN_smooth', type=int, help='')
parser.add_argument('CIMLR_k', type=int, help='')
parser.add_argument('save_name', type=str, help='')
parser.add_argument('n_random_reps', type=int, help='')
args = parser.parse_args()

job_id =  args.job_id-1
data_path = args.data_path
thres = args.thres
data_ref = args.data_ref
diff_in_cov = args.diff_in_cov
svg_pval_thres = args.svg_pval_thres
ndims = args.ndims
scVI_error = args.scVI_error
SNF_metric = args.SNF_metric
WNN_smooth = args.WNN_smooth
CIMLR_k = args.CIMLR_k
save_name = args.save_name
n_random_reps = args.n_random_reps

'set save'
if save_name == "thres":
  save_path = f"{data_ref}/{save_name}_{thres}"
elif save_name == "diff_in_cov":
  save_path = f"{data_ref}/{save_name}_{diff_in_cov}"
elif save_name == "svg_pval_thres":
  save_path = f"{data_ref}/{save_name}_{svg_pval_thres}"
elif save_name == "num_nn":
  save_path = f"{data_ref}/{save_name}_{num_nn}"
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

### load simulation parameters   #### change
import os
os.chdir("/home/liyijun/ST_benchmark_01082020/")
os.chdir(data_path) 
sim_params_df = pd.read_csv('sim_params_df.csv', index_col=0)
method_proc = pd.read_csv('/home/liyijun/ST_benchmark_01082020/method_proc.csv', index_col = 0)
method_name = "MOFA+"
if sim_params_df.iloc[job_id]['spat_prob'] == 1.0:
  spat_lim = int(sim_params_df.iloc[job_id]['spat_prob'])
  rep_id = int(sim_params_df.iloc[job_id]['rep_id'])
  data_name = "sim_spat_" + str(int(sim_params_df.iloc[job_id]['spat_prob'])) + "_rep_" + str(int(sim_params_df.iloc[job_id]['rep_id']))
else:
  spat_lim = sim_params_df.iloc[job_id]['spat_prob']
  rep_id = int(sim_params_df.iloc[job_id]['rep_id'])
  data_name = "sim_spat_" + str(sim_params_df.iloc[job_id]['spat_prob']) + "_rep_" + str(int(sim_params_df.iloc[job_id]['rep_id']))

ent = entry_point()

### read input files ### change
file1 = f"{save_path}/results/SG/sim_spat_{spat_lim}_rep_{rep_id}_SG_bin_{method_proc.at[method_name,'proc']}.csv" 
file2 = f"{save_path}/results/HVG/sim_spat_{spat_lim}_rep_{rep_id}_HVG_{method_proc.at[method_name,'proc']}.csv" 

dat1 = pd.read_csv(file1, index_col=0)
dat2 = pd.read_csv(file2, index_col=0)

col_names = []
for i in range(dat1.shape[1]):
  col_names.append('cell_'+str(i+1))

dat1.columns = col_names
dat2.columns = col_names

dat1.reset_index(inplace=True)
dat2.reset_index(inplace=True)

df1 = pd.wide_to_long(dat1, stubnames = ['cell_'], i = 'index', j ='number')
df2 = pd.wide_to_long(dat2, stubnames = ['cell_'], i = 'index', j ='number')

df1.reset_index(inplace=True)
df2.reset_index(inplace=True)

df1.columns = ['feature','sample','value']
df2.columns = ['feature','sample','value']

df1['view'] = "dat1"
df2['view'] = "dat2"
df = df1.append(df2)
df['group'] = "group0"

ent.set_data_options(scale_groups = False,scale_views = True)
ent.set_data_df(df)

#### set reduced dimension ### change
if ndims==0:
  ndims_f = f"{save_path}/results/scVI/sim_spat_{spat_lim}_rep_{rep_id}_scVI_hyperparams.csv"
  ndims_df = pd.read_csv(ndims_f)
  bottleneck_nodes = ndims_df["n_latent"][0]
else:
  bottleneck_nodes = ndims

ent.set_model_options(factors = bottleneck_nodes, spikeslab_weights = True)

'stochastic runs'
seed_list = [x+12345 for x in list(range(n_random_reps))]
for sd in seed_list:
  ent.set_train_options(iter = 500, convergence_mode = "fast",startELBO = 1, freqELBO = 1, dropR2 = None,gpu_mode = True, verbose = False,seed = sd)
  ent.build()
  ent.run()

  #### save output ### change
  outfile = f"{save_path}/results/{method_name}/sim_spat_{spat_lim}_rep_{rep_id}_{method_name}_dim_{bottleneck_nodes}_bin_{sd}.hdf5"	##change 
  ent.save(outfile)


