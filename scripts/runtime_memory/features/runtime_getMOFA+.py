from mofapy2.run.entry_point import entry_point
import pandas as pd
import numpy as np
import math
import argparse
import time

### load job ID
parser = argparse.ArgumentParser(description='provide arguments')
parser.add_argument('job_id', type=int, help='')
parser.add_argument('data_path', type=str, help='')
parser.add_argument('ndims',type=int, help='')
args = parser.parse_args()

job_id =  args.job_id-1
data_path = args.data_path
ndims = args.ndims

### load simulation parameters   #### change
import os
os.chdir("/home/liyijun/ST_benchmark_01082020/")
os.chdir(data_path) 
params_name = f"data/features.csv"
params_df = pd.read_csv(params_name).drop(columns="Unnamed: 0")
method_proc = pd.read_csv('/home/liyijun/ST_benchmark_01082020/method_proc.csv', index_col = 0)
method_name = "MOFA+"
data_name = params_df["features"][job_id]
rep = params_df["reps"][job_id]
save_path = data_path

### MOFA+ start
### read input files ### change
file1 = f"data/data_features_{data_name}_{method_proc.at[method_name,'proc']}_1_{rep}.csv" 
file2 = f"data/data_features_{data_name}_{method_proc.at[method_name,'proc']}_2_{rep}.csv" 

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

#### set reduced dimension ### change
bottleneck_nodes = ndims

### define data options
sd=12345
tic = time.perf_counter()

ent = entry_point()
ent.set_data_options(scale_groups = False,scale_views = True)
ent.set_data_df(df)
ent.set_model_options(factors = bottleneck_nodes, spikeslab_weights = True)
ent.set_train_options(iter = 500, convergence_mode = "fast",startELBO = 1, freqELBO = 1, dropR2 = None,gpu_mode = True, verbose = False,seed = sd)
ent.build()
ent.run()

toc = time.perf_counter()
elapse = toc-tic

if not os.path.exists(f"results/{method_name}"):
  os.makedirs(f"results/{method_name}", exist_ok=True)

with open(f"results/{method_name}/{data_name}_{method_name}_time_{rep}.csv", 'w') as f:
  f.write("%d" % elapse)
