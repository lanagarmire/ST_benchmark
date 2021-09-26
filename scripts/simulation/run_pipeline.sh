#! /bin/bash

job_num0=$(sbatch --parsable sim_set_parameters.sbat)
job_num1=$(sbatch --parsable --dependency=afterok:$job_num0 sim_generate_data.sbat)
job_num2=$(sbatch --parsable --dependency=afterok:$job_num1 sim_scVI.sbat)
job_num3=$(sbatch --parsable --dependency=afterok:$job_num1 sim_SNF.sbat)
job_num4=$(sbatch --parsable --dependency=afterok:$job_num1 sim_getMOFA+_factor.sbat)
job_num5=$(sbatch --parsable --dependency=afterok:$job_num4 sim_getMOFA+.sbat)
job_num6=$(sbatch --parsable --dependency=afterok:$job_num1 sim_getfull.sbat)
job_num7=$(sbatch --parsable --dependency=afterok:$job_num1 sim_getWNN.sbat)
job_num8=$(sbatch --parsable --dependency=afterok:$job_num1 sim_getCIMLR.sbat)
job_num9=$(sbatch --parsable --dependency=afterok:$job_num2:$job_num3:$job_num5:$job_num6 sim_clus_bin_ground_truth.sbat)
job_num10=$(sbatch --parsable --dependency=afterok:$job_num7:$job_num8:$job_num9 sim_clus_bin_ground_truth_metrics_all.sbat)

# show dependencies in squeue output:
squeue -u $USER -o "%.8A %.4C %.10m %.20E"
