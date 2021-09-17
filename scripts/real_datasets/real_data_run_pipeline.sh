#! /bin/bash

job_num1=$(sbatch --parsable real_data_getHVG_SG.sbat)
job_num2=$(sbatch --parsable --dependency=afterok:$job_num1 real_data_scVI.sbat)
job_num3=$(sbatch --parsable --dependency=afterok:$job_num2 real_data_SNF.sbat)
job_num4=$(sbatch --parsable --dependency=afterok:$job_num3 real_data_getMOFA+_factor.sbat)
job_num5=$(sbatch --parsable --dependency=afterok:$job_num4 real_data_getMOFA+.sbat)
job_num6=$(sbatch --parsable --dependency=afterok:$job_num5 real_data_getfull.sbat)
job_num7=$(sbatch --parsable --dependency=afterok:$job_num6 real_data_getWNN.sbat)
job_num8=$(sbatch --parsable --dependency=afterok:$job_num7 real_data_getCIMLR.sbat)
job_num9=$(sbatch --parsable --dependency=afterok:$job_num8 real_data_clus_bin_ground_truth.sbat)
job_num10=$(sbatch --parsable --dependency=afterok:$job_num9 real_data_clus_bin_ground_truth_metrics_all.sbat)

# show dependencies in squeue output:
squeue -u $USER -o "%.8A %.4C %.10m %.20E"
