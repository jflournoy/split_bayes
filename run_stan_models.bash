#!/bin/bash
#
#SBATCH --job-name=run_stan_model
#SBATCH --output=run_stan_model_%A_%a.log
#
#SBATCH --cpus-per-task=28
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2000
#SBATCH --partition=defq,fat,long,longfat

models=(run_stan_ht_model.r run_stan_dl_model.r run_stan_pu_model.r)

module load R gcc

srun Rscript --verbose ${models[$SLURM_ARRAY_TASK_ID]}
