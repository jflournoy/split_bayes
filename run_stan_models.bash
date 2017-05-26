#!/bin/bash
#
#SBATCH --job-name=run_stan_model
#SBATCH --output=run_stan_model.log
#
#SBATCH --cpus-per-task=28
#SBATCH --ntasks=3
#SBATCH --mem-per-cpu=1000
#SBATCH --partition=defq,fat,long,longfat

module load R gcc

srun Rscript --verbose run_stan_ht_model.r
srun Rscript --verbose run_stan_dl_model.r
srun Rscript --verbose run_stan_pu_model.r
