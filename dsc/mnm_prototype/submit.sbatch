#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=dsc-submitter
#SBATCH --output=dsc-log-%J.out
#SBATCH --error=dsc-log-%J.err

#./fixed_mix.dsc --host dsc_mnm.yml -o mnm_20200430 -s existing -e ignore &> mnm_20200430.log
./fixed_mix.dsc --host dsc_mnm.yml --target simulate_only --n_dataset 20000 -o mnm_sumstats -s existing -e ignore &> mnm_sumstats.log
