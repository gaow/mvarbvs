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

# Run the main DSC
./fixed_mix.dsc --host dsc_mnm.yml -o mnm_20200510 -s existing -e ignore &> mnm_20200510.log
# Run the mthess DSC
#./fixed_mix.dsc --target mthess -o mthess_20200526 -s existing -e ignore -c 38 --n_dataset 200 &> mnm_20200526.log
# Run data generation step only to get data for ED estimate
#./fixed_mix.dsc --host dsc_mnm.yml --target simulate_only --n_dataset 25000 -o mnm_sumstats -s existing -e ignore &> mnm_sumstats.log
