#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=dsc_ukb_eigen
#SBATCH --mail-type=BEGIN,END,FAIL

# load the openmpi module
module load R

# Run the main DSC
<<<<<<< HEAD
# dsc ukb_rss.dsc --host dsc_mnm.yml -o output/ukb_rss -s existing &> mnm_rss_ukb.log
# Run the mthess DSC
#./fixed_mix.dsc --host dsc_mnm.yml --target mthess -o mthess_20200526 -s existing -e ignore --n_dataset 200 &> mnm_20200526.log
=======
dsc ukb_rss.dsc --host dsc_mnm.yml -o output/ukb_rss -s existing &> mnm_rss_ukb.log
>>>>>>> 5a9635b0c369a0a8c58d976d2b6e91fe650a9c3c
# Run data generation step only to get data for ED estimate
# dsc ukb_rss.dsc --host dsc_mnm.yml --target simulate_only -s existing -o output/ukb_rss &> mnm_rss_ukb.log

dsc ukb_rss.dsc --host dsc_mnm.yml --target eigen -s existing -o output/ukb_rss_eigen &> mnm_rss_ukb_eigen.log
