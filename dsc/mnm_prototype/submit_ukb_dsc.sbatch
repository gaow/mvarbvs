#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=dsc_ukb
#SBATCH --mail-type=BEGIN,END,FAIL

# load module
module load R/4.1.0

# Simulate the main DSC
dsc ukb_rss.dsc --target simulate_only -o output/ukb_rss_20220619 -s existing --host dsc_mnm.yml &> mvsusie_rss_ukb.log

# Run the main DSC
# dsc ukb_rss.dsc --host dsc_mnm.yml -o output/ukb_rss_20220429 -s existing &> mvsusie_rss_ukb.log

# Run DSC for smaller number of traits
# dsc ukb_rss.dsc --target small_compare --host dsc_mnm.yml -o output/ukb_rss_small -s existing &> mnm_rss_ukb_small.log

# Run DSC on smaller number of traits using MT-HESS
# dsc ukb_rss.dsc --target full_data --host dsc_mnm.yml -o output/ukb_rss_small --truncate -s existing &> mnm_rss_ukb_small_fulldata.log

# Run data generation step only to get data for ED estimate
# dsc ukb_rss.dsc --host dsc_mnm.yml --target simulate_only -s existing -o output/ukb_rss &> mnm_rss_ukb.log

