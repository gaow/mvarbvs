#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --partition=mstephens
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=dsc_ukb
#SBATCH --mail-type=BEGIN,END,FAIL

# load module
module load R/4.1.0

# Simulate the main DSC
# dsc ukb_rss_main.dsc --target simulate_only -o output/ukb_rss_20220619 -s existing --host dsc_mnm.yml &> mvsusie_rss_ukb.log

# Run the main DSC
# dsc ukb_rss_main.dsc --host dsc_mnm.yml -o output/ukb_rss_20220619 -s existing &> mvsusie_rss_cs_ukb.log

# Run DSC for smaller number of traits 2, indep signal
# dsc ukb_rss_small.dsc --target simulate_small2_indep_only --host dsc_mnm.yml -o output/ukb_rss_small2_indep_20220619 -s existing -e ignore &> mnm_rss_ukb_small2_indep.log
# dsc ukb_rss_small.dsc --target small2_indep_compare --host dsc_mnm.yml -o output/ukb_rss_small2_indep_20220619 -s existing -e ignore &> mnm_rss_ukb_small2_indep.log

# Run DSC for cafeh main simulation
# dsc ukb_rss_main.dsc --target cafeh --host dsc_mnm.yml -o output/ukb_rss_20220619 -s existing &> cafeh_rss_cs_ukb.log


