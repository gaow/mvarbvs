#!/usr/bin/env dsc

%include modules/data_ukb
%include modules/simulate_fixed_mix
%include modules/mnm
%include modules/mnm_rss
%include modules/mnm_suff
%include modules/susie_suff
%include modules/susie_rss
%include modules/paintor
%include modules/score

DSC:
  define:
    simulate: artificial_mixture_ukb, ukb_bloodcells_mixture
    simulate_small: artificial_mixture_small, artificial_mixture_ukb_small2, artificial_mixture_ukb_small2_indep
    mnm: mnm_oracle, mnm_naive
    mnm_rss: mnm_rss_oracle, mnm_rss_naive, mnm_rss_ed, mnm_rss_ed_ddcan, mnm_rss_identity, mnm_rss_shared, mnm_rss_shared_corZ, mnm_rss_naive_corZ, mnm_rss_ed_corZ, mnm_rss_ed_ddcan_corZ, mnm_rss_identity_corZ, mnm_rss_shared_corY, mnm_rss_naive_corY, mnm_rss_ed_corY, mnm_rss_ed_ddcan_corY, mnm_rss_identity_corY
    mnm_suff: mnm_suff_oracle, mnm_suff_naive, mnm_suff_ed, mnm_suff_ed_ddcan, mnm_suff_identity
    susie_uni: susie_suff, susie_rss
  run:
    default: data_ukb * simulate * ((mnm_suff, mnm_rss) * mvsusie_scores, susie_uni * susie_scores)
    small_compare: data_ukb * simulate_small * ((mnm_suff_naive, mnm_rss_naive, mnm_rss_naive_corY, mnm_rss_naive_corZ) * mvsusie_scores, PAINTOR)
    simulate_small_only: data_ukb * simulate_small
    simulate_only: data_ukb * simulate 
  exec_path: code
  output: output/ukb_rss
  global:
    data_file: /project2/mstephens/yuxin/ukb-bloodcells/regions1000_5000.txt
    prior_file: /project2/mstephens/yuxin/ukb-bloodcells/ukb_prior_simulation.rds
    # nullz_file: /project2/mstephens/yuxin/ukb-bloodcells/nullz_cor.rds
    genotype_dir: '/project2/mstephens/yuxin/ukb-bloodcells/genotypes/'
    ld_dir: '/project2/mstephens/yuxin/ukb-bloodcells/LD/'
    ldeigen_dir: '/project2/mstephens/yuxin/ukb-bloodcells/LDeigen/'
    varY_file: '/project2/mstephens/yuxin/ukb-bloodcells/BloodCells.cov.rds'
    # number of dataset to evaluate
    n_dataset: 600
    # number of causal as a global variable, <0 is to use default
    # C: -1
    pve: 0.0005
    
    # Settings for small_compare
    C: 2
    nullz_file: ~/GitHub/mvarbvs/dsc/mnm_prototype/output/ukb_rss_small/nullz_cor.rds
