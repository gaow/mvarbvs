#!/usr/bin/env dsc

%include modules/data_ukb
%include modules/simulate_fixed_mix
%include modules/mnm
%include modules/mnm_rss
%include modules/susie_rss
%include modules/paintor
%include modules/mthess
%include modules/score

DSC:
  define:
    simulate: artificial_mixture_ukb, ukb_bloodcells_mixture
    simulate_small: artificial_mixture_ukb_small4, artificial_mixture_ukb_small2, artificial_mixture_ukb_small2_indep
    mnm: mnm_oracle, mnm_naive
    mnm_rss: mnm_rss_oracle, mnm_rss_shared,  mnm_rss_identity, mnm_rss_naive, mnm_rss_ed, mnm_rss_edscale, mnm_rss_shared_corZ, mnm_rss_identity_corZ, mnm_rss_naive_corZ, mnm_rss_ed_corZ, mnm_rss_edscale_corZ, mnm_rss_shared_corY, mnm_rss_identity_corY, mnm_rss_naive_corY, mnm_rss_ed_corY, mnm_rss_edscale_corY
    susie_uni: susie_rss
  run:
    default: data_ukb * simulate * (mnm_rss * mvsusie_scores, susie_uni * susie_scores)
    small_compare: data_ukb * simulate_small * ((mnm_rss_naive, mnm_rss_naive_corY, mnm_rss_naive_corZ) * mvsusie_scores, PAINTOR)
    full_data: data_ukb * artificial_mixture_ukb_small2_indep * mthess
    simulate_small_only: data_ukb * simulate_small
    simulate_only: data_ukb * simulate 
  exec_path: code
  output: output/ukb_rss
  global:
    data_file: /project2/mstephens/yuxin/ukb-bloodcells/regions1000_5000.txt
    prior_file: /project2/mstephens/yuxin/ukb-bloodcells/analysis_20220619/ukb_prior_simulation_20220619.rds
    nullz_file: /project2/mstephens/yuxin/ukb-bloodcells/nullz_cor_20220429.rds
    genotype_dir: '/project2/mstephens/yuxin/ukb-bloodcells/genotypes/'
    ld_dir: '/project2/mstephens/yuxin/ukb-bloodcells/LD/'
    varY_file: '/project2/mstephens/yuxin/ukb-bloodcells/analysis_20220619/ukbbloodcells_prepare.Ycor.rds'
    # number of dataset to evaluate
    n_dataset: 600
    # number of causal as a global variable, <0 is to use default
    C: -1
    pve: 0.0005
    
    # Settings for small_compare
    # C: 2
    # nullz_file: ~/GitHub/mvarbvs/dsc/mnm_prototype/output/ukb_rss_small/nullz_cor.rds
