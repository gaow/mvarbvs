#!/usr/bin/env dsc

%include modules/data_ukb
%include modules/simulate_fixed_mix
%include modules/mnm
%include modules/mnm_rss
%include modules/susie_rss
%include modules/paintor
%include modules/flashfm
%include modules/cafeh
%include modules/mthess
%include modules/score

DSC:
  define:
    simulate: artificial_mixture_ukb, ukb_bloodcells_mixture
    simulate_small2: artificial_mixture_ukb_small2, artificial_mixture_ukb_small2_indep
    mnm: mnm_oracle, mnm_naive
    mnm_rss: mnm_rss_oracle, mnm_rss_shared,  mnm_rss_identity, mnm_rss_naive, mnm_rss_ed, mnm_rss_edscale, mnm_rss_shared_corZ, mnm_rss_identity_corZ, mnm_rss_naive_corZ, mnm_rss_ed_corZ, mnm_rss_edscale_corZ, mnm_rss_shared_corY, mnm_rss_identity_corY, mnm_rss_naive_corY, mnm_rss_ed_corY, mnm_rss_edscale_corY
    susie_uni: susie_rss
  run:
    default: data_ukb * simulate * (mnm_rss * mvsusie_scores, susie_uni * susie_scores)
    cafeh: data_ukb * simulate * cafeh_prepare * cafeh * cafeh_scores
    small2_compare: data_ukb * artificial_mixture_ukb_small2 * ((mnm_rss_naive, mnm_rss_naive_corY, mnm_rss_naive_corZ) * mvsusie_scores, FLASHFMwithFINEMAP * flashfm_scores, PAINTOR, cafeh_prepare * cafeh * cafeh_scores)
    small2_indep_compare: data_ukb * artificial_mixture_ukb_small2_indep * ((mnm_rss_naive, mnm_rss_naive_corY, mnm_rss_naive_corZ) * mvsusie_scores, FLASHFMwithFINEMAP * flashfm_scores, PAINTOR, cafeh_prepare * cafeh * cafeh_scores)
    small4_compare: data_ukb * artificial_mixture_ukb_small4 * ((mnm_rss_naive, mnm_rss_naive_corY, mnm_rss_naive_corZ) * mvsusie_scores, PAINTOR)
    full_data: data_ukb * artificial_mixture_ukb_small2_indep * mthess
    simulate_small2_only: data_ukb * artificial_mixture_ukb_small2
    simulate_small2_indep_only: data_ukb * artificial_mixture_ukb_small2_indep
    simulate_small3_only: data_ukb * artificial_mixture_ukb_small3
    simulate_small4_only: data_ukb * artificial_mixture_ukb_small4
    simulate_only: data_ukb * simulate 
  exec_path: code
  output: output/ukb_rss
  global:
    data_file: /project2/mstephens/yuxin/ukb-bloodcells/regions1000_5000.txt
    prior_file: /project2/mstephens/yuxin/ukb-bloodcells/analysis_20220619/ukb_prior_simulation_20220619.rds
    # nullz_file: /project2/mstephens/yuxin/ukb-bloodcells/analysis_20220619/nullz_cor_20220619.rds
    genotype_dir: '/project2/mstephens/yuxin/ukb-bloodcells/genotypes/'
    ld_dir: '/project2/mstephens/yuxin/ukb-bloodcells/LD/'
    # varY_file: '/project2/mstephens/yuxin/ukb-bloodcells/analysis_20220619/ukbbloodcells_prepare.Ycor.rds'
    # number of dataset to evaluate
    n_dataset: 600
    # number of causal as a global variable, <0 is to use default
    # C: -1
    pve: 0.0005
    
    # Settings for small2_compare
    # C: -1
    # varY_file: '/project2/mstephens/yuxin/ukb-bloodcells/analysis_20220619/artificial_small2.Ycor.rds'
    # nullz_file: /project2/mstephens/yuxin/ukb-bloodcells/analysis_20220619/artificial_mixture_2_20220619_nullzcor.rds
    # Settings for small2_indep_compare
    C: 2
    varY_file: '/project2/mstephens/yuxin/ukb-bloodcells/analysis_20220619/artificial_small2.Ycor.rds'
    nullz_file: /project2/mstephens/yuxin/ukb-bloodcells/analysis_20220619/artificial_mixture_2_indep_20220619_nullzcor.rds
    # C: 2
    # nullz_file: ~/GitHub/mvarbvs/dsc/mnm_prototype/output/ukb_rss_small/nullz_cor.rds
