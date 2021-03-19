#!/usr/bin/env dsc

# DSC benchmark using fixed mixture patterns, for multivariate fine-mapping

%include modules/data
%include modules/simulate_fixed_mix
%include modules/mnm
%include modules/mnm_rss
%include modules/susie
%include modules/score
%include modules/mthess
%include modules/atlasqtl

DSC:
  define:
    simulate: artificial_mixture, gtex_mixture
    simulate_missing: artificial_mixture_missing, gtex_mixture_missing
    simulate_identity: artificial_mixture_identity, gtex_mixture_identity
    mnm: mnm_oracle, mnm_naive, mnm_ed, mnm_ed_max10, mnm_identity, mnm_shared
    mnm_rss: mnm_rss_oracle, mnm_rss_naive_corZ, mnm_rss_ed_corZ, mnm_rss_identity_corZ, mnm_rss_shared_corZ
  run:
    default: medium_data * simulate * (mnm * mvsusie_scores, susie * susie_scores, mnm_rss_ed_corZ * mvsusie_scores, atlasqtl) #, mthess)
    missing_data: medium_data * simulate_missing * (mnm, mnm_rss) * mvsusie_scores
    simulate_only: full_data * simulate # using command argument --n_dataset 20000 this should simulate 20000 data-sets and generate univariate summary stats, for learning about mixture prior using EB
    mthess: medium_data * artificial_mixture_small * (mnm_oracle, mnm_naive, mnm_identity, mnm_shared, mthess, atlasqtl)
  exec_path: code
  output: finemap_fixed_mixture
  global:
    data_file: ../data/gtex-v8-manifest.txt
    prior_file: ../data/prior_simulation.rds
    nullz_file: ../data/nullz_cor_simulation.rds
    # number of dataset to evaluate
    n_dataset: 500
    # number of causal as a global variable, <0 is to use some of the default options available, see get_n_signal() in regression_simulator.R
    C: -2
    # per eQTL heritability is 0.05 to 0.07 according to https://www.nature.com/articles/ng.3506
    # Here I'm using a total PVE 0.15
    pve: 0.15