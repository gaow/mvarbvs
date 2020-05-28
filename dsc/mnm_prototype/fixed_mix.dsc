#!/usr/bin/env dsc

%include modules/data
%include modules/simulate_fixed_mix
%include modules/mnm
%include modules/mnm_rss
%include modules/score
%include modules/mthess
%include modules/atlasqtl

DSC:
  define:
    simulate: artificial_mixture, gtex_mixture
    simulate_missing: artificial_mixture_missing
    simulate_identity: artificial_mixture_identity, gtex_mixture_identity
    mnm: mnm_oracle, mnm_naive, mnm_ed, mnm_identity, mnm_shared
    mnm_missing: mnm_oracle, mnm_rss_oracle, mnm_rss_naive
  run:
    default: small_data * simulate * (mnm * susie_scores, atlasqtl) #, mthess)
    missingdata: tiny_data * simulate_missing * mnm_missing * susie_scores
    simulate_only: full_data * simulate_identity # using command argument --n_dataset 20000 this should simulate 20000 data-sets and generate univariate summary stats, for learning about mixture prior using EB
    mthess: small_data * artificial_mixture_small * (mnm_oracle, mnm_naive, mnm_identity, mnm_shared, mthess, atlasqtl)
  exec_path: code
  output: finemap_fixed_mixture
  global:
    data_file: ../data/gtex-v8-manifest.txt
    prior_file: ../data/prior_simulation.rds
    # number of dataset to evaluate
    n_dataset: 500
    # number of causal as a global variable, <0 is to use default
    C: -1
    pve: 0.15
