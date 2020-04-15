#!/usr/bin/env dsc

%include modules/data
%include modules/simulate_fixed_mix
%include modules/mnm
%include modules/score
%include modules/atlasqtl

DSC:
  define:
    simulate: artificial_mixture, gtex_mixture
    method: mnm_oracle, mnm_naive, mnm_identity, mnm_shared #, atlasqtl, mnm_ed
    score: susie_scores
  run:
    default: lite_data * simulate * method * score
  exec_path: code
  output: finemap_fixed_mixture
  global:
    data_file: ../data/gtex-v8-manifest.txt
    prior_file: ../data/prior_simulation.rds
    # number of dataset to evaluate
    n_dataset: 300
    # number of causal as a global variable, <0 is to use default
    C: -1
    pve: 0.1