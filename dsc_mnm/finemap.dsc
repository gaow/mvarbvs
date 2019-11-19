#!/usr/bin/env dsc

%include modules/data
%include modules/simulate
%include modules/mnm
%include modules/score

DSC:
  define:
    simulate: high_het
    mnm: mnm_high_het
    score: susie_scores
  run:
    default: oracle_generator * full_data * simulate * mnm * score
  exec_path: code
  output: finemap_output
  global:
    data_file: ../data/gtex-v8-manifest.txt
    # number of dataset to evaluate
    n_dataset: 500
    # number of conditions as a global variable
    R: 5, 45
    # number of causal as a global variable, <0 is to use default
    C: 1
    #grid: (0.1, 0.25, 0.4, 0.6)
    grid: 0.3
    pve: 0.1
    mixture_1: dict(identity=0.1,equal_effects=0.2,singleton=0.2,simple_het_1=0.1,simple_het_2=0.1,simple_het_3=0.1,null=0)
