#!/usr/bin/env dsc

%include modules/data
%include modules/simulate_fixed_mix
%include modules/mnm
%include modules/mnm_rss
%include modules/score

DSC:
  define:
    simulate: artificial_mixture_ukb, ukb_bloodcells_mixture
    mnm: mnm_oracle, mnm_naive
    mnm_rss: mnm_oracle, mnm_rss_oracle, mnm_rss_naive, mnm_rss_ed, mnm_rss_identity, mnm_rss_shared
  run:
    default: data_ukb * simulate * mnm * susie_scores
    simulate_only: data_ukb * simulate # using command argument --n_dataset 20000 this should simulate 20000 data-sets and generate univariate summary stats, for learning about mixture prior using EB
  exec_path: code
  output: output/ukb_rss
  global:
    data_file: /project2/mstephens/yuxin/ukb-bloodcells/regions1000_5000.txt
    prior_file: /project2/mstephens/yuxin/ukb-bloodcells/ukb_prior_simulation.rds
    genotype_dir: '/project2/mstephens/yuxin/ukb-bloodcells/genotypes/'
    ld_dir: '/project2/mstephens/yuxin/ukb-bloodcells/LD/'
    ldeigen_dir: '/project2/mstephens/yuxin/ukb-bloodcells/LDeigen/'
    varY_file: '/project2/mstephens/yuxin/ukb-bloodcells/BloodCells.cov.rds'
    # number of dataset to evaluate
    n_dataset: 600
    # number of causal as a global variable, <0 is to use default
    C: -1
    pve: 0.002
