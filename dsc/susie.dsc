#!/usr/bin/env dsc
%include modules/zzz

DSC:
  define:
      setup: liter_data * summarize_ld
      get_response: simple_lm
  run: 
    run_susie: setup * get_response * fit_susie * plot_sse
    run_varbvs: setup * get_response * fit_varbvs * plot_sse
  output: benchmark
  exec_path: modules
  global:
    data_file: ~/Documents/GTExV8/Thyroid.Lung.FMO2.filled.rds
    dap_g_data: ../data/sim.r.geno.rds
