#!/usr/bin/env dsc
%include modules/zzz

DSC:
  define:
      setup: liter_data * summarize_ld
      get_response: simple_lm
      fit_susie: fit_susie01, fit_susie02, fit_susie04, fit_susie05 #, fit_susie_auto
  run: 
    run_susie: setup * get_response * fit_susie * (plot_sse, plot_susie)
    run_varbvs: setup * get_response * fit_varbvs * plot_sse
    run_finemap: setup * get_response * get_sumstats * fit_finemap * plot_finemap
  output: benchmark
  exec_path: modules
  global:
    data_file: "~/Documents/GTExV8/Thyroid.Lung.FMO2.filled.rds"
    dap_g_data: "../data/sim.r.geno.rds"
