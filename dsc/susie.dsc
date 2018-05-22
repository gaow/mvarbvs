#!/usr/bin/env dsc
%include modules/zzz

DSC:
  define:
      get_data: liter_data 
      get_response: simple_lm
  run: 
    run_susie: get_data * summarize_ld * get_response * fit_susie * (plot_sse, plot_susie), 
        get_data * summarize_ld * get_response * get_sumstats * (fit_dap_z, fit_finemap)
    run_varbvs: get_data * summarize_ld * get_response * fit_varbvs * plot_sse
  output: benchmark
  exec_path: modules
  global:
    data_file: ../data/gtex-manifest.txt
    dap_g_data: ../data/dap-manifest.txt
