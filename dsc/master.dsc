#!/usr/bin/env dsc
%include modules/zzz

DSC:
  define:
    get_data: full_data, lite_data, liter_data, two_effect
    get_response: base_sim, original_Y
    fit: (init_mnm * fit_mnm * plot_sse), 
        (fit_susie * plot_sse), 
        (fit_varbvs * plot_sse), 
        (fit_finemap * plot_finemap), 
        (fit_dap * plot_dap),
        (fit_caviar * plot_caviar)
  run:
    setup: liter_data * summarize_ld
    benchmark: full_data * summarize_ld * get_Y * get_sumstats * fit
    debug_mnm_1: lite_data * summarize_ld * get_Y * get_sumstats * init_mnm * fit_mnm_debug
    debug_mnm_2: liter_data * summarize_ld * get_Y * get_sumstats * init_mnm * fit_mnm_debug   
    caviar: lite_data * summarize_ld * get_Y * get_sumstats * (fit_caviar * plot_caviar)
  output: benchmark
  exec_path: modules
  global:
    data_file: ~/Documents/GTExV8/Thyroid.Lung.FMO2.filled.rds
