#!/usr/bin/env dsc

%include modules/setup
%include modules/fit
%include modules/evaluate

DSC:
  define:
    get_data: full_data, lite_data, two_effect
    get_Y: original_Y
    fit: (init_mnm * fit_mnm), fit_susie, fit_varbvs, 
        (fit_finemap * plot_finemap), 
        (fit_dap * plot_dap)
  run:
    benchmark: full_data * get_Y * get_sumstats * fit
    debug_mnm_1: lite_data * get_Y * get_sumstats * init_mnm * fit_mnm_debug
    debug_mnm_2: liter_data * get_Y * get_sumstats * init_mnm * fit_mnm_debug   
    caviar: lite_data * get_Y * get_sumstats * (fit_caviar * plot_caviar)
  output: benchmark
  exec_path: modules
  global:
    data_file: ~/Documents/GTExV8/Thyroid.Lung.FMO2.filled.rds
