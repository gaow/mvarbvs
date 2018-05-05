#!/usr/bin/env dsc

%include modules/setup
%include modules/fit
%include modules/evaluate

DSC:
  define:
    get_Y: original_Y
    fit: (init_mnm * fit_mnm), fit_susie, fit_varbvs, fit_finemap, fit_dap #, fit_dap_mv
  run:
    first_pass: get_data * get_Y * get_sumstats * fit
  output: benchmark
  exec_path: modules
  global:
    data_file: ~/Documents/GTExV8/Thyroid.Lung.FMO2.filled.rds
