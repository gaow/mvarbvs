#!/usr/bin/env dsc

%include module/setup
%include module/fit
%include module/evaluate

DSC:
  define:
    get_Y: original_Y
    init: init_mnm
    fit: fit_mnm, fit_susie, fit_varbvs, fit_finemap, fit_dap, fit_dap_mv
  run:
    first_pass: get_data * get_Y * get_sumstats * init * fit
    dap: get_data * get_Y * get_sumstats * init * fit_dap
  output: mnm_model
  exec_path: module
  global:
    data_file: ~/Documents/GTExV8/Thyroid.Lung.FMO2.filled.rds