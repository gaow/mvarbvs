#!/usr/bin/env dsc
%include modules/zzz

DSC:
  run: full_data * summarize_ld * original_Y * (fit_susie * plot_susie, fit_dap * plot_dap)
  exec_path: modules
  global:
    data_file: ../data/gtex-manifest.txt
