#!/usr/bin/env dsc
%include modules/zzz

DSC:
  run: 
    run_susie: liter_data * summarize_ld * simple_lm * fit_susie * (plot_sse, plot_susie)
    run_varbvs: liter_data * summarize_ld * simple_lm * fit_varbvs * plot_sse
    run_comparison: liter_data * summarize_ld * lm_less * get_sumstats * (fit_susie * plot_susie, fit_dap * plot_dap, fit_caviar * plot_caviar, fit_finemap * plot_finemap)
    hard_case: full_data * lm_less03 * get_sumstats * (fit_susie10 * plot_susie, fit_dap * plot_dap)
  exec_path: modules
  global:
    data_file: ../data/gtex-manifest.txt
    dap_g_data: ../data/dap-manifest.txt
