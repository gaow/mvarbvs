#!/usr/bin/env dsc
%include modules/zzz

DSC:
  define:
    setup: (full_data, lite_data, liter_data) * summarize_ld
    get_response: base_sim, original_Y
  run:
    run_dap_z_uni: dap_g_data * original_Y * get_sumstats * fit_dap_z * plot_dap
  output: benchmark
  exec_path: modules
  global:
    data_file: "~/Documents/GTExV8/Thyroid.Lung.FMO2.filled.rds"
    dap_g_data: "../data/sim.r.geno.rds"
