#!/usr/bin/env dsc
%include modules/zzz

DSC:
  define:
      setup: (full_data, lite_data, liter_data) * summarize_ld
      get_response: base_sim, original_Y
  run: 
    debug_mnm: setup * get_response * get_sumstats * init_mnm * fit_mnm_debug * plot_sse
    run_mnm: setup * get_response * get_sumstats * init_mnm * fit_mnm * plot_sse
    run_susie: setup * get_response * fit_susie * plot_sse
    run_varbvs: setup * get_response * fit_varbvs * plot_sse
    run_susie_uni: dap_g_data * original_Y * fit_susie * plot_sse
    run_varbvs_uni: dap_g_data * original_Y * fit_varbvs * plot_sse
  output: benchmark
  exec_path: modules
  global:
    data_file: "~/Documents/GTExV8/Thyroid.Lung.FMO2.filled.rds"
    dap_g_data: "../data/sim.r.geno.rds"
