#!/usr/bin/env dsc

%include modules/data

DSC:
  run:
    default: small_data * simulate
  exec_path: code
  global:
    data_file: ../data/gtex-v8-manifest.txt
    prior_file: ../data/prior_simulation.rds
    # number of dataset to evaluate
    n_dataset: 500
    # number of causal as a global variable, <0 is to use default
    C: -1
    # This is per eQTL heritability, with is_pve_total set to FALSE in regression_simulator.R
    # see Figure 1 of this paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4028055/
    pve: 0.15

simulate: misc.R + regression_simulator.R + \
          R(res = mr_mash_sim(X, Y, missing_Y, scale_Y, prior_file, eff_mode, n_signal, is_pve_total, var_Y, residual_mode)) 
  # here, X has already been centered and scaled
  X: $X
  Y: $Y
  var_Y: $var_Y
  missing_Y: FALSE
  scale_Y: TRUE
  prior_file: "${prior_file}"
  is_pve_total: FALSE
  # set signal to <0 to use a default setting
  eff_mode: "artificial_mixture_50", "gtex_mixture"
  n_signal: ${C}
  residual_mode: "identity", "var_Y"
  # per-condition PVE
  pve: ${pve}
  $Y: res$Y 
  $B: res$B
  $V: res$V
  $causal_variables: res$causal_variables
  $causal_responses: res$causal_responses
  $Sigma: res$Sigma