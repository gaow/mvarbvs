simulate_y_base: regression_simulator.R + regression.R + \
                R(res = simulate_main(X, Y, missing_Y, scale_Y, prior_file, eff_mode, n_signal, var_Y, residual_mode, save_summary_stats))
  @CONF: R_libs=susieR
  # here, X has already been centered and scaled
  X: $X
  Y: $Y
  var_Y: $var_Y
  missing_Y: TRUE, FALSE
  scale_Y: TRUE
  save_summary_stats: TRUE
  prior_file: "${prior_file}"
  # set signal to <0 to use a default setting
  eff_mode: "artificial_mixture_50", "gtex_mixture"
  n_signal: ${C}
  residual_mode: "identity", "var_Y"
  # per-condition PVE
  pve: ${pve}
  $Y: res$Y
  $J: ncol(X)
  $R: res$n_traits
  $meta: list(true_coef=res$true_coef, residual_variance=res$residual_variance,
              original_Y=Y, Y_sd=res$Y_sd, prior=res$prior)
  $sumstats: res$sumstats

artificial_mixture(simulate_y_base):
    eff_mode: "artificial_mixture_50"
    residual_mode: "var_Y"
    missing_Y: FALSE

gtex_mixture(simulate_y_base):
    eff_mode: "gtex_mixture"
    residual_mode: "var_Y"
    missing_Y: FALSE

artificial_mixture_identity(simulate_y_base):
    eff_mode: "artificial_mixture_50"
    residual_mode: "identity"
    missing_Y: FALSE

gtex_mixture_identity(simulate_y_base):
    eff_mode: "gtex_mixture"
    residual_mode: "identity"
    missing_Y: FALSE