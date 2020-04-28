simulate_y_base: regression_simulator.R + \
                R(res = simulate_main(X, Y, missing_Y, scale_Y, prior_file, eff_mode, n_signal, var_Y, residual_mode))
  @CONF: R_libs=susieR
  X: $X
  Y: $Y
  var_Y: $var_Y
  missing_Y: TRUE, FALSE
  scale_Y: TRUE
  prior_file: "${prior_file}"
  # set signal to <0 to use a default setting
  eff_mode: "artificial_mixture_50", "gtex_mixture"
  n_signal: ${C}
  residual_mode: "identity", "var_Y"
  # per-condition PVE
  pve: ${pve}
  $Y: res$Y
  $J: ncol(res$X)
  $R: res$n_traits
  $meta: list(true_coef=res$true_coef, X_csd=res$X_csd, X_mean=res$X_mean,
              residual_variance=res$residual_variance, original_Y=Y, prior=res$prior)

artificial_mixture(simulate_y_base):
    eff_mode: "artificial_mixture_50"
    missing_Y: FALSE

gtex_mixture(simulate_y_base):
    eff_mode: "gtex_mixture"
    missing_Y: FALSE