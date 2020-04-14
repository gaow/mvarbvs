simulate_y: regression_simulator.R + \
                R(res = simulate_main(X, meta, eff_mode, n_traits, n_signal, residual_mode))
  @CONF: R_libs=susieR
  X: $X
  Y: $Y
  n_traits: ${R}
  # set signal to <0 to use a default setting
  n_signal: ${C}
  eff_mode: "mixture_01"
  residual_mode: "identity"
  # per-condition PVE
  pve: ${pve}
  $Y: res$Y
  $R: ncol(res$Y)
  $J: ncol(res$X)
  $pve_out: conf['pve']
  $meta: dict(true_coef=res['true_coef'], X_csd = res['X_csd'], X_mean = res['X_mean'],
              residual_variance=res['residual_variance'], original_Y=Y)