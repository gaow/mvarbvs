fit_susie: fit_susie.R
  # Prior variance of nonzero effects.
  @CONF: R_libs = susieR@stephenslab/susieR
  maxI: 200
  maxL: 5
  estimate_residual_variance: FALSE, TRUE
  prior_var: 0.05, 0.1, 0.2, 0.4
  auto: FALSE
  data: $data
  $posterior: posterior
  $fitted: fitted
