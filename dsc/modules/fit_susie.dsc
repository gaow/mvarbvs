
fit_susie02: fit_susie.R
  # Prior variance of nonzero effects.
  @CONF: R_libs = susieR@stephenslab/susieR
  maxI: 50
  maxL: 5
  estimate_residual_variance: FALSE
  prior_var: 0.2
  auto: FALSE
  data: $data
  $posterior: posterior
  $fitted: fitted

fit_susie01(fit_susie02):
  prior_var: 0.1

fit_susie04(fit_susie02):
  prior_var: 0.4

fit_susie05(fit_susie02):
  estimate_residual_variance: TRUE

fit_susie_auto(fit_susie02):
  auto: TRUE
