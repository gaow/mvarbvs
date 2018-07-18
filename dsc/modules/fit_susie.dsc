fit_susie: fit_susie.R
  # Prior variance of nonzero effects.
  @CONF: R_libs = susieR@stephenslab/susieR
  maxI: 200
  maxL: 5
  estimate_residual_variance: FALSE, TRUE
  prior_var: 0, 0.05, 0.1, 0.2, 0.4
  data: $data
  $posterior: posterior
  $fitted: fitted

fit_susie_auto: fit_susie.R
  @CONF: R_libs = susieR@stephenslab/susieR
  data: $data
  prior_var: "auto"
  $posterior: posterior
  $fitted: fitted

fit_susie01(fit_susie):
  maxL: 1

fit_susie10(fit_susie):
  maxL: 10
