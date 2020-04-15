mnm_base: misc.R + mnm.R
  @CONF: R_libs = (mmbr@gaow/mmbr, flashier@willwerscheid/flashier)
  X: $X
  Y: $Y
  meta: $meta
  eff_mode: 'identity', 'low_het', 'mid_het', 'high_het', 'shared', 'singleton', 'singleton_1', 'mixture_1'
  L: 1, 2, 10
  missing_Y: TRUE, FALSE
  resid_method: 'oracle', 'diag', 'flash'
  $result: result

mnm_high_het(mnm_base):
  eff_mode: 'high_het'
  missing_Y: FALSE
  resid_method: 'oracle', 'flash'

mnm_high_het_missing(mnm_high_het):
  missing_Y: TRUE
  resid_method: 'oracle', 'flash'

mnm_low_het(mnm_base):
  eff_mode: 'low_het'
  missing_Y: FALSE
  resid_method: 'oracle', 'flash'

mnm_mid_het(mnm_base):
  eff_mode: 'mid_het'
  missing_Y: FALSE
  resid_method: 'oracle', 'flash'

mnm_shared(mnm_base):
  eff_mode: 'shared'
  missing_Y: FALSE
  resid_method: 'oracle', 'flash'

mnm_singleton(mnm_base):
  eff_mode: 'singleton'
  missing_Y: FALSE
  resid_method: 'oracle', 'flash'

mnm_singleton_first(mnm_base):
  eff_mode: 'singleton_1'
  missing_Y: FALSE
  resid_method: 'oracle', 'flash'

mnm_mixture01(mnm_base):
  eff_mode: 'mixture_1'
  missing_Y: FALSE
  resid_method: 'oracle', 'flash'

mnm_oracle(mnm_base):
  eff_mode: 'oracle'
  missing_Y: FALSE
  resid_method: 'oracle', 'flash'

mnm_identity(mnm_base):
  eff_mode: 'identity'
  missing_Y: FALSE
  resid_method: 'oracle', 'flash'

mnm_naive(mnm_base):
  eff_mode: 'naive'
  missing_Y: FALSE
  resid_method: 'oracle', 'flash'