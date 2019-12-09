mnm_base: mnm.R
  @CONF: R_libs = (mmbr@gaow/mmbr, flashier@willwerscheid/flashier)
  X: $X
  Y: $Y
  meta: $meta
  cfg: $configurations
  eff_mode: 'identity', 'low_het', 'mid_het', 'high_het', 'shared', 'singleton', 'mixture_1'
  alpha: 0, 1
  L: 1, 2, 10
  missing_Y: TRUE, FALSE
  resid_method: 'oracle', 'diag', 'flash'
  $result: result

mnm_high_het(mnm_base):
  eff_mode: 'high_het'
  missing_Y: FALSE
  resid_method: 'oracle', 'diag', 'flash'

mnm_high_het(mnm_high_het):
  missing_Y: TRUE
  resid_method: 'diag', 'flash'