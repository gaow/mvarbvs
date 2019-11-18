mnm_base: mnm.R
  @CONF: R_libs = (mmbr@gaow/mmbr, flashier@willwerscheid/flashier)
  X: $X
  Y: $Y
  meta: $meta
  cfg: $configurations
  eff_mode: 'identity', 'low_het', 'mid_het', 'high_het', 'shared', 'singleton', 'mixture_1'
  alpha: 0, 1
  L: 2, 1
  missing_Y: TRUE, FALSE
  resid_method: 'oracle', 'diag', 'flash'
  $result: result

mnm_identity(mnm_base):
  eff_mode: 'identity'

mnm_low_het(mnm_base):
  eff_mode: 'low_het'

mnm_mid_het(mnm_base):
  eff_mode: 'mid_het'

mnm_high_het(mnm_base):
  eff_mode: 'high_het'

mnm_shared(mnm_base):
  eff_mode: 'shared'

mnm_singleton(mnm_base):
  eff_mode: 'singleton'

mnm_mixture01(mnm_base):
  eff_mode: 'mixture_1'