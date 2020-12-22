mnm_suff_base: misc.R + mnm_suff.R
  @CONF: R_libs = (mmbr@stephenslab/mmbr, flashier@willwerscheid/flashier)
  suffstats: $suffstats  
  ld: $ld
  meta: $meta
  eff_mode: 'oracle','identity', 'shared', 'naive', 'ED'
  L: 10
  resid_method: 'oracle', 'covY', 'diag'
  n_thread: 5
  $result: result
  $prior: prior

mnm_suff_shared(mnm_suff_base):
  eff_mode: 'shared'
  resid_method: 'oracle', 'covY', 'diag'

mnm_suff_oracle(mnm_suff_base):
  eff_mode: 'oracle'
  resid_method: 'oracle', 'covY', 'diag'

mnm_suff_naive(mnm_suff_base):
  eff_mode: 'naive'
  resid_method: 'oracle', 'covY', 'diag'

mnm_suff_ed(mnm_suff_base):
  eff_mode: 'ED'
  resid_method: 'oracle', 'covY', 'diag'

mnm_suff_identity(mnm_suff_base):
  eff_mode: 'identity'
  resid_method: 'oracle', 'covY', 'diag'

