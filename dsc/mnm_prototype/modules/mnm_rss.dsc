mnm_rss_base: misc.R + mnm_rss.R
  @CONF: R_libs = (mmbr@stephenslab/mmbr, flashier@willwerscheid/flashier)
  sumstats: $sumstats  
  ld: $ld
  ldeigen: $ldeigen
  meta: $meta
  Y: $Y
  eff_mode: 'oracle','identity', 'shared', 'naive', 'ED'
  L: 10
  resid_method: 'oracle', 'identity', 'nullz', 'varY'
  n_thread: 5
  $result: result

mnm_rss_shared(mnm_rss_base):
  eff_mode: 'shared'
  resid_method: 'oracle', 'identity', 'nullz', 'varY'

mnm_rss_oracle(mnm_rss_base):
  eff_mode: 'oracle'
  resid_method: 'oracle', 'identity', 'nullz', 'varY'

mnm_rss_naive(mnm_rss_base):
  eff_mode: 'naive'
  resid_method: 'oracle', 'identity', 'nullz', 'varY'

mnm_rss_ed(mnm_rss_base):
  eff_mode: 'ED'
  resid_method: 'oracle', 'identity', 'nullz', 'varY'

mnm_rss_identity(mnm_rss_base):
  eff_mode: 'identity'
  resid_method: 'oracle', 'identity', 'nullz', 'varY'

