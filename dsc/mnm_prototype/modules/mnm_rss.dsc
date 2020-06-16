mnm_rss_base: mnm_rss.R
  @CONF: R_libs = (mmbr@stephenslab/mmbr, flashier@willwerscheid/flashier)
  sumstats: $sumstats  
  ld: $ld
  meta: $meta
  Y: $Y
  eff_mode: 'identity', 'low_het', 'mid_het', 'high_het', 'shared', 'singleton', 'singleton_1', 'mixture_1'
  L: 10
  resid_method: 'oracle'
  prior_scale: 'oracle'
  n_thread: 5
  $result: result
  $prior: prior

mnm_rss_high_het(mnm_rss_base):
  eff_mode: 'high_het'

mnm_rss_low_het(mnm_rss_base):
  eff_mode: 'low_het'

mnm_rss_mid_het(mnm_rss_base):
  eff_mode: 'mid_het'

mnm_rss_shared(mnm_rss_base):
  eff_mode: 'shared'
  prior_scale: 'oracle','simulated'
  resid_method: 'oracle'

mnm_rss_singleton(mnm_rss_base):
  eff_mode: 'singleton'

mnm_rss_singleton_first(mnm_rss_base):
  eff_mode: 'singleton_1'

mnm_rss_mixture01(mnm_rss_base):
  eff_mode: 'mixture_1'

mnm_rss_oracle(mnm_rss_base):
  eff_mode: 'oracle'
  prior_scale: 'oracle','simulated'
  resid_method: 'oracle'

mnm_rss_naive(mnm_rss_base):
  eff_mode: 'naive'
  prior_scale: 'oracle','simulated'
  resid_method: 'oracle'

mnm_rss_ed(mnm_rss_base):
  eff_mode: 'ED'
  prior_scale: 'oracle','simulated'
  resid_method: 'oracle'

mnm_rss_identity(mnm_rss_base):
  eff_mode: 'identity'
  prior_scale: 'oracle','simulated'
  resid_method: 'oracle'

