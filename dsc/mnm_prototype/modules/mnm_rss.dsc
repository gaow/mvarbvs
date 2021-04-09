mnm_rss_base: misc.R + mnm_rss.R
  @CONF: R_libs = (mvsusieR@stephenslab/mvsusieR, flashier@willwerscheid/flashier)
  sumstats: $sumstats  
  suffstats: $suffstats
  ld: $ld
  meta: $meta
  prior: 'oracle','identity', 'shared', 'naive', 'ED'
  L: 10
  maxiter: 1000
  resid_method: 'oracle', 'identity', 'nullz', 'corY'
  n_thread: 4
  nullz_file: "${nullz_file}" 
  $result: result
  
mnm_rss_oracle(mnm_rss_base):
  prior: 'oracle'
  resid_method: 'oracle', 'identity', 'nullz', 'corY'

mnm_rss_shared(mnm_rss_base):
  prior: 'shared'
  resid_method: 'oracle'

mnm_rss_naive(mnm_rss_base):
  prior: 'naive'
  resid_method: 'oracle'

mnm_rss_ed(mnm_rss_base):
  prior: 'ED'
  resid_method: 'oracle'

mnm_rss_identity(mnm_rss_base):
  prior: 'identity'
  resid_method: 'oracle'

mnm_rss_shared_corZ(mnm_rss_base):
  prior: 'shared'
  resid_method: 'nullz'

mnm_rss_naive_corZ(mnm_rss_base):
  prior: 'naive'
  resid_method: 'nullz'

mnm_rss_ed_corZ(mnm_rss_base):
  prior: 'ED'
  resid_method: 'nullz'

mnm_rss_identity_corZ(mnm_rss_base):
  prior: 'identity'
  resid_method: 'nullz'

mnm_rss_shared_corY(mnm_rss_base):
  prior: 'shared'
  resid_method: 'corY'

mnm_rss_naive_corY(mnm_rss_base):
  prior: 'naive'
  resid_method: 'corY'

mnm_rss_ed_corY(mnm_rss_base):
  prior: 'ED'
  resid_method: 'corY'

mnm_rss_identity_corY(mnm_rss_base):
  prior: 'identity'
  resid_method: 'corY'
  
