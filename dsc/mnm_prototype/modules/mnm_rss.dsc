mnm_rss_base: misc.R + mnm_rss.R
  @CONF: R_libs = (mmbr@stephenslab/mmbr, flashier@willwerscheid/flashier)
  sumstats: $sumstats  
  suffstats: $suffstats
  ld: $ld
  ldeigen: $ldeigen
  meta: $meta
  prior: 'oracle','identity', 'shared', 'naive', 'ED'
  L: 10
  resid_method: 'oracle', 'identity', 'nullz', 'varY'
  n_thread: 5
  nullz_file: '/project2/mstephens/yuxin/ukb-bloodcells/nullz_cor.rds'
  $result: result

mnm_rss_shared(mnm_rss_base):
  prior: 'shared'
  resid_method: 'oracle', 'identity', 'nullz', 'varY'

mnm_rss_oracle(mnm_rss_base):
  prior: 'oracle'
  resid_method: 'oracle', 'identity', 'nullz', 'varY'

mnm_rss_naive(mnm_rss_base):
  prior: 'naive'
  resid_method: 'oracle', 'identity', 'nullz', 'varY'

mnm_rss_ed(mnm_rss_base):
  prior: 'ED'
  resid_method: 'oracle', 'identity', 'nullz', 'varY'

mnm_rss_identity(mnm_rss_base):
  prior: 'identity'
  resid_method: 'oracle', 'identity', 'nullz', 'varY'

