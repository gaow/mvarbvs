susie_rss: susie_rss.R
  @CONF: R_libs = susieR
  sumstats: $sumstats  
  suffstats: $suffstats
  ld: $ld
  L: 10
  $fitted: res$fitted
  $posterior: res$posterior
  
susie_rss_cs: susie_rss_cs.R
  ld: $ld
  res: $fitted
  $sets: sets