susie_rss: susie_rss.R
  @CONF: R_libs = susieR
  sumstats: $sumstats  
  ld: $ld
  L: 10
  estimate_residual_variance: TRUE, FALSE
  $fitted: res$fitted
  $posterior: res$posterior