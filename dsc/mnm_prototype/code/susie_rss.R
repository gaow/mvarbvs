library(susieR)

susie_rss_analyze = function(z, R, n, L) {
  fit = tryCatch(susie_rss(z, R, n=n, L=L,
                           max_iter = 1000,
                           estimate_residual_variance = TRUE),
                 error = function(e) list(sets = NULL, pip=NULL))
  return(fit)
}

susie_rss_multiple = function(Z, R, n, L) {
  fitted = list()
  posterior = list()
  if (is.null(dim(Z))) Z = matrix(ncol=1, Z)
  for (r in 1:ncol(Z)) {
    fitted[[r]] = susie_rss_analyze(Z[,r], R, n, L)
    fitted[[r]]$cs_corr = susieR:::get_cs_correlation(fitted[[r]], Xcorr=R)
    if(is.null(fitted[[r]]$sets))
      posterior[[r]] = NULL
    else
      posterior[[r]] = summary(fitted[[r]])
  }
  return(list(fitted=fitted, posterior=posterior))
}

library(data.table);
Z = sumstats$bhat/sumstats$sbhat;
if (is.character(ld)) {
  R = readRDS(ld)
} else {
  R = ld
}
n = suffstats$N
res = susie_rss_multiple(Z, R, n, L)

