library(susieR)

susie_rss_analyze = function(z, R, L, estimate_residual_variance) {
  fit = tryCatch(susie_rss(z, R, L=L,
                           max_iter = 1000, estimate_residual_variance = estimate_residual_variance),
                 error = function(e) list(sets = NULL, pip=NULL))
  return(fit)
}

susie_rss_multiple = function(Z, R, L, estimate_residual_variance) {
  fitted = list()
  posterior = list()
  if (is.null(dim(Z))) Z = matrix(ncol=1, Z)
  for (r in 1:ncol(Z)) {
    fitted[[r]] = susie_rss_analyze(Z[,r], R, L, estimate_residual_variance)
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
R = readRDS(ld);
res = susie_rss_multiple(Z, R, L, estimate_residual_variance)

