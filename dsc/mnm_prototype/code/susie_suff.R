library(susieR)

susie_suff_analyze = function(XtX, Xty, yty, n , L) {
  fit = tryCatch(susie_suff_stat(XtX = XtX, Xty = Xty, yty = yty, n = n, L=L,
                           max_iter = 1000),
                 error = function(e) list(sets = NULL, pip=NULL))
  return(fit)
}

susie_suff_multiple = function(XtX, XtY, YtY, n, L) {
  fitted = list()
  posterior = list()
  if (is.null(dim(XtY))) XtY = matrix(ncol=1, XtY)
  for (r in 1:ncol(XtY)) {
    fitted[[r]] = susie_suff_analyze(XtX, XtY[,r], YtY[r,r], n, L)
    fitted[[r]]$cs_corr = susieR:::get_cs_correlation(fitted[[r]], Xcorr=cov2cor(XtX))
    if(is.null(fitted[[r]]$sets))
      posterior[[r]] = NULL
    else
      posterior[[r]] = summary(fitted[[r]])
  }
  return(list(fitted=fitted, posterior=posterior))
}

library(data.table);

XtX = readRDS(ld) * suffstats$N
res = susie_suff_multiple(XtX, as.matrix(suffstats$XtY), as.matrix(suffstats$YtY), suffstats$N, L)