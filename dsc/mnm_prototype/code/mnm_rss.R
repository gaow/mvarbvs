library(data.table)
Z = as.matrix(sumstats$bhat/sumstats$sbhat)
Z[is.na(Z)] = 0
if(prior == 'EDscale'){
  priorU = meta$prior[['ED']]
  priorU$xUlist = lapply(priorU$xUlist, function(x){
    if(max(diag(x))!=0){
      x = x/max(diag(x))
    }
    return(x)
  })
}else{
  priorU = meta$prior[[prior]]
}
if(resid_method == 'identity'){
  resid_Z = diag(ncol(Z))
}else if(resid_method == 'oracle'){
  resid_Z = meta$residual_variance
}else if(resid_method == 'nullz'){
  resid_Z = readRDS(nullz_file)[[meta$eff_mode]]
}else if(resid_method == 'corY'){
  resid_Z = cov2cor(as.matrix(suffstats$YtY) / (suffstats$N-1))
}

# LD info can be provided either as R data objects or R RDS files
if (is.character(ld)) {
  LD = readRDS(ld)
} else {
  LD = ld
}

m_init = mvsusieR::create_mixture_prior(mixture_prior = list(matrices=priorU$xUlist, weights=priorU$pi), 
                                 null_weight=priorU$null_weight, max_mixture_len=-1)
result = mvsusieR::mvsusie_rss(Z, LD, N=suffstats$N, L=L, prior_variance=m_init, residual_variance=resid_Z, 
                          compute_objective=T, estimate_prior_variance=T, estimate_prior_method='EM', 
                          precompute_covariances=T, n_thread=n_thread, max_iter=maxiter)
result$cs_corr = susieR:::get_cs_correlation(result, Xcorr=LD)

