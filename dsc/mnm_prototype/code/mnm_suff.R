library(data.table)
prior = meta$prior[[eff_mode]]

if(eff_mode %in% c('oracle', 'naive', 'identity', 'shared')){
  prior$xUlist = lapply(1:length(prior$xUlist), function(i) prior$xUlist[[i]]*suffstats$N)
}

suffstats$XtY = as.matrix(suffstats$XtY)
suffstats$YtY = as.matrix(suffstats$YtY)

if (resid_method == 'covY') {
  resid_Y <- suffstats$YtY / (suffstats$N-1)
} else if (resid_method == 'diag') {
  resid_Y <- diag(diag(suffstats$YtY / (suffstats$N-1)))
} else {
  resid_Y <- meta$residual_variance
}

XtX = readRDS(ld) * (suffstats$N - 1)

m_init = mmbr::create_mash_prior(mixture_prior = list(matrices=prior$xUlist, weights=prior$pi), 
                                 null_weight=prior$null_weight, max_mixture_len=-1)
result = mmbr::msusie_suff_stat(XtX, suffstats$XtY, suffstats$YtY, suffstats$N, L=L, 
                                prior_variance=m_init, residual_variance=resid_Y, 
                                compute_objective=TRUE, estimate_residual_variance=F, 
                                estimate_prior_variance=T, estimate_prior_method='EM', 
                                precompute_covariances=T, n_thread=n_thread, max_iter=1000)
result$cs_corr = susieR:::get_cs_correlation(result, Xcorr=cov2cor(XtX))