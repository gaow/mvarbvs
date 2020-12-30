library(data.table)
Z = as.matrix(sumstats$bhat/sumstats$sbhat)
Z[is.na(Z)] = 0
prior = meta$prior[[prior]]
if(resid_method == 'identity'){
  resid_Z = diag(ncol(Z))
}else if(resid_method == 'oracle'){
  resid_Z = meta$residual_variance
}else if(resid_method == 'nullz'){
  resid_Z = readRDS(nullz_file)[[meta$eff_mode]]
}else if(resid_method == 'corY'){
  resid_Z = cov2cor(suffstats$YtY / (suffstats$N-1))
}

LD = readRDS(ld)
ldeigen = readRDS(ldeigen)

m_init = mmbr::create_mash_prior(mixture_prior = list(matrices=prior$xUlist, weights=prior$pi), 
                                 null_weight=prior$null_weight, max_mixture_len=-1)
result = mmbr::msusie_rss(Z, LD, eigenR = ldeigen, L=L, prior_variance=m_init, residual_variance=resid_Z, 
                          compute_objective=TRUE, estimate_residual_variance=F, 
                          estimate_prior_variance=T, estimate_prior_method='EM', 
                          precompute_covariances=T, n_thread=n_thread, max_iter=1000)
result$pip_conditions = mmbr:::mmbr_get_pip_per_condition(result, m_init)
