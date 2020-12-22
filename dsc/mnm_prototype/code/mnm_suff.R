library(data.table)
prior = meta$prior[[eff_mode]]

if (resid_method == 'covY') {
  resid_Y <- suffstats$YtY / (suffstats$N-1)
} else if (resid_method == 'diag') {
  resid_Y <- diag(diag(suffstats$YtY / (suffstats$N-1)))
} else {
  resid_Y <- meta$residual_variance
}

XtX = as.matrix(fread(ld)) * suffstats$N

m_init = mmbr::create_mash_prior(mixture_prior = list(matrices=xUlist, weights=prior$pi), 
                                 null_weight=prior$null_weight, max_mixture_len=-1)
result = mmbr::msusie_suff_stat(XtX, suffstats$XtY, suffstats$YtY, suffstats$N, L=L, 
                                prior_variance=m_init, residual_variance=resid_Z, 
                                compute_objective=TRUE, estimate_residual_variance=F, 
                                estimate_prior_variance=T, estimate_prior_method='EM', 
                                precompute_covariances=T, n_thread=n_thread, max_iter=1000)
result$pip_conditions = mmbr:::mmbr_get_pip_per_condition(result, m_init)
