library(data.table)
Z = as.matrix(sumstats$bhat/sumstats$sbhat)
Z[is.na(Z)] = 0
prior = meta$prior[[eff_mode]]
if (resid_method == 'all') {
    resid_Z <- cor(Z)
} else if (resid_method == 'null') {
    max_absz = apply(abs(Z),1, max)
    nullish = which(max_absz < 2)
    if(length(nullish)<=ncol(Z)){
      warning("not enough null data to estimate null correlation.")
      resid_Z <- diag(ncol(Z))
    }else{
      nullish_z = Z[nullish,]
      resid_Z <- cor(nullish_z)
      qrZ = qr(resid_Z)
      if (qrZ$rank != ncol(Z)) {
         warning("not enough null independent data to estimate null correlation.") 
         resid_Z =  diag(ncol(Z))
      }
    }
} else if (resid_method == 'identity'){
    resid_Z = diag(ncol(Z))
} else if (resid_method == 'flash'){
    resid_Z <- cov2cor(compute_cov_flash(Y))
} else{
    resid_Z <- meta$residual_variance
}
if(any(is.na(Y))){
    n_per_condition = crossprod(!is.na(Y))
    ### adjusted residual variance
    resid_Z = resid_Z * cov2cor(n_per_condition)
}
if(prior_scale == 'oracle'){
  xUlist = lapply(prior$xUlist, function(U) U * nrow(Y))
}else if(prior_scale == 'simulated'){
  xUlist = prior$xUlist
}
m_init = mmbr::create_mash_prior(mixture_prior = list(matrices=xUlist, weights=prior$pi), null_weight=prior$null_weight, max_mixture_len=-1)
result = mmbr::msusie_rss(Z, ld, L=L, prior_variance=m_init, residual_variance=resid_Z, compute_objective=TRUE, estimate_residual_variance=F, estimate_prior_variance=T, estimate_prior_method='EM', precompute_covariances=T, n_thread=n_thread, max_iter=1000)
#result$pip_conditions = mmbr::mmbr_get_pip_per_condition(result, m_init)
