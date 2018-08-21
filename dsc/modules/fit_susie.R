fitted <- list()
for (r in 1:ncol(data$Y)) {
  if ('Z' %in% names(data)) {
      data$Y[,r] = residuals(lm(data$Y[,r]~Z, na.action=na.exclude))
  }
  if (prior_var == 'auto') {
      fitted[[r]] <- susieR::susie_auto(data$X,data$Y[,r],L_max=100,tol=1e-3)
  } else if (prior_var == 0) {
      fitted[[r]] <- susieR::susie(data$X,data$Y[,r],
                               L=maxL,
                               max_iter=maxI,
                               estimate_residual_variance=estimate_residual_variance,
                               estimate_prior_variance=TRUE,
                               tol=1e-3)
  } else {
      fitted[[r]] <- susieR::susie(data$X,data$Y[,r],
                               L=maxL,
                               max_iter=maxI,
                               estimate_residual_variance=estimate_residual_variance,
                               scaled_prior_variance=prior_var,
                               tol=1e-3)
  }
  fitted[[r]]$lfsr <- susieR::susie_get_lfsr(fitted[[r]])
  fitted[[r]]$n_in_CI <- susieR:::n_in_CS(fitted[[r]])
  fitted[[r]]$in_CI <- susieR:::in_CS(fitted[[r]])
}

posterior <- list(PosteriorMean=do.call(cbind, lapply(1:length(fitted), function(i) susieR:::coef.susie(fitted[[i]]))),
                  lfsr=do.call(cbind, lapply(1:length(fitted), function(i) fitted[[i]]$lfsr)),
                  alpha=lapply(1:length(fitted), function(i) fitted[[i]]$alpha),
                  in_CI=lapply(1:length(fitted), function(i) fitted[[i]]$in_CI),
                  n_in_CI=do.call(cbind, lapply(1:length(fitted), function(i) fitted[[i]]$n_in_CI)),
                  niter = sapply(1:length(fitted), function(i) fitted[[i]]$niter)
                  )
