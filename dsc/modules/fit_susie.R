fitted <- list()
for (r in 1:ncol(data$Y)) {
  if (auto)
      fitted[[r]] <- susieR::susie_auto(data$X,data$Y[,r])
  else
      fitted[[r]] <- susieR::susie(data$X,data$Y[,r],
                               L=maxL,
                               max_iter=maxI,
                               estimate_residual_variance = estimate_residual_variance, 
                               prior_variance=prior_var, 
                               intercept=FALSE)
  fitted[[r]]$lfsr <- susieR:::lfsr_fromfit(fitted[[r]])
  fitted[[r]]$n_in_CI <- susieR:::n_in_CI(fitted[[r]])
  fitted[[r]]$in_CI <- susieR:::in_CI(fitted[[r]])
}

posterior <- list(PosteriorMean=do.call(cbind, lapply(1:length(fitted), function(i) susieR:::coef.susie(fitted[[i]]))),
                  lfsr=do.call(cbind, lapply(1:length(fitted), function(i) fitted[[i]]$lfsr)),
                  alpha=do.call(cbind, lapply(1:length(fitted), function(i) fitted[[i]]$alpha)),
                  in_CI=lapply(1:length(fitted), function(i) fitted[[i]]$in_CI),
                  n_in_CI=do.call(cbind, lapply(1:length(fitted), function(i) fitted[[i]]$n_in_CI))
                  )
