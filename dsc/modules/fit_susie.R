fitted <- list()
for (r in 1:ncol(data$Y)) {
  fitted[[r]] <- susieR::susie(data$X,data$Y[,r],
                               L=maxL,
                               max_iter=maxI,
                               estimate_residual_variance=estimate_s2)
  fitted[[r]]$lfsr <- susieR:::lfsr_fromfit(fitted[[r]])
  fitted[[r]]$n_in_CI <- susieR:::n_in_CI(fitted[[r]])
  fitted[[r]]$in_CI <- susieR:::in_CI(fitted[[r]])
}

posterior <- list(PosteriorMean=do.call(cbind, lapply(1:length(fitted), function(i) susieR:::coef.susie(fitted[[i]]))),
                  lfsr=do.call(cbind, lapply(1:length(fitted), function(i) fitted[[i]]$lfsr)),
                  alpha=do.call(cbind, lapply(1:length(fitted), function(i) fitted[[i]]$alpha)),
                  n_in_CI=do.call(cbind, lapply(1:length(fitted), function(i) fitted[[i]]$n_in_CI)),
                  in_CI= do.call(cbind, lapply(1:length(fitted), function(i) fitted[[i]]$in_CI))
                  )
