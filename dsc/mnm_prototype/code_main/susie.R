library(susieR)
fitted <- list()
posterior <- list()
for (r in 1:ncol(Y)) {
  fitted[[r]] <- susie(X,Y[,r],
                       L=L,
                       max_iter=1000,
                       estimate_residual_variance=TRUE,
                       estimate_prior_variance=TRUE)
  fitted[[r]]$cs_corr = susieR:::get_cs_correlation(fitted[[r]], X=X)
  posterior[[r]] <- summary(fitted[[r]])
}
