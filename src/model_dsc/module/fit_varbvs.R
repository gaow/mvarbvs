fitted <- list()
for (r in 1:ncol(Y)) {
  sigma <- var(Y[,r])
  fitted[[r]] <- varbvs::varbvsnorm(X,Y[,r],sigma,sa,logodds,alpha0,mu0,update.order = 1:p,
                                    update.sigma = FALSE,update.sa = FALSE,tol = 1e-6,
                                    verbose = FALSE, maxiter=maxI)
}

post_mean <- do.call(cbind, lapply(1:length(fitted), function(i) fitted[[i]]$alpha * fitted[[i]]$mu))
lfdr <- do.call(cbind, lapply(1:length(fitted), function(i) 1 - fitted[[i]]$alpha))
posterior <- list(PosteriorMean=post_mean, lfdr=lfdr)
