library("susieR")
fitted <- list()
n <- nrow(data$X)
for (r in 1:ncol(data$Y)) {
  fitted[[r]] <- susie(data$X,data$Y[,r],sa2=sa^2, L=maxL,max_iter=maxI)
  pos_prob <- pnorm(0,mean=t(fitted[[r]]$mu),sd=sqrt(fitted[[r]]$mu2-fitted[[r]]$mu^2))
  neg_prob <- 1-pos_prob
  fitted[[r]]$lfsr <- 1-colSums(fitted[[r]]$alpha * t(pmax(pos_prob,neg_prob)))
}
post_mean <- do.call(cbind, lapply(1:length(fitted), function(i) susieR:::coef.susie(fitted[[i]])))
lfsr <- do.call(cbind, lapply(1:length(fitted), function(i) fitted[[i]]$lfsr))
posterior <- list(PosteriorMean=post_mean, lfsr=lfsr)
