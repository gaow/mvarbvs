update_mash_model <- function(X, Y, V, fitted_g) {
  ## result contains 'PosteriorMean' 'PosteriorSD' 'lfdr' 'NegativeProb' 'lfsr'
  reg <- mm_regression(X, Y)
  mash_data <- mashr::mash_set_data(reg[1,,], Shat = reg[2,,], V = as.matrix(V))
  return(mashr::mash(mash_data, g = fitted_g, fixg = TRUE, outputlevel=3))
}

update_mnmash_model <- function(X, Y, V, fitted_g, fitted) {
  ## "fitted" include p_alpha, alpha, mu and Xr
  maxL = ncol(fitted$alpha)
  for (l in 1:maxL) {
    ## remove the lth effect
    fitted$Xr <- fitted$Xr - X %*% (fitted$alpha[,l] * fitted$mu[[l]])
    ## update mash model
    mout <- update_mash_model(X, Y - fitted$Xr, V, fitted_g)
    ## update fitted values
    fitted$mu[[l]] <- mout$result$PosteriorMean
    fitted$s[[l]] <- mout$result$PosteriorCov
    fitted$eb[[l]] <- mout$result$elbo_base
    fitted$lfsr[[l]] <- mout$result$lfsr
    fitted$neg[[l]] <- mout$result$NegativeProb
    l10bf <- mashr::get_log10bf(mout)
    alpha_post <- exp((l10bf - max(l10bf)) * log(10)) * fitted$p_alpha
    fitted$alpha[,l] <- alpha_post / sum(alpha_post)
    ## add back the updated lth effect
    fitted$Xr <- fitted$Xr + X %*% (fitted$alpha[,l] * fitted$mu[[l]])
  }
  return(fitted)
}

## Initialize storage for results
p_alpha <- rep(1, ncol(data$X)) / ncol(data$X)
alpha <- matrix(0, ncol(data$X), maxL)
mu <- lapply(1:maxL, function(i) matrix(0, ncol(data$X), ncol(data$Y)))
Xr <- matrix(0, nrow(data$Y), ncol(data$Y))
fitted <- list(p_alpha=p_alpha, alpha=alpha, mu=mu, s=list(), Xr=Xr, eb=list(), lfsr=list(), neg=list())
fitted_track <- list()

## Fit m&m model
for (i in 1:maxI) {
  fitted <- update_mnmash_model(data$X, data$Y, data$V, model$fitted_g, fitted)
  fitted_track[[i]] <- fitted
}

## Compute posterior mean and covariances
post_mean <- matrix(0, ncol(data$X), ncol(data$Y))
for (l in 1:maxL) {
  post_mean <- post_mean + fitted$mu[[l]] * fitted$alpha[,l]
}
post_cov <- array(0, dim=c(ncol(data$Y), ncol(data$Y), ncol(data$X)))
for (j in 1:ncol(data$X)) {
  for (l in 1:maxL) {
    post_cov[,,j] <- post_cov[,,j] + (fitted$mu[[l]][j,] %*% t(fitted$mu[[l]][j,]) + fitted$s[[l]][,,j]) * fitted$alpha[j,l]
  }
  post_cov[,,j] <- post_cov[,,j] - post_mean[j,] %*% t(post_mean[j,])
}

## Compute lfsr
lfsr <- do.call(rbind, lapply(1:maxL, function(l) colSums(fitted$alpha[,l] * fitted$lfsr[[l]])))
posterior <- list(PosteriorMean=post_mean,
                  PosteriorCov=post_cov,
                  alpha = fitted$alpha,
                  lfsr=lfsr,
                  n_in_CI=susieR:::n_in_CI(t(fitted$alpha)),
                  in_CI=susieR:::in_CI(t(fitted$alpha))
                  )
