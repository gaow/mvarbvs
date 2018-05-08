## M&M ash module core update
mnm_update_model <- function(X, Y, V, fitted_g, fitted, get_kl = FALSE) {
  ## "fitted" include p_alpha, alpha, mu and Xr
  maxL = ncol(fitted$alpha)
  for (l in 1:maxL) {
    ## remove the lth effect
    fitted$Xr <- fitted$Xr - X %*% (fitted$alpha[,l] * fitted$mu[[l]])
    ## update mash model
    reg <- mm_regression(X, Y - fitted$Xr)
    mash_data <- mashr::mash_set_data(reg[1,,], Shat = reg[2,,], V = V)
    mout <- mashr::mash(mash_data, g = fitted_g, fixg = TRUE, outputlevel=3)
    ## update fitted values
    fitted$mu[[l]] <- mout$result$PosteriorMean
    fitted$s[[l]] <- mout$result$PosteriorCov
    fitted$lfsr[[l]] <- mout$result$lfsr
    fitted$neg[[l]] <- mout$result$NegativeProb
    l10bf <- mashr::get_log10bf(mout)
    alpha_post <- exp((l10bf - max(l10bf)) * log(10)) * fitted$p_alpha
    fitted$alpha[,l] <- alpha_post / sum(alpha_post)
    ## add back the updated lth effect
    fitted$Xr <- fitted$Xr + X %*% (fitted$alpha[,l] * fitted$mu[[l]])
    if (get_kl) {
        # Justified by A.46 of FLASH paper
        # Here KL is denoted as (13.28) of BDA 3
        fitted$kl[l] <- -1 * mout$loglik + compute_sse_Eloglik(reg[1,,], reg[2,,], V,
                                                               fitted$Sigma, 
                                                               fitted$alpha[,l],
                                                               mout$result$PosteriorMean,
                                                               mout$result$PosteriorCov)
    }
  }
  return(fitted)
}

## Compute posterior mean and covariances
mnm_compute_posterior_matrices = function(fitted, J, R, L) {
    post_mean <- matrix(0, J, R)
    for (l in 1:L) {
      post_mean <- post_mean + fitted$mu[[l]] * fitted$alpha[,l]
    }
    post_cov <- array(0, dim=c(R, R, J))
    for (j in 1:J) {
      for (l in 1:L) {
        post_cov[,,j] <- post_cov[,,j] + (fitted$mu[[l]][j,] %*% t(fitted$mu[[l]][j,]) + fitted$s[[l]][,,j]) * fitted$alpha[j,l]
      }
      post_cov[,,j] <- post_cov[,,j] - post_mean[j,] %*% t(post_mean[j,])
    }
    return(list(PosteriorMean = post_mean, PosteriorCov = post_cov))
}

## Initialize storage for results
data$X <- as.matrix(data$X)
data$Y <- as.matrix(data$Y)
maxL <- min(maxL, ncol(data$X))
p_alpha <- rep(1, ncol(data$X)) / ncol(data$X)
alpha <- matrix(0, ncol(data$X), maxL)
mu <- lapply(1:maxL, function(i) matrix(0, ncol(data$X), ncol(data$Y)))
Xr <- matrix(0, nrow(data$Y), ncol(data$Y))
fitted <- list(p_alpha=p_alpha, alpha=alpha, mu=mu, s=list(), Xr=Xr, kl=vector(), lfsr=list(), neg=list(), Sigma=V)
fitted_track <- list()
Vcorr <- cov2cor(V)
## For ELBO
XtX <- t(data$X) %*% data$X
## Fit m&m model
for (i in 1:maxI) {
  fitted <- mnm_update_model(data$X, data$Y, Vcorr, model$fitted_g, fitted, get_elbo)
  if (get_elbo) {
      post_mat = mnm_compute_posterior_matrices(fitted, ncol(data$X), ncol(data$Y), maxL)
      fitted$Sigma = compute_mnm_residual_covariance(data$X, data$Y, XtX,
                                                     post_mat$PosteriorMean, 
                                                     post_mat$PosteriorCov)
      fitted$post_loglik = compute_mnm_Eloglik(data$X, data$Y, 
                                          XtX, fitted$Sigma,
                                          post_mat$PosteriorMean, 
                                          post_mat$PosteriorCov)
      fitted$elbo = fitted$post_loglik - sum(fitted$kl)
  }
  fitted_track[[i]] <- fitted
}

post_mat = mnm_compute_posterior_matrices(fitted, ncol(data$X), ncol(data$Y), maxL)

## Compute lfsr
lfsr <- do.call(rbind, lapply(1:maxL, function(l) colSums(fitted$alpha[,l] * fitted$lfsr[[l]])))
posterior <- list(PosteriorMean=post_mat$PosteriorMean,
                  PosteriorCov=post_mat$PosteriorCov,
                  alpha = fitted$alpha,
                  lfsr=lfsr,
                  n_in_CI=susieR:::n_in_CI(t(fitted$alpha)),
                  in_CI=susieR:::in_CI(t(fitted$alpha))
                  )
