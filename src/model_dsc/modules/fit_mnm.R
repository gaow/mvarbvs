source("utils.R")

update_mash_model <- function(X, Y, V, fitted_g) {
  ## result contains 'PosteriorMean' 'PosteriorSD' 'lfdr' 'NegativeProb' 'lfsr'
  reg <- mm_regression(X, Y)
  mash_data <- mashr::mash_set_data(reg[1,,], Shat = reg[2,,], V = as.matrix(V))
  return(mashr::mash(mash_data, g = fitted_g, fixg = TRUE, outputlevel=99))
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
fitted <- list(p_alpha=p_alpha, alpha=alpha, mu=mu, Xr=Xr)
fitted_track <- list()
for (i in 1:maxI) {
  fitted <- update_mnmash_model(data$X, data$Y, data$V, model$fitted_g, fitted)
  fitted_track[[i]] <- fitted
}
