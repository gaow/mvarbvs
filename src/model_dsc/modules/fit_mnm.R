source("utils.R")

update_mash_model <- function(X, Y, V, fitted_g) {
  ## result contains 'PosteriorMean' 'PosteriorSD' 'lfdr' 'NegativeProb' 'lfsr'
  reg <- mm_regression(X, Y)
  mash_data <- mashr::mash_set_data(reg[1,,], Shat = reg[2,,], V = as.matrix(V))
  return(mashr::mash(mash_data, g = fitted_g, fixg = TRUE, outputlevel=99))
}

## Initialize storage for results
p_alpha <- rep(1, ncol(data$X)) / ncol(data$X)
alpha <- matrix(0, ncol(data$X), maxL)
mu <- lapply(1:maxL, function(i) matrix(0, ncol(data$X), ncol(data$Y)))
Xr <- matrix(0, nrow(data$Y), ncol(data$Y))
for (l in 1:maxL) {
  ## remove the lth effect
  Xr <- Xr - data$X %*% (alpha[,l] * mu[[l]])
  ## update mash model
  mout <- update_mash_model(data$X, data$Y - Xr, data$V, model$fitted_g)
  ## update fitted values
  mu[[l]] <- mout$result$PosteriorMean
  l10bf <- mashr::get_log10bf(mout)
  alpha_post <- exp((l10bf - max(l10bf)) * log(10)) * p_alpha
  alpha[,l] <- alpha_post / sum(alpha_post)
  ## add back the updated lth effect
  Xr <- Xr + data$X %*% (alpha[,l] * mu[[l]])
}
model$alpha <- alpha
model$mu <- mu
model$Xr <- Xr
