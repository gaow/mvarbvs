source("utils.R")

update_mash_model <- function(X, Y, V, fitted_g) {
  reg <- univariate_regression(X, Y)
  mash_data <- mashr::mash_set_data(as.matrix(reg$betahat), Shat = as.matrix(reg$sebetahat), V = as.matrix(V))
  ## result contains 'PosteriorMean' 'PosteriorSD' 'lfdr' 'NegativeProb' 'lfsr'
  return(mashr::mash_compute_posterior_matrices(fitted_g, mash_data))
}

model <- update_mash_model(data$X, data$Y, data$V, model$fitted_g)
