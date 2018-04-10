# Initialize model data: priors and init values
source("utils.R")

if (Sigma == 'empirical') {
  data$V = cor(data$Y)
} else {
  ## FIXME: add other methods to compute Sigma
  data$V = cor(data$Y)
}
reg = univariate_regression(data$X, data$Y)
mash_data = mashr::mash_set_data(as.matrix(reg$betahat), Shat = as.matrix(reg$sebetahat), V = as.matrix(data$V))
if (U == 'auto') {
  U = mashr::cov_canonical(mash_data)
} else {
  ## FIXME: add other methods to get U
  U = mashr::cov_canonical(mash_data)
}
if (p == 'auto') {
  data$fitted_g = mashr::mash(mash_data, Ulist=U, outputlevel = 1)$fitted_g
} else {
  ## FIXME: need to use pre-fitted pi on larger data from mash procedure
  data$fitted_g = list(pi=p, Ulist=U, grid=grid, usepointmass=F)
}
