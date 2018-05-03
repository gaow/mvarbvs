# Initialize model data: priors and init values

if (Sigma == 'empirical') {
  data$V = cor(data$Y)
} else {
  ## FIXME: add other methods to compute Sigma
  data$V = cor(data$Y)
}
mash_data = mashr::mash_set_data(reg[1,,], Shat = reg[2,,], V = as.matrix(data$V))
if (U == 'auto') {
  U = mashr::cov_canonical(mash_data)
} else {
  ## FIXME: add other methods to get U
  U = mashr::cov_canonical(mash_data)
}
model = list()
if (p == 'auto') {
  model$fitted_g = mashr::mash(mash_data, Ulist=U, outputlevel=1, usepointmass=TRUE)$fitted_g
} else {
  ## FIXME: need to use pre-fitted pi on larger data from mash procedure
  model$fitted_g = list(pi=p, Ulist=U, grid=grid, usepointmass=TRUE)
}
