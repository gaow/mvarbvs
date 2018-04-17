## Perform univariate regression for each column of Y on each column of X
univariate_regression = function(X, y, Z=NULL, return_residue=FALSE) {
  if (!is.null(Z)) {
    y = .lm.fit(Z, y)$residuals
  }
  calc_stderr = function(X, residuals) {
    # S = (X'X)^-1 \Sigma
    sqrt(diag(sum(residuals^2) / (nrow(X) - 2) * chol2inv(chol(t(X) %*% X))))
  }
  output = do.call(rbind,
                   lapply(1:ncol(X), function(i) {
                     g = .lm.fit(cbind(1, X[,i]), y)
                     return(c(coef(g)[2], calc_stderr(cbind(1, X[,i]), g$residuals)[2]))
                   })
                   )
  if (return_residue) {
    return(list(betahat = output[,1], sebetahat = output[,2],
                residuals = y))
  } else {
    return(list(betahat = output[,1], sebetahat = output[,2]))
  }
}

library(abind)
mm_regression = function(X, Y, Z=NULL) {
  reg = lapply(seq_len(ncol(Y)), function (i) simplify2array(univariate_regression(X, Y[,i])))
  reg = do.call(abind, c(reg, list(along=0)))
  # return array: out[1,,] is betahat, out[2,,] is shat
  return(aperm(reg, c(3,2,1)))
}

#find how many variables in the 95% CI
# x is a probability vector
n_in_CI_x = function(x){
  sum(cumsum(sort(x,decreasing = TRUE))<0.95)+1
}

# return binary vector indicating if each point is in CI
# x is a probability vector
in_CI_x = function(x){
  n = n_in_CI_x(x)
  o = order(x,decreasing=TRUE)
  result = rep(0,length(x))
  result[o[1:n]] = 1
  return(result)
}

# Return binary matrix indicating which variables are in CI of each
# of effect
in_CI = function(alpha){
  t(apply(alpha,1,in_CI_x))
}

n_in_CI = function(alpha){
  apply(alpha,1,n_in_CI_x)
}
