## Perform univariate regression for each column of Y on each column of X
univariate_regression = function(X, y, Z=NULL, return_residue=FALSE) {
  if (!is.null(Z)) {
    y = .lm.fit(Z, y)$residuals
  }
  calc_stderr = function(X, residuals) { sqrt(diag(sum(residuals^2) / (nrow(X) - 2) * chol2inv(chol(t(X) %*% X)))) }
  output = do.call(rbind,
                   lapply(c(1:ncol(X)), function(i) {
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
