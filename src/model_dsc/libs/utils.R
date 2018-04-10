## Perform univariate regression for each column of Y on each column of X
univariate_regression = function(X, y, Z = NULL){
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
  return(list(betahat = output[,1], sebetahat = output[,2],
              residuals = y))
}
