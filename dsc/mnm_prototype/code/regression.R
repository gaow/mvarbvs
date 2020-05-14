library(abind)
mm_regression = function(X, Y, Z=NULL,center=TRUE,scale=TRUE) {
  if (!is.null(Z)) {
      Z = as.matrix(Z)
  }
  reg = lapply(seq_len(ncol(Y)), function (i) simplify2array(susieR:::univariate_regression(X, Y[,i], Z, center, scale)))
  reg = do.call(abind, c(reg, list(along=0)))
  # return array: out[1,,] is betahat, out[2,,] is shat
  out = aperm(reg, c(3,2,1))
  out = list(bhat = out[1,,], sbhat=out[2,,])
  if (!is.null(colnames(X))) {
    rownames(out$bhat) = colnames(X)
    rownames(out$sbhat) = colnames(X)
  }
  if (!is.null(colnames(Y))) {
    colnames(out$bhat) = colnames(Y)
    colnames(out$sbhat) = colnames(Y)
  }
  return(out)
}