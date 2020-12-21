library(abind)
mm_regression = function(X, Y, Z=NULL,center=TRUE,scale=TRUE) {
  if (!is.null(Z)) {
      Z = as.matrix(Z)
  }
  if(any(is.na(Y))){
    reg = lapply(seq_len(ncol(Y)), function (i) simplify2array(susieR:::univariate_regression(X, Y[,i], Z, center, scale)))
    reg = do.call(abind, c(reg, list(along=0)))
    # return array: out[1,,] is betahat, out[2,,] is shat
    out = aperm(reg, c(3,2,1))
    out = list(bhat = out[1,,], sbhat=out[2,,])
  }else{
    out = univariate_regression_Y(X, Y, Z, center, scale)
  }
  
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

library(data.table)
mm_regression_plink = function(geno.file, Y, Y_sample, prefix) {
  Y_names = colnames(Y)
  Y = cbind(Y_sample, Y_sample, Y)
  if(is.null(Y_names)){
    Y_names = 1:(ncol(Y)-2)
  }
  colnames(Y) = c('FID', 'IID', Y_names)
  Y = as.matrix(Y)
  write.table(Y, paste0(prefix, 'pheno'), quote=FALSE, row.names = FALSE)
  cmd = paste0('/project2/mstephens/software/plink-2.0/plink2 --pfile ', geno.file, 
               ' --glm hide-covar no-x-sex omit-ref --pheno ', paste0(prefix, 'pheno'),
               ' --threads 4 --memory 28000 ',
               ' --out ', prefix)
  system(cmd)
  
  bhat = c()
  sbhat = c()
  for(i in 1:(ncol(Y)-2)){
    out = fread(paste0(prefix,'.',Y_names[i],'.glm.linear'))
    bhat = cbind(bhat, out$BETA)
    sbhat = cbind(sbhat, out$SE)
    file.remove(paste0(prefix,'.',Y_names[i],'.glm.linear'))
  }
  rownames(bhat) = rownames(sbhat) = out$ID
  colnames(bhat) = colnames(sbhat) = Y_names
  file.remove(paste0(prefix, 'pheno'))
  file.remove(paste0(prefix, '.log'))
  return(list(bhat=bhat, sbhat=sbhat))
}

univariate_regression_Y = function (X, Y, Z = NULL, center = TRUE,
                                  scale = TRUE) {
  # no missing values in X, Y, Z
  n = nrow(X)
  if(!('scaled:scale' %in% names(attributes(X)))){
    cm = colMeans(X,na.rm = TRUE)
    csd = susieR:::compute_colSds(X)
    csd[csd == 0] = 1
    attr(X,"scaled:center") = cm
    attr(X,"scaled:scale") = csd
  }
  cm = attr(X,"scaled:center")
  csd = attr(X,"scaled:scale")

  if (!center)
      cm = rep(0,length = length(cm))
  if (!scale) 
      csd = rep(1,length = length(cm))
  if (center) {
    Y = t(t(Y) - colMeans(Y))
  } 

  if (!is.null(Z)) {
    if (center)
      Z = scale(Z,center = TRUE,scale = scale)
    ZtY = crossprod(Z, Y)
    dz = colSums(Z^2)
    Zcoef = ZtY/dz
    Y = Y - Z %*% Zcoef
  }
  XtY = crossprod(X, Y)/csd
  dx = (colSums(X^2) - nrow(X) * cm^2)/(csd^2)
  betahat = XtY/dx
  sb = betahat/csd
  sebetahat = do.call(cbind, lapply(1:ncol(Y), function(i){
    yhat = t(t(X)*sb[,i]) # N by J
    if(center){
      yhat <- t(t(yhat) - colMeans(yhat))
    }
    Yresid = Y[,i] - yhat
    rm(yhat)
    sqrt(colSums(Yresid^2)/((n-2) *dx))
  }))
  # sebetahat = do.call(rbind, lapply(1:ncol(X), function(i){
  #   yhat = tcrossprod(X[,i], sb[i,])
  #   if(center){
  #     yhat <- t(t(yhat) - colMeans(yhat))
  #   }
  #   Yresid = Y - yhat
  #   sqrt(colSums(Yresid^2)/((n-2) *dx[i]))
  # }))
  return(list(bhat = betahat,sbhat = sebetahat))
}

mm_sufficient = function(X,Y,center=TRUE,scale=TRUE, computeXtX = FALSE){
  if(!('scaled:scale' %in% names(attributes(X)))){
    cm = colMeans(X,na.rm = TRUE)
    csd = susieR:::compute_colSds(X)
    csd[csd == 0] = 1
    attr(X,"scaled:center") = cm
    attr(X,"scaled:scale") = csd
  }
  cm = attr(X,"scaled:center")
  csd = attr(X,"scaled:scale")

  if(center){
    Y = sweep(Y, 2, colMeans(Y), '-')
  }else{
    cm = rep(0,length = length(cm))
  }
  if (!scale)
    csd = rep(1,length = length(cm))
  
  XtY = crossprod(X, Y)/csd
  YtY = crossprod(Y)
  N = nrow(Y)

  if(computeXtX){
    sx <- colSums(X)
    XtX <- (crossprod(X) - tcrossprod(sx) / nrow(X)) / tcrossprod(csd)
  }else{
    XtX <- NULL
  }
  return(list(XtX = XtX, XtY=XtY, YtY=YtY, N=N))
}

