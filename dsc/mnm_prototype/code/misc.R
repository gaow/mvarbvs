get_center <- function(k,n) {
  ## For given number k, get the range k surrounding n/2
  ## but have to make sure it does not go over the bounds
  if (is.null(k)) {
      return(1:n)
  }
  start = floor(n/2 - k/2)
  end = floor(n/2 + k/2)
  if (start<1) start = 1
  if (end>n) end = n
  return(start:end)
}

subset_N <- function(gene, N_sub){
    ###If we want to subset individuals
    if(!is.null(N_sub)){
        Ntot <- nrow(gene$X)
        to_keep <- sort(sample(x=c(1:Ntot), size=N_sub, replace=F))

        gene$X <- gene$X[to_keep, ]
        gene$y <- gene$y[to_keep]
        gene$Z <- gene$Z[to_keep, ]

        return(gene)
    } else { ###If we don't want to subset individuals
        return(gene)
    }
}

compute_cov_diag <- function(Y){
    covar <- diag(apply(Y, 2, var, na.rm=T))
    return(covar)
}

compute_cov_flash <- function(Y, error_cache = NULL){
    covar <- diag(ncol(Y))
    tryCatch({
    fl <- flashier::flash(Y, var.type = 2, prior.family = c(flashier::prior.normal(), flashier::prior.normal.scale.mix()), backfit = TRUE, verbose.lvl=0)
    if(fl$n.factors==0){
      covar <- diag(fl$residuals.sd^2)
    } else {
      fsd <- sapply(fl$fitted.g[[1]], '[[', "sd")
      covar <- diag(fl$residuals.sd^2) + crossprod(t(fl$flash.fit$EF[[2]]) * fsd)
    }
    if (nrow(covar) == 0) {
      covar <- diag(ncol(Y))
      stop("Computed covariance matrix has zero rows")
    }
    }, error = function(e) {
      if (!is.null(error_cache)) {
        saveRDS(list(data=Y, message=warning(e)), error_cache)
        warning("FLASH failed. Using Identity matrix instead.")
        warning(e)
      } else {
        stop(e)
      }
    })
    s <- apply(Y, 2, sd, na.rm=T)
    if (length(s)>1) s = diag(s)
    else s = matrix(s,1,1)
    covar <- s%*%cov2cor(covar)%*%s
    return(covar)
}

create_missing <- function(Y1, Y0) {
    if (ncol(Y0) < ncol(Y1)) {
        for (i in 1:(ncol(Y1) - ncol(Y0))) {
            Y0 = cbind(Y0, sample(Y0[, sample.int(ncol(Y0), size=1)]))
        }
    }
    if (ncol(Y0) > ncol(Y1)) Y0 = Y0[,sample(1:ncol(Y1))]
    res = Y1
    res[which(is.na(Y0))] = NA
    # it is possible that some rows of Y1 are made all NA
    # here we have to make up for it
    na_rows = which(apply(res, 1, function(x) all(is.na(x))))
    for (i in na_rows) {
	    non_na = sample.int(ncol(res), size=1)
        res[i,non_na] = Y1[i,non_na]
    }
    return(res)
}

###Functions to compute MAF and missing genotype rate
compute_maf <- function(geno){
   f <- mean(geno,na.rm = TRUE)/2
   return(min(f, 1-f))
}

compute_missing <- function(geno){
  miss <- sum(is.na(geno))/length(geno)
  return(miss)
}

mean_impute <- function(geno){
  f <- apply(geno, 2, function(x) mean(x,na.rm = TRUE))
  for (i in 1:length(f)) geno[,i][which(is.na(geno[,i]))] <- f[i]
  return(geno)
}

is_zero_variance <- function(x) {
  if (length(unique(x[!is.na(x)]))==1) return(T)
  else return(F)
}

### Filter X matrix
filter_X <- function(X, missing_rate_thresh, maf_thresh) {
    rm_col <- which(apply(X, 2, compute_missing) > missing_rate_thresh)
    if (length(rm_col)) X <- X[, -rm_col]
    rm_col <- which(apply(X, 2, compute_maf) < maf_thresh)
    if (length(rm_col)) X <- X[, -rm_col]
    X <- mean_impute(X)
    rm_col <- which(apply(X, 2, is_zero_variance))
    if (length(rm_col)) X <- X[, -rm_col]
    return(X)
}

get_genotype <- function(geno_file){
  library(data.table)
  library(Matrix)
  geno <- fread(paste0(geno_file, '.raw.gz'),sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  class(geno) <- "data.frame"
  # Extract the genotypes.
  X <- as(as.matrix(geno[-(1:6)]), 'dgCMatrix')
  cm = colMeans(X,na.rm = TRUE)
  csd = susieR:::compute_colSds(X)
  csd[csd == 0] = 1
  attr(X,"scaled:center") = cm
  attr(X,"scaled:scale") = csd
  attr(X,"sample_names") = geno[,2]
  return(X)
}


