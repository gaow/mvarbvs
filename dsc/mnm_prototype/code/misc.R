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

center_scale <- function(X) return(susieR:::set_X_attributes(as.matrix(X), center=TRUE, scale = TRUE))

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
      if (length(unique(x))==1) return(T)
      else return(F)
    }

    ### Filter X matrix
    filter_X <- function(X, missing_rate_thresh, maf_thresh) {
        rm_col <- which(apply(X, 2, compute_missing) > missing_rate_thresh)
        if (length(rm_col)) X <- X[, -rm_col]
        rm_col <- which(apply(X, 2, compute_maf) < maf_thresh)
        if (length(rm_col)) X <- X[, -rm_col]
        rm_col <- which(apply(X, 2, is_zero_variance))
        if (length(rm_col)) X <- X[, -rm_col]
        return(mean_impute(X))
    }