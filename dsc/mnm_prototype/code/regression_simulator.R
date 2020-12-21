# input is prob of having each of n effect variables
# output is the number of effect variables for current data-set
get_n_signal = function(p = c(0.3,0.3,0.2,0.1,0.1,0,0,0,0,0)) {
    one_hot = rmultinom(1, 1, prob=p)[,1]
    return (which(one_hot == 1))
}

# input:
# J: number of total effects
# R: number of total conditions
# n: number of effect variables
# U: prior mixture object
# w: weight for each component
get_effects = function(J, R, n, U, w) {
    # Generate B under multivariate normal mixture
    # beta ~ \pi_0\delta_0 + \sum \pi_i N(0, U_i)
    if (J<n) stop("J should be greater than n")
    b = matrix(0, J, R)
    which_b = sample(1:J,n)
    Ub = list()
    for (j in which_b) {
        # each effect comes from a mixture component with given prob
        dist_index = which(rmultinom(1, 1, prob=w)[,1] == 1)
        b[j,] = MASS::mvrnorm(n = 1, mu=rep(0,R), Sigma=U[[dist_index]])
        # keep track of which matrix the effect came from
        Ub[[length(Ub) + 1]] = U[[dist_index]]
    }
    return(list(b=b,Ub=Ub))
}

get_residual_correlation = function(residual, R) {
    if (is.null(residual) || nrow(residual) != R) {
        return(diag(R))
    } else {
        return(cov2cor(residual))
    }
}

# is_pve_total: the input PVE is the total PVE for all variables
# If set to FALSE then it is per variable PVE, and total PVE for that effect is number of effect variables * the PVE.
get_y = function(X, b, residual_corr, pve, is_pve_total=FALSE, max_pve=0.8) {
    if('scaled:scale' %in% names(attributes(X)) ){
      sb = b/attr(X,"scaled:scale")
      Xsb <- X %*% sb
      yhat <- t(t(Xsb) - colMeans(Xsb))
    }else{
      yhat = X %*% b
    }
    genetic_var = apply(yhat, 2, var)
    if (is_pve_total) {
        pve = rep(pve, ncol(b))
    } else {
        pve = pve * apply(b, 2, function(x) length(which(x!=0)))
        pve[which(pve>max_pve)] = max_pve
    }
    # this is vector of residual variance scalar
    sigma = sqrt(genetic_var / pve - genetic_var)
    sigma[which(is.infinite(sigma) | is.nan(sigma))] = 0
    # some sigma can be zero due to lack of effect variables in that column of b
    # we have to set the corresponding residual variance to a reasonable default number
    # maybe the smallest of the rest of sigma?
    if (all(sigma == 0)) {
        sigma = rep(1, length(sigma))
    } else {
        sigma[which(sigma == 0)] = min(sigma[which(sigma!=0)])
    }
    residual_var = t(t(residual_corr * sigma) * sigma)
    y = yhat + MASS::mvrnorm(n = nrow(X), mu=rep(0,nrow(residual_var)), Sigma=residual_var)
    return(list(y = y, residual_var = residual_var))
}

get_prior = function(U,prior) {
    return(list(oracle=list(xUlist = U, pi = prior$w, null_weight = 0),
        identity = list(xUlist = list(identity=diag(nrow(U[[1]]))), pi=1, null_weight=0),
        shared = list(xUlist = list(matrix(1,nrow(U[[1]]),nrow(U[[1]]))), pi=1, null_weight=0),
        naive = list(xUlist = mmbr:::create_cov_canonical(nrow(U[[1]]))),
        ED = list(xUlist = prior$ED$U, pi = prior$ED$w, null_weight=0)))
}

# main simulation function where
# X is genotype, 
# U and w are mash mixture prior components and weights
# pve is per variable pve or total pve, see `is_pve_total` in `get_y`
# n is the number of true effects, default to NULL to use my default setting (see `get_n_signal()`)
# residual is diagonal by default value NULL
# scale_y is a boolean indicating whether or not the output variable y is to be scaled or not
mash_sim = function(X, U, w, pve, is_pve_total=FALSE, n=NULL, residual=NULL,scale_y=TRUE) {
    if (is.null(n) || n < 1) n = get_n_signal()
    R = nrow(U[[1]])
    if (R == 1) stop("This simulator is not meant for univariate data")
    for (i in 2:length(U)) {
        if (nrow(U[[i]])!=R) stop("Prior dimension are inconsistent in the given mixture")
    }
    effects = get_effects(ncol(X), R, n, U, w)
    b = effects$b
    residual_corr = get_residual_correlation(residual, R)
    y_sim = get_y(X, b, residual_corr, pve, is_pve_total)
    if (scale_y) {
        # y * diag(1/sd_y) = x * (b * diag(1/sd_y)) + (e * diag(1/sd_y))
        sd_y = apply(y_sim$y, 2, sd)
        y_sim$y = t(t(y_sim$y) / sd_y)
        y_sim$residual_var = t(y_sim$residual_var / sd_y) / sd_y
        b = t(t(b) / sd_y)
        for (i in 2:length(U)) {
            # b ~ N(0, diag(1/sd_y) * U * diag(1/sd_y))
            U[[i]] = t(U[[i]] / sd_y) / sd_y
        }
    } else {
        sd_y = NA
    }
    # FIXME: format data for R = 1
    return(list(Y=y_sim$y, true_coef=b, n_signal=n, n_traits=R, 
                residual_variance=y_sim$residual_var, U=U, 
                true_U=effects$Ub, Y_sd=sd_y))
}

simulate_main = function(X, Y, missing_Y, scale_Y, prior_file, prior, pve, is_pve_total, n_signal, var_Y, residual_mode, save_summary_stats, plink, prefix='data', save_suff_stats) {
    if(!is.matrix(X)){
      geno.file = X
      X = get_genotype(geno.file)
    }
    prior_data = readRDS(prior_file)
    if (residual_mode == 'identity') residual = NULL
    else residual = var_Y
    res = mash_sim(X, prior_data[[prior]]$U, prior_data[[prior]]$w, pve, is_pve_total, n_signal, residual, scale_Y)
    if (missing_Y && !is.null(Y)) res$Y = create_missing(res$Y, Y)
    res$prior = get_prior(res$U, prior_data[[prior]])
    if (save_summary_stats){
      if(plink){
        res$sumstats = mm_regression_plink(geno.file, res$Y, attr(X, 'sample_names'), prefix)
      }else{
        res$sumstats = mm_regression(X, res$Y)
      }
    }
    if (save_suff_stats) res$suff = mm_sufficient(X, res$Y) # XtX is not involved
    return(res)
}

# A wrapper function to mash_sim to work with Fabio's mr_mash DSC code
# missing_Y: whether or not to create missing values in Y based on original Y input
# prior_file: RDS file for prior database 
# prior: a string of a key in the RDS file to extract the prior to be used
# See `mash_sim` for documentation on other input parameters
mr_mash_sim = function(X, Y, missing_Y, scale_Y, prior_file, prior, n_signal, is_pve_total, var_Y, residual_mode) {
    prior_data = readRDS(prior_file)
    if (residual_mode == 'identity') residual = NULL
    else residual = var_Y
    res = mash_sim(X,prior_data[[prior]]$U, prior_data[[prior]]$w, pve, n_signal, is_pve_total, residual, scale_Y)
    if (missing_Y && !is.null(Y)) res$Y = create_missing(res$Y, Y)
    causal_variables = which(apply(res$true_coef, 1, function(x) any(x) != 0))
    causal_responses = list()
    for (i in 1:length(causal_variables)) {
        causal_responses[[i]] = which(res$true_coef[causal_variables[i],] != 0)
    }
    return(list(X=X, Y=res$Y, 
              B=res$true_coef, 
              V=res$residual_variance, 
              intercepts=rep(0,ncol(res$Y)), # because we have not simulated intercept term in Y. 
              causal_variables= causal_variables, 
              causal_responses= causal_responses,
              Sigma=res$true_U, 
              Gamma=NA))
}
