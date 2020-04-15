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
    for (j in which_b) {
        # each effect comes from a mixture component with given prob
        dist_index = which(rmultinom(1, 1, prob=w)[,1] == 1)
        b[j,] = MASS::mvrnorm(n = 1, mu=rep(0,R), Sigma=U[[dist_index]])
    }
    return(b)
}

get_residual_correlation = function(residual, R) {
    if (is.null(residual)) {
        return(diag(R))
    } else {
        if (nrow(residual) != R) stop("residual dimension mismatch")
        return(cov2cor(residual))
    }
}

# is_pve_variable_avg: the input PVE is the average PVE per variable.
# That is, total PVE for that effect is number of effect variables * the PVE.
get_y = function(X, b, residual_corr, pve, is_pve_variable_avg=TRUE, max_pve=0.8) {
    yhat = X %*% b
    genetic_var = apply(yhat, 2, var)
    if (is_pve_variable_avg) {
        pve = pve * apply(b, 2, function(x) length(which(x!=0)))
        pve[which(pve>max_pve)] = max_pve
    } else {
        pve = rep(pve, ncol(b))
    }
    # this is vector of residual variance scalar
    sigma = sqrt(genetic_var / pve - genetic_var)
    sigma[which(is.infinite(sigma))] = 0
    # some sigma can be zero due to lack of effect variables in that column of b
    # we have to set the corresponding residual variance to a reasonable default number
    # maybe the smallest of the rest of sigma?
    if (all(sigma == 0)) {
        sigma = rep(1, length(sigma))
    } else {
        sigma[which(sigma == 0)] = min(sigma[which(sigma!=0)])
    }
    sigma = diag(sigma)
    residual_var = sigma %*% residual_corr %*% sigma
    y = yhat + MASS::mvrnorm(n = nrow(X), mu=rep(0,nrow(residual_var)), Sigma=residual_var)
    return(list(y = y, residual_var = residual_var))
}

mash_sim = function(X, J, U, w, pve, n=NULL, residual=NULL) {
    if (is.null(n) || n < 1) n = get_n_signal()
    R = nrow(U[[1]])
    if (R == 1) stop("This simulator is not meant for univariate data")
    for (i in 2:length(U)) {
        if (nrow(U[[i]])!=R) stop("Prior dimension are inconsistent in the given mixture")
    }
    b = get_effects(J, R, n, U, w)
    residual_corr = get_residual_correlation(residual, R)
    y_sim = get_y(X, b, residual_corr, pve)
    # FIXME: format data for R = 1
    return(list(Y=y_sim$y, true_coef=b, n_signal=n, n_traits=R, residual_variance=y_sim$residual_var))
}

get_prior = function(U,w) {
    return(list(oracle=list(xUlist = U, pi = w, null_weight = 0),
    identity = list(xUlist = list(identity=diag(nrow(U[[1]]))), pi=1, null_weight=0),
    shared = list(xUlist = list(matrix(1,nrow(U[[1]]),nrow(U[[1]]))), pi=1, null_weight=0)))
}

simulate_main = function(X, prior_file, prior, n_signal, var_Y, residual_mode) {
    prior_data = readRDS(prior_file)
    X = susieR:::set_X_attributes(X)
    if (residual_mode == 'identity') residual = NULL
    else residual = var_Y
    res = mash_sim(X, ncol(X), prior_data[[prior]]$U, prior_data[[prior]]$w, pve, n_signal, residual)
    res$X = X
    res$prior = get_prior(prior_data[[prior]]$U, prior_data[[prior]]$w)
    res$X_mean = attributes(X)[["scaled:center"]]
    res$X_csd = attributes(X)[["scaled:scale"]]
    return(res)
}