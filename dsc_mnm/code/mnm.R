compute_cov_diag <- function(Y){
    covar <- diag(apply(Y, 2, var, na.rm=T))
    return(covar)
}

compute_cov_flash <- function(Y){
    fl <- flashier::flash(Y, var.type = 2, prior.family = c(flashier::prior.normal(), flashier::prior.normal.scale.mix()), backfit = TRUE, verbose.lvl=0)
    if(fl$n.factors==0){
      covar <- diag(fl$residuals.sd^2)
    } else {
      fsd <- sapply(fl$fitted.g[[1]], '[[', "sd")
      covar <- diag(fl$residuals.sd^2) + crossprod(t(fl$flash.fit$EF[[2]]) * fsd)
    }
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

prior = cfg[[as.character(ncol(Y))]][[eff_mode]]
if (missing_Y) Y = create_missing(Y, meta$original_Y)
if (resid_method == 'flash') {
    resid_Y <- compute_cov_flash(Y)
} else if (resid_method == 'diag') {
    resid_Y <- compute_cov_diag(Y)
} else {
    resid_Y <- meta$residual_variance 
}
m_init = mmbr:::MashInitializer$new(NULL, NULL, xUlist=append(list(matrix(0,ncol(Y),ncol(Y))), prior$xUlist), prior_weights=prior$pi, null_weight=prior$null_weight, alpha=alpha, top_mixtures=-1)
result = mmbr::msusie(X, Y, L=L, prior_variance=m_init, residual_variance=resid_Y, compute_objective=T, estimate_residual_variance=F, estimate_prior_variance=F)
#result$pip_conditions = mmbr::mmbr_get_pip_per_condition(result, m_init)