#' @title Residual covariance for a M&M fit
#' B is J X R matrix of M&M output
#' SM is J X R X R matrix of M&M output with 2nd moment option on
#' XtX is just precomputed t(X) %*% X
compute_mnm_residual_covariance = function(X, Y, XtX, B, SM) {
    # out = t(Y) %*% Y - 2 * t(B) %*% t(X) %*% Y # + E[B^TX^TXB]
    # E[B^TX^TXB] is not easy to compute properly
    # use MLE for now
    return(t(Y - X%*%B) %*% (Y - X%*%B) / nrow(X))
} 
                                            
#' @title expected loglikelihood for a M&M fit
# https://gaow.github.io/mvarbvs/writeup/20171215_MNMModel_Finemap.html
#' S is simply XtX pre-computed
#' Sigma is current estimate of residual variance
#' B is J X R matrix of M&M output
#' SM is J X R X R matrix of M&M output with 2nd moment option on
compute_mnm_Eloglik = function(X,Y,S,Sigma,B,SM){
    inv_Sigma = solve(Sigma)
    det_Sigma = det(Sigma)
    N = nrow(Y)
    R = ncol(Y)
    t0 = vector()
    for (j in 1:ncol(X)) {
        t0[j] = S[j,j] * sum(inv_Sigma * SM[,,j])
    }
    t1 = sum(diag(inv_Sigma %*% t(B) %*% S %*% B)) + sum(t0) -
            2 * sum(diag(Y %*% inv_Sigma %*% t(B) %*% t(X))) +
            sum(diag(Y %*% inv_Sigma %*% t(Y)))
    out = -0.5 * N * R * log(2 * pi) - 0.5 * N * log(det_Sigma) - 0.5 * t1
    return(out)
}

#' @title posterior expected loglikelihood for a MASH problem
## E[log(\hat{B}|B, Shat)]
## Need posterior mean and posterior second moment from MASH
## do not use any computational trick here because this is 
## just for sanity check
compute_mash_Eloglik = function(betahat, Shat, b, b2) {
    inv_Shat = solve(Shat)
    det_Shat = det(Shat)
    res = nrow(b) * log(2*pi) + log(det_Shat) +
            (t(betahat) %*% inv_Shat %*% betahat -
            2 * t(betahat) %*% inv_Shat %*% b +
            sum(diag(inv_Shat %*% b2)))
    return(-0.5 * res)
}

#' @title sum of MASH posterior expected loglikelihood
#' Bhat is J x R matrix of MASH input
#' SDhat is J X R matrix to be expanded with V, turning into J X R X R
#' V is R X R matrix of MASH input
#' Sigma is residual variance
#' alpah is a J vector of weights
#' B is J X R matrix of MASH output
#' SM is J X R X R matrix of MASH output with 2nd moment option on
compute_sse_Eloglik = function(Bhat, SDhat, V, Sigma, alpha, B, SM) {
    ## FIXME: I think it is wrong here because it is not single effect model
    ## where J effects should NOT be factorized.
    ## But otherwise isn't it a matrix normal density with both row and column covariances?
    res = vector()
    for (j in 1:nrow(Bhat)) {
        ## Is R X R
        Shat = SDhat[j,] * t(V * SDhat[j,]) # faster than diag(SDhat[j,]) %*% V %*% diag(SDhat[j,])
        ## Is R X 1
        B_j = B[j,] * alpha[j]
        ## 2nd moment, R X R
        B2_j = (B[j,] %*% t(B[j,]) + SM[,,j]) * alpha[j]
        res[j] = compute_mash_Eloglik(Bhat[j,], Shat, B_j, B2_j)
    }
    return(sum(res))
}
