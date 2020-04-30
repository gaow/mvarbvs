# this is a temp solution to the problem with atlasqtl where perfectly correlated variables are removed
# even though this patch attempts to address to the issue, unfortunately it failed to divide PIP by the number of occurances.
# https://github.com/hruffieux/atlasqtl/commit/627385114a0559fdc8f957254db045437aaa85fb
rmvd_recover = function(dat, all_snps) {
    # first column is remained, 2nd removed
    coll_x = cbind(names(dat$rmvd_coll_x), dat$rmvd_coll_x)
    rownames(coll_x) = coll_x[,2]
    # make mapping -- remained : removed
    res = list()
    for (i in 1:nrow(coll_x)) {
        if (is.null(res[[coll_x[i,1]]])) res[[coll_x[i,1]]] = coll_x[i,2]
        else res[[coll_x[i,1]]] = append(res[[coll_x[i,1]]], coll_x[i,2])
    }
    # count for each remained how many removed is coll with it; plus itself
    for (name in names(res)) {
        res[[name]] = length(res[[name]]) + 1
    }
    # create full result with removed being NA
    new_pip = matrix(all_snps,length(all_snps),1)
    colnames(new_pip) = c("snp")
    reported_pip = cbind(rownames(dat$gam_vb), dat$gam_vb)
    colnames(reported_pip) = c("snp", colnames(dat$gam_vb))
    d = merge(new_pip, reported_pip, by='snp', all = T)
    if (!all(d$snp == all_snps)) stop("Merged table rows are messed up")
    ## life is hard ...
    d$snp = NULL
    d = as.data.frame(lapply(d, function(x) as.numeric(as.character(x))))
    rownames(d) = all_snps
    # fill PIPs for the removed
    for (i in 1:nrow(d)) {
        if (all(is.na(d[i,]))) {
            snp = as.character(all_snps[i])
            # an invariant site
            if (snp %in% dat$rmvd_cst_x) {
                d[i,] = 0
            } else {
                coll_snp = coll_x[snp,][1]
                d[i,] = d[coll_snp,] / res[[coll_snp]]
            }
        }
    }
    # adjust PIP for the remained
    for (coll_snp in names(res)) {
        d[coll_snp,] = d[coll_snp,] / res[[coll_snp]]
    }
    return(as.matrix(d))
}

# prior are mean and variance of expected number of effect variables
# set variance to 9 to indicate \pm 3 and use mean of actual non-zero effects
# I think this is quite fair
p0 = c(mean(colSums(meta$true_coef != 0)), 9)
result = atlasqtl::atlasqtl(Y = Y, X = X, p0 = p0, user_seed = DSC_REPLICATE)
result$gam_vb_completed = rmvd_recover(result, colnames(X))
# FIXME: have to discuss this
result$pip = as.vector(1 - apply(1 - result$gam_vb_completed, 1, prod))
