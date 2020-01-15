library(atlasqtl)
# see `?atlasqtl`
pat = meta$true_coef != 0
p0 = c(mean(colSums(pat)), 10)
result = atlasqtl::atlasqtl(Y = Y, X = X, p0 = p0, user_seed = DSC_REPLICATE)
# FIXME: this can be wrong
result$pip = as.vector(1 - apply(1 - result$gam_vb, 1, prod))