load_data = function(genotype_file, expr_file, geno_table, expr_table) {
    geno <- h5read(genotype_file, geno_table)
    gdata = geno$block0_values
    colnames(gdata) = geno$axis1
    rownames(gdata) = geno$axis0
    expr <- h5read(expr_file, expr_table)
    edata = expr$block0_values
    colnames(edata) = tools::file_path_sans_ext(expr$axis1)
    rownames(edata) = apply(sapply(strsplit(expr$axis0, "-"), `[`, c(1,2)), 2, function(x) paste(x, collapse = '-'))
    edata = edata[, basename(geno_table)]
    return(list(X=gdata,y=edata))
}

univariate_lm = function(X,y,Z=NULL){
  P = dim(X)[2]
  output = matrix(0,nrow = P,ncol = 2)
  for(i in 1:P){
    if (is.null(Z))  g = summary(lm(y~X[,i]))
    else  g = summary(lm(y~X[,i] + Z))
    output[i,] = g$coefficients[2,1:2]
  }
  return(list(x = output[,1], s = output[,2], lik=list(name="normal")))
}
