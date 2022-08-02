library(R2HESS)
if(!is.matrix(X)){
  geno.file = X
  X = get_genotype(geno.file)
}
cat('run config\n')
config <- r2hess.makeConfig(data.frame(as.matrix(X)), data.frame(as.matrix(Y)), data_dir, r=ncol(Y))
cat('run mthess\n')
r2hess.run(config)
pip_file = paste0(data_dir, '/Marg_Prob_Incl.txt')
result = list(pip_conditions = t(read.table(pip_file, skip=1)))