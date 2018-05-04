write_finemap_sumstats <- function(z,LD_file,k,prefix='data') {
  n = nrow(y)
  cfg = list(z=paste0(prefix,".z"),
             ld=LD_file,
             snp=paste0(prefix,".snp"),
             config=paste0(prefix,".config"),
             k=paste0(prefix,".k"),
             log=paste0(prefix,".log"),
             meta=paste0(prefix,".manifest"))
  write.table(z,cfg$z,quote=F,col.names=F)
  write.table(t(k),cfg$k,quote=F,col.names=F,row.names=F)
  write("z;ld;snp;config;k;log;n-ind",file=cfg$meta)
  write(paste(cfg$z, cfg$ld, cfg$snp, cfg$config, cfg$k, cfg$log, n, sep=";"),
        file=cfg$meta,append=TRUE)
  return(cfg)
}

#' @param k prior. set to -999 to use default prior
finemap <- function(X,y,z,LD_file,sa=0.4,k=c(0,0,0,1),niter=1000000,prefix='data') {
  ## default version puts all weight on k=4
  ## and gives results more similar to the VB approach
  ## set k = -999 to use default prior
  prior_k = ifelse(k == -999, '', '--prior-k')
  cfg = write_finemap_sumstats(X,y,z,LD_file,k,prefix)
  cmd = paste("finemap --sss --in-files", cfg$meta, prior_k,
              '--n-iterations', niter, '--prior-std', sa, '--regions 1')
  write(cmd, stderr())
  system(cmd)
  return(read.table(cfg$snp,header=TRUE,sep=" "))
}

finemapM <- function(zscore,ld,sa=0.4,k=c(0,0,0,1),niter=1000000,prefix='data') {
  LD_file <- paste0(prefix,".ld")
  write.table(ld,LD_file,quote=F,col.names=F,row.names=F)
  return(do.call(cbind, lapply(1:ncol(Y), function(r) finemap(zscore[,r],LD_file,sa,k,niter,prefix))))
}
