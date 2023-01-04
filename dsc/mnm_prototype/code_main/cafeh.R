# cafeh

run_cafeh <- function(bhat, se, LD_file, n, prefix) {
  cfg = list(bhat=paste0(prefix,".bhat"),
             shat=paste0(prefix,".shat"),
             LD=paste0(prefix,".ld"),
             n=paste0(prefix,".samplesize"))
  
  n = matrix(n, ncol=ncol(bhat), nrow=nrow(bhat))
  colnames(bhat) = colnames(se) = colnames(n) = paste0('t', 1:ncol(bhat))
  rownames(bhat) = rownames(se) = rownames(n) = paste0('rs', 1:nrow(bhat))
  write.table((t(bhat)), cfg$bhat, sep = "\t", quote = FALSE)
  write.table((t(se)), cfg$shat, sep = "\t", quote = FALSE)
  LD = readRDS(LD_file)
  rownames(LD) = colnames(LD) = rownames(bhat)
  write.table(LD, cfg$LD, sep = "\t", quote = FALSE)
  write.table(t(n), cfg$n, sep = "\t", quote = FALSE)

  return(cfg)
}


