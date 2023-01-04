#' CAVIAR I/O
write_mscaviar_cfg <- function(z, prefix) {
  cfg = list(z=paste0(prefix,".z.txt"),
             ld=paste0(prefix,".ldfiles.txt"),
             set=paste0(prefix,"_set"),
             post=paste0(prefix,"_post"),
             log=paste0(prefix,".log"))
  filenames = numeric(ncol(zscore))
  for(r in 1:ncol(zscore)){
    filenames[r] = paste0(prefix, r,".z")
    write.table(zscore[,r], filenames[r], quote = F, col.names=F)
  }
  write.table(filenames, cfg$z, quote = F, col.names=F, row.names = F)
  write.table(rep(LD_file, ncol(zscore)), cfg$ld, quote = F, col.names=F, row.names = F)
  return(cfg)
}

#' Run CAVIAR
#' https://github.com/fhormoz/caviar

run_mscaviar <- function(z, LD_file, args = "", prefix="data")
{
  cfg = write_caviar_sumstats(z, prefix)
  cmd = paste("MsCAVIAR", "-z", cfg$z, "-l", cfg$ld, "-o", prefix, args)
  dscrutils::run_cmd(cmd)
  if(!all(file.exists(cfg$post, cfg$set, cfg$log))) {
      stop("Cannot find one of the post, set, and log files")
  }
  
  # remove unused files
  file.remove(cfg$z)
  
  log <- readLines(cfg$log)

  # read output tables
  snp <- read.delim(cfg$post)  
  stopifnot(ncol(snp) == 3)
  names(snp) <- c("snp", "snp_prob_set", "snp_prob")
  snp$snp <- as.character(snp$snp)
  snp <- rank_snp(snp)

  # `set` of snps
  set <- readLines(cfg$set)
  set_ordered <- left_join(data_frame(snp = set), snp, by = "snp") %>% 
    arrange(rank) %$% snp
  return(list(snp=snp, set=set_ordered))
}

rank_snp <- function(snp) {
  snp <- arrange(snp, -snp_prob) %>%
    mutate(
        rank = seq(1, n()),
        snp_prob_cumsum = cumsum(snp_prob) / sum(snp_prob)) %>%
    select(rank, snp, snp_prob, snp_prob_cumsum, snp_prob_set)
  return(snp)    
}

finemap_mcaviar <- function(zscore, LD_file, args, prefix, parallel = FALSE) {
  if (is.null(dim(zscore))) {
      zscore = matrix(ncol=1,zscore)
  }
  filenames = numeric(ncol(zscore)) 
  for(r in 1:ncol(zscore)){ 
    filenames[r] = paste0(prefix, r,".z")
    write.table(zscore[,r], filenames[r], quote = F, col.names=F, sep='\t')
  }
  write.table(filenames, paste0(prefix,".z.txt"), quote = F, col.names=F, row.names = F)
  write.table(rep(LD_file, ncol(zscore)), paste0(prefix,".ldfiles.txt"), quote = F, col.names=F, row.names = F)
  # write.table(rep('test.ld', ncol(zscore)), paste0(prefix,".ldfiles.txt"), quote = F, col.names=F, row.names = F)
  run_mscaviar(paste0(prefix,".z.txt"), paste0(prefix,".ldfiles.txt"), n=n)
}
