#' PAINTOR I/O
library(data.table)
write_paintor_cfg <- function(Z, LD_file, prefix) {
  tmp = unlist(strsplit(prefix, '/'))
  dir = paste0(tmp[-length(tmp)], collapse = '/')
  cfg = list(input=paste0(prefix, '.files'),
             dir=dir,
             ld=paste0(prefix,".ld"),
             post=paste0(prefix,".results"),
             logBF=paste0(prefix,".LogBF"),
             enrichment=paste0(prefix,".Enrichment"))
  write.table(Z, prefix, quote = F, row.names=F, sep=' ')
  ld = readRDS(LD_file)
  write.table(ld, cfg$ld, quote=F,col.names=F,row.names=F, sep=' ')
  write.table(tmp[length(tmp)], cfg$input, quote=F,col.names=F,row.names=F)
  anno = data.frame(dummy = rep(1, nrow(Z)))
  write.table(anno, paste0(prefix, '.annotations'), quote=F,row.names=F)
  return(cfg)
}

#' Run PAINTOR
#' https://github.com/gkichaev/PAINTOR_V3.0

run_paintor <- function(Z, LD_file, args = "", prefix="data"){
  cfg = write_paintor_cfg(Z, LD_file, prefix)
  cmd = paste0("code/PAINTOR -input ",cfg$input, ' -in ', cfg$dir, ' -out ', cfg$dir, ' ',
               '-Zhead ', paste0(colnames(Z), collapse = ','), ' ', 
               '-LDname ', paste0(rep('ld', ncol(Z)), collapse = ','), ' ',
               '-annotations dummy -Lname ', cfg$logBF, ' -Gname ', cfg$enrichment, ' ',
               args)
  dscrutils::run_cmd(cmd)
  if(!file.exists(cfg$post)){
    stop("Cannot find post file")
  }
  
  # remove unused files
  file.remove(prefix)
  file.remove(cfg$input)
  file.remove(cfg$ld)
  file.remove(paste0(prefix, '.annotations'))
  
  # read output tables
  snp <- fread(cfg$post)
  pip = snp$Posterior_Prob

  return(pip)
}

finemap_paintor <- function(Z, LD_file, args, prefix, parallel = FALSE) {
  if (is.null(dim(Z))) {
    Z = matrix(ncol=1,Z)
  }
  pip = run_paintor(Z, LD_file, args, prefix)
  return(pip)
}
