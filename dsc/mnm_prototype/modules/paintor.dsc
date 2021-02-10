# PAINTOR
PAINTOR: paintor.R + R(Z = as.matrix(sumstats$bhat/sumstats$sbhat);
                       Z[is.na(Z)] = 0;
                       pip = finemap_paintor(Z, LD_file, args, prefix=cache))
  sumstats: $sumstats  
  LD_file: $ld
  args: '-enumerate 2'
  cache: file(PAINTOR)
  $pip: pip
  
