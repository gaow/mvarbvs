# PAINTOR
PAINTOR: paintor.R + R(b = $(meta)$true_coef;
                       nc = sum(rowSums($(meta)$true_coef!=0) > 0);
                       if(nc > 2) nc = 2;
                       args = paste0('-enumerate ', nc);
                       Z = as.matrix(sumstats$bhat/sumstats$sbhat);
                       Z[is.na(Z)] = 0;
                       pip = finemap_paintor(Z, LD_file, args, prefix=cache))
  sumstats: $sumstats  
  LD_file: $ld
  cache: file(PAINTOR)
  $pip: pip
  
