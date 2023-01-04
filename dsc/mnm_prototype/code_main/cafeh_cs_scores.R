
cafeh_cs_scores = function(sets, sets_purity, true_coef) {
  if (is.null(dim(true_coef))) beta_idx = which(true_coef!=0)
  else beta_idx = which(apply(true_coef, 1, sum) != 0)
  
  pure_cs_idx = which(sets_purity > 0.5)
  if(length(pure_cs_idx) == 0){
    size = NA
    total = 0
    purity = NA
  }else{
    pure_cs = sets[pure_cs_idx]
    pure_cs_purity = sets_purity[pure_cs_idx]
    
    size = sapply(pure_cs,length)
    purity = unlist(pure_cs_purity)
    total = length(pure_cs)
  }
  valid = 0
  if (total > 0) {
    for (i in 1:total) {
      if (any((pure_cs[[i]]+1) %in% beta_idx))
        valid = valid + 1
    }
  }
  
  return(list(n_signal = length(beta_idx),
              n_condition_signal = colSums(true_coef != 0),
              total=total, valid=valid, size=size, purity=purity))
}

cafeh_cs_scores_multiple = function(sets, sets_purity, true_coef){
  total = valid = size = purity = list()
  for (idx in 1:length(sets)){
    out = cafeh_cs_scores(sets[[idx]], sets_purity[[idx]], true_coef)
    total[[idx]] = out$total
    valid[[idx]] = out$valid
    size[[idx]] = out$size
    purity[[idx]] = out$purity
  }
  names(total) = names(valid) = names(size) = names(purity) = names(sets)
  n_signal = out$n_signal
  n_condition_signal = out$n_condition_signal
  return(list(n_signal = n_signal,
              n_condition_signal = n_condition_signal,
              total=total, valid=valid, size=size, purity=purity))
}
