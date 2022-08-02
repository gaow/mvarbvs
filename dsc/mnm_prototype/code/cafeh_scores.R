
cafeh_scores = function(sets, sets_purity, pip, trait_pip, true_coef) {
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
              total=total, valid=valid, size=size, purity=purity,
              pip = c(pip), trait_pip = t(trait_pip)))
}
