
#' @title Compare SuSiE fits to truth
#' @param sets a list of susie CS info from susie fit
#' @param pip probability for p variables
#' @param true_coef true regression coefficients
#' @return total the number of total CS
#' @return valid the number of CS that captures a true signal
#' @return size an array of size of CS
#' @return purity an array of purity of CS
#' @return top the number of CS whose highest PIP is the true causal
susie_cs_scores = function(sets, pip, true_coef) {
  if (is.null(dim(true_coef))) beta_idx = which(true_coef!=0)
  else beta_idx = which(apply(true_coef, 1, sum) != 0)
  if(is.null(sets)){
    return(list(total=NA, valid=NA, size=NA, purity=NA, avgr2=NA))
  }
  cs = sets$cs
  if (is.null(cs)) {
    size = NA
    total = 0
    purity = NA
    avgr2 = NA
  } else {
    size = sapply(cs,length)
    purity = as.vector(sets$purity[,1])
    avgr2 = as.vector(sets$purity[,2])^2
    total = length(cs)
  }
  valid = 0
  if (total > 0) {
    for (i in 1:total){
      if (any(cs[[i]]%in%beta_idx)) valid=valid+1
    }
  }
  return(list(total=total, valid=valid, size=size, purity=purity, avgr2=avgr2))
}

susie_cs_scores_multiple = function(res, sets, truth) {
  total = valid = size = purity = list()
  objective = converged = expected = vector()
  for (r in 1:length(res)) {
    total_r = valid_r = size_r = purity_r = list()
    for (idx in 1:length(sets[[r]])){
      out = susie_cs_scores(sets[[r]][[idx]], res[[r]]$pip, truth[,r])
      total_r[[idx]] = out$total
      valid_r[[idx]] = out$valid
      size_r[[idx]] = out$size
      purity_r[[idx]] = out$purity
    }
    names(total_r) = names(valid_r) = names(size_r) = names(purity_r) = names(sets[[r]])
    total[[r]] = total_r
    valid[[r]] = valid_r
    size[[r]] = size_r
    purity[[r]] = purity_r
    
    expected[r] = sum(truth[,r] !=0)
    if(is.null(susieR::susie_get_objective(res[[r]]))){
      objective[r] = NA
      converged[r] = NA
    }else{
      objective[r] = susieR::susie_get_objective(res[[r]])
      converged[r] = res[[r]]$converged
    }
  }
  return(list(total=total, valid=valid, size=size, purity=purity, expected = expected,
              objective=objective, converged=converged))
}
