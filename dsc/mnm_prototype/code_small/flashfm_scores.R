
flashfm_scores = function(m, true_coef, ld) {
  J = nrow(true_coef)
  if (is.character(ld)) {
    LD = readRDS(ld)
  } else {
    LD = ld
  }
  groups = m$snpGroups$groups.flashfm
  groups = lapply(groups, function(x) as.numeric(gsub('rs', '', gsub('id', '', x))))
  groups_purity = lapply(groups, function(x){
    if (length(x) == 1){
      return(1)
    }else{
      value = abs(Rfast::upper_tri(LD[x,x]))
      min(value)
    }
  })

  condition_pip = lapply(m$mpp.pp$MPP, function(x){
    pip = x[,2]
    full_pip = numeric(J)
    full_pip[as.numeric(gsub('rs', '', gsub('id', '', names(pip))))] = pip
    return(full_pip)
  })
  condition_pip = matrix(unlist(condition_pip), ncol=length(condition_pip))
  
  condition_groups = matrix(NA, length(groups), ncol(true_coef)) # L by R
  condition_group_status = matrix(NA, length(groups), ncol(true_coef)) # L by R, 9 for FP, 1 for TP, -1 for FN, 0 for TN
  condition_group_size = matrix(NA, length(groups), ncol(true_coef))
  condition_group_purity = matrix(NA, length(groups), ncol(true_coef))
  rownames(condition_groups) = rownames(condition_group_status) = 
    rownames(condition_group_size) = rownames(condition_group_purity) = names(groups)
  
  for(r in 1:ncol(true_coef)){
    for(gp in names(groups)){
      if(gp %in% rownames(m$mpp.pp$MPPg[[r]])){
        condition_groups[gp, r] = m$mpp.pp$MPPg[[r]][gp, 2] > 0.95
      }else{
        condition_groups[gp, r] = 0
      }
      if (condition_groups[gp, r] == 1) {
        condition_group_size[gp, r] = length(groups[[gp]])
        condition_group_purity[gp, r] = groups_purity[[gp]]
      }
      if (condition_groups[gp, r] == 1 &&
          any(groups[[gp]] %in% which(true_coef[, r] != 0))) {
        # TP
        condition_group_status[gp, r] = 1
      } else if (condition_groups[gp, r] == 1 &&
                 !any(groups[[gp]] %in% which(true_coef[, r] != 0))) {
        # FP
        condition_group_status[gp, r] = 9
      } else if (condition_groups[gp, r] == 0 &&
                 any(groups[[gp]] %in% which(true_coef[, r] != 0))) {
        # FN
        condition_group_status[gp, r] = -1
      } else {
        # TN
        condition_group_status[gp, r] = 0
      }
    }
  }
  
  total_condition_group = colSums(condition_groups) # length R vector
  false_positive_condition_group = colSums(condition_group_status == 9, na.rm = T)
  true_positive_condition_group = colSums(condition_group_status == 1, na.rm = T)
  false_neg_condition_group = colSums(condition_group_status == -1, na.rm = T)
  condition_group_size = mtx_to_list(condition_group_size) # length R list
  condition_group_purity = mtx_to_list(condition_group_purity)

  return(
    list(
      n_signal = sum(rowSums(true_coef != 0) > 0),
      n_condition_signal = colSums(true_coef != 0),
      total_cond_discoveries = total_condition_group,
      false_pos_cond_discoveries = false_positive_condition_group,
      true_cond_discoveries = true_positive_condition_group,
      false_neg_cond_discoveries = false_neg_condition_group,
      size_cond_cs = condition_group_size,
      purity_cond_cs = condition_group_purity,
      condition_pip = condition_pip
    )
  )
}