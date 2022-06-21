library(KScorrect)
#' @title Check if produced confidence sets have overlaps
#' @param cs a list a susie confidence sets from susie fit
#' @return number of overlaps in terms of both SNPs and CS
check_overlap = function(cs) {
  if (length(cs) == 0) {
    return(list(snp = 0, cs = 0))
  } else {
    overlaps_snp = 0
    overlaps_cs = 0
    for (i in 1:length(cs)) {
      for (j in 1:i) {
        if (i == j) next
        overlap = intersect(cs[[i]], cs[[j]])
        overlaps_snp = overlaps_snp + length(overlap)
        overlaps_cs = overlaps_cs + (length(overlap) > 0)
      }
    }
    return(list(snp = overlaps_snp, cs = overlaps_cs))
  }
}

mtx_to_list = function(m){
  lapply(seq_len(ncol(m)), function(i){
    x = m[,i]
    idx = which(!is.na(x))
    if(length(idx) == 0){
      return(NA)
    }else{
      return(x[idx])
    }
  })
}

#' @title Compare SuSiE fits to truth
#' @param sets a list of susie CS info from susie fit
#' @param pip probability for p variables
#' @param true_coef true regression coefficients
#' @return total the number of total CS
#' @return valid the number of CS that captures a true signal
#' @return size an array of size of CS
#' @return purity an array of purity of CS
#' @return top the number of CS whose highest PIP is the true causal
mvsusie_scores = function(m, true_coef, lfsr_cutoff = 0.05) {
  sets = m$sets
  pip = m$pip
  if (is.null(dim(true_coef)))
    beta_idx = which(true_coef != 0)
  else
    beta_idx = which(rowSums(true_coef != 0) > 0)
  cs = sets$cs
  cs_significant = m$single_effect_lfsr < lfsr_cutoff # L by R
  condition_cs_status = matrix(NA, nrow(cs_significant), ncol(cs_significant)) # 9 for FP, 1 for TP, -1 for FN, 0 for TN
  condition_cs_size = matrix(NA, nrow(cs_significant), ncol(cs_significant))
  condition_cs_purity = matrix(NA, nrow(cs_significant), ncol(cs_significant))
  condition_cs_avgr2 = matrix(NA, nrow(cs_significant), ncol(cs_significant))
  cs_idx = unlist(sets$cs_index)
  if (is.null(cs)) {
    size = NA
    total = 0
    purity = NA
    avgr2 = NA
  } else {
    size = sapply(cs, length)
    purity = as.vector(sets$purity[, 1])
    avgr2 = as.vector(sets$purity[, 2]) ^ 2
    total = length(cs)
  }
  valid = 0
  top_hit = 0
  if (total > 0) {
    for (i in 1:total) {
      if (any(cs[[i]] %in% beta_idx))
        valid = valid + 1
      set.idx = cs[[i]]
      highest.idx = which.max(pip[set.idx])
      if (set.idx[highest.idx] %in% beta_idx)
        top_hit = top_hit + 1
    }
    # per condition effect estimate
    if (!(is.null(dim(true_coef)))) {
      for (r in 1:ncol(true_coef)) {
        for (e in 1:nrow(cs_significant)) {
          if (e %in% cs_idx) {
            cs_condition = cs[[which(cs_idx == e)]]
            if (cs_significant[e, r] == 1) {
              condition_cs_size[e, r] = length(cs_condition)
              condition_cs_purity[e, r] = sets$purity[which(cs_idx == e), 1]
              condition_cs_avgr2[e, r] = sets$purity[which(cs_idx == e), 2]^2
            }
            if (cs_significant[e, r] == 1 &&
                any(cs_condition %in% which(true_coef[, r] != 0))) {
              # TP
              condition_cs_status[e, r] = 1
            } else if (cs_significant[e, r] == 1 &&
                       !any(cs_condition %in% which(true_coef[, r] != 0))) {
              # FP
              condition_cs_status[e, r] = 9
            } else if (cs_significant[e, r] == 0 &&
                       any(cs_condition %in% which(true_coef[, r] != 0))) {
              # FN
              condition_cs_status[e, r] = -1
            } else {
              # TN
              condition_cs_status[e, r] = 0
            }
          }
        }
      }
    }
  }
  total_condition_cs = colSums(cs_significant) # length R vector
  false_positive_condition_cs = colSums(condition_cs_status == 9, na.rm = T)
  true_positive_condition_cs = colSums(condition_cs_status == 1, na.rm = T)
  false_neg_condition_cs = colSums(condition_cs_status == -1, na.rm = T)
  condition_cs_size = mtx_to_list(condition_cs_size) # length R list
  condition_cs_purity = mtx_to_list(condition_cs_purity)
  condition_cs_avgr2 = mtx_to_list(condition_cs_avgr2)
  
  # effect size estimates
  # for each true effect, look at its b1 and b2 estimate, sum them over, compute stderr,
  # and see the percentile of the true parameter with respect to the estimate
  estimate_diff = vector()
  if (!(is.null(dim(true_coef)))) {
    for (r in 1:ncol(true_coef)) {
      true_idx = which(true_coef[, r] != 0)
      truth = true_coef[, r][true_idx]
      if (length(truth) == 0)
        next
      mu = m$b1[, true_idx, r]
      mu2 = m$b2[, true_idx, r]
      if (is.null(dim(mu))) {
        mu = matrix(mu, 1, length(mu))
        mu2 = matrix(mu2, 1, length(mu2))
      }
      sd_est = sqrt(mu2 - mu ^ 2)
      quant_diff = 0
      for (i in 1:length(truth)) {
        quant_diff = quant_diff + abs(KScorrect::pmixnorm(truth[i], mu[i, ], sd_est[i, ], 
                                                          rep(1 /length(mu[i, ]), length(mu[i, ]))) - 0.5)
      }
      estimate_diff = append(estimate_diff, quant_diff / length(truth))
    }
  }
  overlaps = check_overlap(cs)
  cs_corr = ifelse(all(!is.na(m$cs_corr)), max(abs(m$cs_corr[upper.tri(m$cs_corr)])), NA)
  return(
    list(
      total = total,
      valid = valid,
      size = size,
      purity = purity,
      avgr2 = avgr2,
      top = top_hit,
      overlap_var = overlaps$snp,
      overlap_cs = overlaps$cs,
      n_signal = length(beta_idx),
      included_signal = sum(beta_idx %in% unlist(cs)),
      total_cond_discoveries = total_condition_cs,
      false_pos_cond_discoveries = false_positive_condition_cs,
      true_cond_discoveries = true_positive_condition_cs,
      false_neg_cond_discoveries = false_neg_condition_cs,
      size_cond_cs = condition_cs_size,
      purity_cond_cs = condition_cs_purity,
      avgr2_cond_cs = condition_cs_avgr2,
      avg_diff_eff_size_percentile = mean(estimate_diff),
      cs_correlation = cs_corr,
      converged = m$convergence$converged
    )
  )
}