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

#' @title Compare SuSiE fits to truth
#' @param sets a list of susie CS info from susie fit
#' @param pip probability for p variables
#' @param true_coef true regression coefficients
#' @return total the number of total CS
#' @return valid the number of CS that captures a true signal
#' @return size an array of size of CS
#' @return purity an array of purity of CS
#' @return top the number of CS whose highest PIP is the true causal
susie_scores = function(m, true_coef, lfsr_cutoff = 0.05) {
  sets = m$sets
  pip = m$pip
  if (is.null(dim(true_coef))) beta_idx = which(true_coef != 0)
  else beta_idx = which(apply(true_coef, 1, sum) != 0)
  cs = sets$cs
  cs_significant = mmbr::mmbr_get_cs_lfsr(m) < lfsr_cutoff
  condition_cs_status = matrix(NA, nrow(cs_significant), ncol(cs_significant)) # 9 for FP, 1 for TP, -1 for FN, 0 for TN
  cs_idx = unlist(sets$cs_index)
  if (is.null(cs)) {
    size = 0
    total = 0
    purity = 0
  } else {
    size = sapply(cs,length)
    purity = as.vector(sets$purity[,1])
    total = length(cs)
  }
  valid = 0
  top_hit = 0
  if (total > 0) {
    for (i in 1:total) {
      if (any(cs[[i]] %in% beta_idx)) valid = valid+1
      set.idx = cs[[i]]
      highest.idx = which.max(pip[set.idx])
      if (set.idx[highest.idx] %in% beta_idx) top_hit = top_hit+1
    }
      # per condition effect estimate
      if (!(is.null(dim(true_coef)))) {
        for (r in 1:ncol(true_coef)) {
          for (e in 1:nrow(cs_significant)) {
            if (e %in% cs_idx) {
              cs_condition = cs[[which(cs_idx == e)]]
              if (cs_significant[e,r] == 1 && any(cs_condition %in% which(true_coef[,r] != 0))) {
                condition_cs_status[e,r] = 1
              } else if (cs_significant[e,r] == 1 && !any(cs_condition %in% which(true_coef[,r] != 0))) {
                condition_cs_status[e,r] = 9
              } else if (cs_significant[e,r] == 0 && length(which(true_coef[,r] != 0)) > 0) {
                condition_cs_status[e,r] = -1
              } else {
                condition_cs_status[e,r] = 0
              }
            }
          }
        }
      }
  }
  false_positive_condition_cs = length(which(condition_cs_status == 9))
  true_positive_condition_cs = length(which(condition_cs_status == 1))
  false_neg_condition_cs = length(which(condition_cs_status == -1))
  # effect size estimates
  # for each true effect, look at its b1 and b2 estimate, sum them over, compute stderr,
  # and see the percentile of the true parameter with respect to the estimate
  estimate_diff = vector()
  if (!(is.null(dim(true_coef)))) {
    for (r in 1:ncol(true_coef)) {
      true_idx = which(true_coef[,r] != 0)
      truth = true_coef[,r][true_idx]
      if (length(truth) == 0) next
      mu = m$b1[,true_idx,r]
      mu2 = m$b2[,true_idx,r]
      if (is.null(dim(mu))) {
        mu = matrix(mu, length(mu), 1)
        mu2 = matrix(mu2, length(mu2), 1)
      }
      sd_est = sqrt(mu2 - mu^2)
      quant_diff = 0
      for (i in 1:length(truth)) {
        quant_diff = quant_diff + abs(KScorrect::pmixnorm(truth[i], mu[i,], sd_est[i,], rep(1/length(mu[i,]), length(mu[i,]))) - 0.5)
      }
      estimate_diff = append(estimate_diff, quant_diff / length(truth))
    }
  }
  overlaps = check_overlap(cs)
  return(list(total=total, valid=valid, size=size, purity=purity, top=top_hit,
              overlap_var = overlaps$snp, overlap_cs = overlaps$cs,
              n_signal=length(beta_idx),
              included_signal = sum(beta_idx %in% unlist(cs)),
              false_pos_cond_discoveries = false_positive_condition_cs,
              true_cond_discoveries = true_positive_condition_cs,
              false_neg_cond_discoveries = false_neg_condition_cs,
              avg_diff_eff_size_percentile = mean(estimate_diff),
              converged = (m$convergence$delta >=0 && m$convergence$converged)))
}