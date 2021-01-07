mvsusie_scores: mvsusie_scores.R + R(sc = mvsusie_scores($(result), $(meta)$true_coef, 0.05))
    @CONF: R_libs = KScorrect
    $n_causal: sc$n_signal
    $total: sc$total
    $valid: sc$valid
    $size: median(sc$size)
    $purity: median(sc$purity)
    $top: sc$top
    $overlap_cs: sc$overlap_cs
    $overlap_var: sc$overlap_var
    $included_causal: sc$included_signal
    $false_pos_cond_discoveries: sc$false_pos_cond_discoveries
    $false_neg_cond_discoveries: sc$false_neg_cond_discoveries
    $true_cond_discoveries: sc$true_cond_discoveries
    $avg_diff_eff_size_percentile: sc$avg_diff_eff_size_percentile
    $converged: sc$converged
    
susie_scores: susie_scores.R + R(sc = susie_scores_multiple($(fitted), $(meta)$true_coef))
    $total: sc$total
    $valid: sc$valid
    $size: sc$size
    $purity: sc$purity
    $avgr2: sc$avgr2
    $top: sc$top
    $objective: sc$objective
    $converged: sc$converged
    $pip: sc$pip
