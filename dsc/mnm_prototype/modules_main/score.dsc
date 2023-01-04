mvsusie_scores: misc.R + mvsusie_scores.R + R(sc = mvsusie_scores($(result), $(meta)$true_coef, 0.05))
    @CONF: R_libs = KScorrect
    $n_causal: sc$n_signal
    $n_condition_causal: sc$n_condition_signal
    $total: sc$total
    $valid: sc$valid
    $size: sc$size
    $purity: sc$purity
    $top: sc$top
    $overlap_cs: sc$overlap_cs
    $overlap_var: sc$overlap_var
    $included_causal: sc$included_signal
    $total_cond_discoveries: sc$total_cond_discoveries
    $false_pos_cond_discoveries: sc$false_pos_cond_discoveries
    $false_neg_cond_discoveries: sc$false_neg_cond_discoveries
    $true_cond_discoveries: sc$true_cond_discoveries
    $size_cond_cs: sc$size_cond_cs
    $purity_cond_cs: sc$purity_cond_cs
    $avg_diff_eff_size_percentile: sc$avg_diff_eff_size_percentile
    $converged: sc$converged
    $cs_correlation: sc$cs_correlation 

mvsusie_cs_scores: misc.R + mvsusie_cs_scores.R + R(sc = mvsusie_cs_scores_multiple($(result), $(sets), $(meta)$true_coef, 0.05))
    @CONF: R_libs = KScorrect
    $n_causal: sc$n_signal
    $n_condition_causal: sc$n_condition_signal
    $total: sc$total
    $valid: sc$valid
    $size: sc$size
    $purity: sc$purity
    $included_causal: sc$included_signal
    $total_cond_discoveries: sc$total_cond_discoveries
    $false_pos_cond_discoveries: sc$false_pos_cond_discoveries
    $false_neg_cond_discoveries: sc$false_neg_cond_discoveries
    $true_cond_discoveries: sc$true_cond_discoveries
    $size_cond_cs: sc$size_cond_cs
    $purity_cond_cs: sc$purity_cond_cs
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
    $cs_correlation: sc$cs_correlation 
    $pip: sc$pip
    $n_causal: sc$expected

susie_cs_scores: susie_cs_scores.R + R(sc = susie_cs_scores_multiple($(fitted), $(sets), $(meta)$true_coef))
    $total: sc$total
    $valid: sc$valid
    $size: sc$size
    $purity: sc$purity
    $objective: sc$objective
    $converged: sc$converged
    $n_causal: sc$expected

flashfm_scores: misc.R + flashfm_scores.R + R(sc = flashfm_scores($(result), $(meta)$true_coef, $(ld)))
    $n_causal: sc$n_signal
    $n_condition_causal: sc$n_condition_signal
    $total_cond_discoveries: sc$total_cond_discoveries
    $false_pos_cond_discoveries: sc$false_pos_cond_discoveries
    $false_neg_cond_discoveries: sc$false_neg_cond_discoveries
    $true_cond_discoveries: sc$true_cond_discoveries
    $size_cond_cs: sc$size_cond_cs
    $purity_cond_cs: sc$purity_cond_cs
    $pip: sc$condition_pip

cafeh_scores: cafeh_scores.R + R(sc = cafeh_scores($(cafeh_cs), $(cafeh_purity), $(cafeh_pip), $(cafeh_trait_pip), $(meta)$true_coef))
    $n_causal: sc$n_signal
    $n_condition_causal: sc$n_condition_signal
    $total: sc$total
    $valid: sc$valid
    $size: sc$size
    $purity: sc$purity
    $pip: sc$pip
    $trait_pip: sc$trait_pip

cafeh_cs_scores: cafeh_cs_scores.R + R(sc = cafeh_cs_scores_multiple($(cs), $(purity), $(meta)$true_coef))
    $n_causal: sc$n_signal
    $n_condition_causal: sc$n_condition_signal
    $total: sc$total
    $valid: sc$valid
    $size: sc$size
    $purity: sc$purity
