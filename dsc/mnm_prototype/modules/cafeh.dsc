cafeh_prepare: cafeh.R + R(cfg = run_cafeh(sumstats$bhat, sumstats$sbhat, ld, suffstats$N, prefix))
  sumstats: $sumstats
  suffstats: $suffstats
  ld: $ld
  prefix: file(cafeh)
  $bhat_file: cfg$bhat
  $shat_file: cfg$shat
  $ld_file: cfg$LD
  $n_file: cfg$n
  
cafeh: cafeh.py
  LD: $ld_file
  bhat: $bhat_file
  shat: $shat_file
  n_file: $n_file
  m_name: file(cafeh)
  $cafeh_pip: pip
  $cafeh_trait_pip: trait_pip
  $cafeh_cs: credible_sets
  $cafeh_purity: purity
  