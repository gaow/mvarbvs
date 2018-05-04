# workhorse(s) for finemapping-m

# Module input
# ============
# $data: full data; or
# $sumstats: summary statistics; or / and
# $ld: LD information

# Module output
# =============
# $fitted: for diagnostics
# $posterior: for inference

init_mnm: init_mnm.R
  # mashr comes from `dev` branch on github
  @CONF: R_libs = mashr
  data: $data
  reg: $sumstats
  # FIXME: these quantities are to be computed seperately and globally using mashr procedure
  # See http://stephenslab.github.io/gtex-eqtls/analysis/20171002_MASH_V8.html
  Sigma: empirical
  (U, grid, p): (auto, (0.9,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.02), auto)
  $data: data
  $model: model

fit_mnm: regression.R + fit_mnm.R
  @CONF: R_libs = mashr
  maxL: 5
  maxI: 10
  data: $data
  model: $model
  $fitted: fitted_track
  $posterior: posterior

fit_susie: fit_susie.R
  # Prior variance of nonzero effects.
  @CONF: R_libs = susieR@stephenslab/susieR
  maxL: 5
  maxI: 50
  data: $data
  $posterior: posterior
  $fitted: fitted

fit_varbvs(fit_susie): setup_varbvs.R + fit_varbvs.R
  @CONF: R_libs = varbvs@pcarbo/varbvs/varbvs-R
  sa: 1

fit_finemap: fit_finemap.R + \
             R(posterior = finemapM(
                         data$X,
                         data$Y,
                         sumstats[1,,]/sumstats[2,,],
                         ld,
                         sa, k,
                         prefix=cache))
  sumstats: $sumstats
  ld: $ld
  k: -9, (0,0,0,1)
  sa: 0.4
  cache: file(FM)
  $posterior: posterior

fit_dap: fit_dap.py
  data: $data
  $posterior: posterior

fit_dap_ss: fit_dap.py
  sumstats: $sumstats
  $posterior: posterior